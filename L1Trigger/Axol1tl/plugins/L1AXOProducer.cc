// FWCore includes
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// File writing includes
#include "TTree.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// L1T includes
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

// ROOT
#include <TTree.h>

// hls & hls4ml includes
#include "ap_fixed.h"
#include "hls4ml/emulator.h"

#include <iostream>

class L1AXOProducer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit L1AXOProducer(const edm::ParameterSet& cfg);
  ~L1AXOProducer() override = default;

private:
  virtual void beginJob() override;
  void analyze(edm::Event const &, edm::EventSetup const &) override;
  virtual void endJob() override;

  edm::Service<TFileService> fs_;
  TTree *tree_;
  float *anomaly_score;

  edm::EDGetToken muToken;
  edm::EDGetToken egToken;
  edm::EDGetToken tauToken;
  edm::EDGetToken jetToken;
  edm::EDGetToken sumToken;

  unsigned nMu;
  unsigned nEG;
  unsigned nTau;
  unsigned nJet;
  unsigned nNNIn;

  typedef ap_fixed<16,6,AP_RND_CONV,AP_SAT> scale_t;
  typedef ap_fixed<16,6,AP_RND_CONV,AP_SAT> bias_t;
  // hls4ml emulator model path
  std::string model_so_path;
  std::vector<scale_t> scale;
  std::vector<bias_t> bias;

};

L1AXOProducer::L1AXOProducer(const edm::ParameterSet& cfg){
  // consume
  muToken = consumes<l1t::MuonBxCollection>(cfg.getParameter<edm::InputTag>("muToken"));
  egToken = consumes<l1t::EGammaBxCollection>(cfg.getParameter<edm::InputTag>("egToken"));
  tauToken = consumes<l1t::TauBxCollection>(cfg.getParameter<edm::InputTag>("tauToken"));
  jetToken = consumes<l1t::JetBxCollection>(cfg.getParameter<edm::InputTag>("jetToken"));
  sumToken = consumes<l1t::EtSumBxCollection>(cfg.getParameter<edm::InputTag>("etSumToken"));
  nMu = cfg.getParameter<unsigned>("nMu");
  nEG = cfg.getParameter<unsigned>("nEg");
  nTau = cfg.getParameter<unsigned>("nTau");
  nJet = cfg.getParameter<unsigned>("nJet");
  // total number of inputs to NN
  nNNIn = 3 * (1 + nMu + nEG + nTau + nJet);

  // store the path to the .so file
  model_so_path = cfg.getParameter<std::string>("model_so_path");

  anomaly_score = new float;
  usesResource(TFileService::kSharedResource);
  tree_ = fs_->make<TTree>("L1AXOTree", "L1AXOTree");
  tree_->Branch("anomaly_score", anomaly_score, "f");
  
}

void L1AXOProducer::analyze(edm::Event const &iEvent, edm::EventSetup const &iSetup) {
  using namespace edm;
  // get input collections
  // BXVector: first index is BX, second index is object
  edm::Handle<BXVector<l1t::Muon>> muons;
  edm::Handle<BXVector<l1t::EGamma>> egammas;
  edm::Handle<BXVector<l1t::Tau>> taus;
  edm::Handle<BXVector<l1t::Jet>> jets;
  edm::Handle<BXVector<l1t::EtSum>> sums;
  iEvent.getByToken(muToken, muons);
  iEvent.getByToken(egToken, egammas);
  iEvent.getByToken(tauToken, taus);
  iEvent.getByToken(jetToken, jets);
  iEvent.getByToken(sumToken, sums);

  // The unscaled inputs are hwInts
  // ap_fixed<14,13> is wide enough for all the ET, pT, eta, phi
  ap_fixed<18,13>* X = new ap_fixed<18,13>[nNNIn];
  // initialize to zeros
  for(unsigned i = 0; i < nNNIn; i++){
    X[i] = 0;
  }

  // fill the inputs
  unsigned ix = 0;
  // sums first, just find the MET
  for(unsigned i = 0; i < sums->size(0); i++){
    if(sums->at(0, i).getType() == l1t::EtSum::EtSumType::kMissingEt){
      X[ix++] = ((float)sums->at(0,i).hwPt()) / 2;
      X[ix++] = 0; // eta
      X[ix++] = sums->at(0,i).hwPhi();
    }
  }
  // egammas next
  for(unsigned i = 0; i < std::min(nEG, egammas->size(0)); i++){
    X[ix++] = ((float)egammas->at(0, i).hwPt()) / 2;
    X[ix++] = egammas->at(0, i).hwEta();
    X[ix++] = egammas->at(0, i).hwPhi();
  }
  // muons next
  for(unsigned i = 0; i < std::min(nMu, muons->size(0)); i++){
    X[ix++] = ((float)muons->at(0, i).hwPt()) / 2;
    X[ix++] = muons->at(0, i).hwEta();
    X[ix++] = muons->at(0, i).hwPhi();
  }
  // jets next
  for(unsigned i = 0; i < std::min(nJet, jets->size(0)); i++){
    X[ix++] = ((float)jets->at(0, i).hwPt()) / 2;
    X[ix++] = jets->at(0, i).hwEta();
    X[ix++] = jets->at(0, i).hwPhi();
  }
  // taus next
  for(unsigned i = 0; i < std::min(nTau, taus->size(0)); i++){
    X[ix++] = ((float)taus->at(0, i).hwPt()) / 2;
    X[ix++] = taus->at(0, i).hwEta();
    X[ix++] = taus->at(0, i).hwPhi();
  }

  // load the NN emulator object
  hls4mlEmulator::ModelLoader loader(model_so_path);
  std::shared_ptr<hls4mlEmulator::Model> model = loader.load_model();

  //ap_ufixed<18, 14> y; // output object
  //std::array<ap_fixed<10, 7, AP_RND_CONV, AP_SAT>, 8>
  std::pair<std::array<ap_fixed<10, 7, AP_RND_CONV, AP_SAT>, 8>,
            ap_ufixed<18, 14>> y;

  // run the actual inference
  model->prepare_input(X);
  model->predict();
  model->read_result(&y);
  std::cout << "y = " << (y.second).to_float() << std::endl;

  // write the result to the output
  // note cast from the ap_fixed emulated type to float for convenience
  *anomaly_score = (y.second).to_float();
  tree_->Fill();
}

void L1AXOProducer::beginJob(){
}

void L1AXOProducer::endJob(){
}

DEFINE_FWK_MODULE(L1AXOProducer);

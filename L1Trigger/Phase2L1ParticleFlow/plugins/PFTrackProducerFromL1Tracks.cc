#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Phase2L1ParticleFlow/interface/PFTrack.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TPFUtils.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/ParametricResolution.h"


namespace l1tpf {
    class PFTrackProducerFromL1Tracks : public edm::stream::EDProducer<> {
        public:
            explicit PFTrackProducerFromL1Tracks(const edm::ParameterSet&) ;
            ~PFTrackProducerFromL1Tracks() {}

        private:
            edm::EDGetTokenT<std::vector<l1t::PFTrack::L1TTTrackType>> TrackTag_;
            int nParam_;
            float fBz_;
            l1tpf::ParametricResolution resolCalo_, resolTrk_;

            virtual void produce(edm::Event&, const edm::EventSetup&) override;

            virtual void beginRun(edm::Run const &, edm::EventSetup const & iSetup) override {
                edm::ESHandle<MagneticField> magneticField;
                iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
                fBz_ = magneticField->inTesla(GlobalPoint(0,0,0)).z();
            }

    }; // class
} // namespace

l1tpf::PFTrackProducerFromL1Tracks::PFTrackProducerFromL1Tracks(const edm::ParameterSet & iConfig) :
    TrackTag_(consumes<std::vector<l1t::PFTrack::L1TTTrackType>>(iConfig.getParameter<edm::InputTag>("L1TrackTag"))),
    nParam_(iConfig.getParameter<unsigned int>("nParam")),
    resolCalo_(iConfig.getParameter<edm::ParameterSet>("resolCalo")),
    resolTrk_(iConfig.getParameter<edm::ParameterSet>("resolTrack"))
{
    produces<l1t::PFTrackCollection>();
}


void 
l1tpf::PFTrackProducerFromL1Tracks::produce(edm::Event & iEvent, const edm::EventSetup &) 
{
  std::unique_ptr<l1t::PFTrackCollection> out(new l1t::PFTrackCollection());
  
  // https://github.com/skinnari/cmssw/blob/80c19f1b721325c3a02ee0482f72fb974a4c3bf7/L1Trigger/TrackFindingTracklet/test/L1TrackNtupleMaker.cc
  edm::Handle<std::vector<l1t::PFTrack::L1TTTrackType>> htracks;
  iEvent.getByToken(TrackTag_, htracks);
  const auto & tracks = *htracks;

  for (unsigned int i = 0, n = tracks.size(); i < n; ++i) {
      const auto & tk = tracks[i]; 

      float pt   = tk.getMomentum(nParam_).perp();
      float eta  = tk.getMomentum(nParam_).eta();
      float phi  = tk.getMomentum(nParam_).phi();
      float z0   = tk.getPOCA(nParam_).z(); //cm
      int charge = tk.getRInv() > 0 ? +1 : -1;

      reco::Candidate::PolarLorentzVector p4p(pt, eta, phi, 0.137); // pion mass
      reco::Particle::LorentzVector p4(p4p.X(), p4p.Y(), p4p.Z(), p4p.E());
      reco::Particle::Point vtx(0.,0.,z0);

      auto caloetaphi = l1tpf::propagateToCalo(p4, math::XYZTLorentzVector(0.,0.,z0,0.), charge, fBz_);

      float trkErr = resolTrk_(pt, std::abs(eta));
      float caloErr = resolCalo_(pt, std::abs(eta));
      int quality = 1;
      out->emplace_back(charge, p4, vtx, 
                        l1t::PFTrack::TrackRef(htracks,i), 
                        nParam_,
                        caloetaphi.first, caloetaphi.second,
                        trkErr, caloErr, quality);
  }
  iEvent.put(std::move(out));
}
using l1tpf::PFTrackProducerFromL1Tracks;
DEFINE_FWK_MODULE(PFTrackProducerFromL1Tracks);

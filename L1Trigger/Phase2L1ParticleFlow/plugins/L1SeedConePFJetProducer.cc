
#include <vector>
#include <numeric>

////////////////////
// FRAMEWORK HEADERS
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/PFJet.h"
#include "DataFormats/Math/interface/deltaR.h"

// bitwise emulation headers
#include "L1Trigger/Phase2L1ParticleFlow/interface/L1SeedConePFJetEmulator.h"

class L1SeedConePFJetProducer : public edm::global::EDProducer<> {
public:
  explicit L1SeedConePFJetProducer(const edm::ParameterSet&);
  ~L1SeedConePFJetProducer() override;

private:
  /// ///////////////// ///
  /// MANDATORY METHODS ///
  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;
  /// ///////////////// ///

  float _coneSize;
  unsigned _nJets;
  bool _HW;
  bool _debug;
  edm::EDGetTokenT<std::vector<l1t::PFCandidate>> _l1PFToken;

  std::vector<l1t::PFJet> processEvent_SW(std::vector<edm::Ptr<l1t::PFCandidate>>& parts) const;
  std::vector<l1t::PFJet> processEvent_HW(std::vector<edm::Ptr<l1t::PFCandidate>>& parts) const;

  l1t::PFJet makeJet_SW(const std::vector<edm::Ptr<l1t::PFCandidate>>& parts) const;

  static std::vector<L1SCJetEmu::Particle> convertEDMToHW(std::vector<edm::Ptr<l1t::PFCandidate>>& edmParticles);
  static std::vector<l1t::PFJet> convertHWToEDM(std::vector<L1SCJetEmu::Jet> hwJets);
};

L1SeedConePFJetProducer::L1SeedConePFJetProducer(const edm::ParameterSet& cfg)
    : _coneSize(cfg.getParameter<double>("coneSize")),
      _nJets(cfg.getParameter<unsigned>("nJets")),
      _HW(cfg.getParameter<bool>("HW")),
      _debug(cfg.getParameter<bool>("debug")),
      _l1PFToken(consumes<std::vector<l1t::PFCandidate>>(cfg.getParameter<edm::InputTag>("L1PFObjects"))) {
  produces<l1t::PFJetCollection>();
}

void L1SeedConePFJetProducer::produce(edm::StreamID /*unused*/,
                                      edm::Event& iEvent,
                                      const edm::EventSetup& iSetup) const {
  std::unique_ptr<l1t::PFJetCollection> newPFJetCollection(new l1t::PFJetCollection);

  edm::Handle<l1t::PFCandidateCollection> l1PFCandidates;
  iEvent.getByToken(_l1PFToken, l1PFCandidates);

  std::vector<edm::Ptr<l1t::PFCandidate>> particles;
  for (unsigned i = 0; i < (*l1PFCandidates).size(); i++) {
    particles.push_back(edm::Ptr<l1t::PFCandidate>(l1PFCandidates, i));
  }

  std::vector<l1t::PFJet> jets;
  if (_HW) {
    jets = processEvent_HW(particles);
  } else {
    jets = processEvent_SW(particles);
  }

  std::sort(jets.begin(), jets.end(), [](l1t::PFJet i, l1t::PFJet j) { return (i.pt() > j.pt()); });
  newPFJetCollection->swap(jets);
  iEvent.put(std::move(newPFJetCollection));
}

/////////////
// DESTRUCTOR
L1SeedConePFJetProducer::~L1SeedConePFJetProducer() {}

l1t::PFJet L1SeedConePFJetProducer::makeJet_SW(const std::vector<edm::Ptr<l1t::PFCandidate>>& parts) const {
  l1t::PFCandidate seed = *parts.at(0);

  auto sumpt = [](float a, const edm::Ptr<l1t::PFCandidate>& b) { return a + b->pt(); };

  // Sum the pt
  float pt = std::accumulate(parts.begin(), parts.end(), 0., sumpt);

  // pt weighted d eta
  std::vector<float> pt_deta;
  pt_deta.resize(parts.size());
  std::transform(parts.begin(), parts.end(), pt_deta.begin(), [&seed, &pt](const edm::Ptr<l1t::PFCandidate>& part) {
    return (part->pt() / pt) * (part->eta() - seed.eta());
  });
  // Accumulate the pt weighted etas. Init to the seed eta, start accumulating at begin()+1 to skip seed
  float eta = std::accumulate(pt_deta.begin() + 1, pt_deta.end(), seed.eta());

  // pt weighted d phi
  std::vector<float> pt_dphi;
  pt_dphi.resize(parts.size());
  std::transform(parts.begin(), parts.end(), pt_dphi.begin(), [&seed, &pt](const edm::Ptr<l1t::PFCandidate>& part) {
    return (part->pt() / pt) * reco::deltaPhi(part->phi(), seed.phi());
  });
  // Accumulate the pt weighted phis. Init to the seed phi, start accumulating at begin()+1 to skip seed
  float phi = std::accumulate(pt_dphi.begin() + 1, pt_dphi.end(), seed.phi());

  l1t::PFJet jet(pt, eta, phi);
  for (auto it = parts.begin(); it != parts.end(); it++) {
    jet.addConstituent(*it);
  }

  return jet;
}

std::vector<l1t::PFJet> L1SeedConePFJetProducer::processEvent_SW(std::vector<edm::Ptr<l1t::PFCandidate>>& work) const {
  // The floating point algorithm simulation
  std::sort(work.begin(), work.end(), [](edm::Ptr<l1t::PFCandidate> i, edm::Ptr<l1t::PFCandidate> j) {
    return (i->pt() > j->pt());
  });
  std::vector<l1t::PFJet> jets;
  jets.reserve(_nJets);
  while (!work.empty() && jets.size() < _nJets) {
    // Take the first (highest pt) candidate as a seed
    edm::Ptr<l1t::PFCandidate> seed = work.at(0);
    // Get the particles within a _coneSize of the seed
    std::vector<edm::Ptr<l1t::PFCandidate>> particlesInCone;
    std::copy_if(
        work.begin(), work.end(), std::back_inserter(particlesInCone), [&](const edm::Ptr<l1t::PFCandidate>& part) {
          return reco::deltaR<l1t::PFCandidate, l1t::PFCandidate>(*seed, *part) <= _coneSize;
        });
    jets.push_back(makeJet_SW(particlesInCone));
    // remove the clustered particles
    work.erase(std::remove_if(work.begin(),
                              work.end(),
                              [&](const edm::Ptr<l1t::PFCandidate>& part) {
                                return reco::deltaR<l1t::PFCandidate, l1t::PFCandidate>(*seed, *part) <= _coneSize;
                              }),
               work.end());
  }
  return jets;
}

std::vector<l1t::PFJet> L1SeedConePFJetProducer::processEvent_HW(std::vector<edm::Ptr<l1t::PFCandidate>>& work) const {
  // The fixed point emulator
  // Convert the EDM format to the hardware format, and call the standalone emulator
  using namespace L1SCJetEmu;
  Config config(_debug, _coneSize, _nJets);
  std::vector<Particle> particles = convertEDMToHW(work);
  std::vector<Jet> jets = emulateEvent(particles, config);
  return convertHWToEDM(jets);
}

std::vector<L1SCJetEmu::Particle> L1SeedConePFJetProducer::convertEDMToHW(std::vector<edm::Ptr<l1t::PFCandidate>>& edmParticles){
  using namespace L1SCJetEmu;
  std::vector<Particle> hwParticles;
  std::for_each(edmParticles.begin(), edmParticles.end(), [&](edm::Ptr<l1t::PFCandidate>& edmParticle){
    hwParticles.push_back(L1SCJetEmu::Particle(edmParticle->pt(),
                                               edmParticle->eta() * etaphi_base,
                                               edmParticle->phi() * etaphi_base));
  });
  return hwParticles;
}

std::vector<l1t::PFJet> L1SeedConePFJetProducer::convertHWToEDM(std::vector<L1SCJetEmu::Jet> hwJets){
  using namespace L1SCJetEmu;
  std::vector<l1t::PFJet> edmJets;
  std::for_each(hwJets.begin(), hwJets.end(), [&](Jet jet){
    edmJets.push_back(l1t::PFJet(jet.hwPt,
                                 float(jet.hwEta) / etaphi_base,
                                 float(jet.hwPhi) / etaphi_base));
  });
  return edmJets;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1SeedConePFJetProducer);

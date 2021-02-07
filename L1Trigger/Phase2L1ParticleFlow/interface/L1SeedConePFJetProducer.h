#ifndef L1Trigger_Phase2L1ParticleFlow_L1SeedConePFJetProducer_h
#define L1Trigger_Phase2L1ParticleFlow_L1SeedConePFJetProducer_h

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
  l1t::PFJet makeJet_HW(const std::vector<l1t::PFCandidate>& parts) const;
};

/////////////
// DESTRUCTOR
L1SeedConePFJetProducer::~L1SeedConePFJetProducer() {}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1SeedConePFJetProducer);

#endif

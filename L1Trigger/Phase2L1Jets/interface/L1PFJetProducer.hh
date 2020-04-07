#ifndef L1PFTAU_PRDC_H
#define L1PFTAU_PRDC_H

#include <vector>

////////////////////
// FRAMEWORK HEADERS
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFJet.h"

using namespace l1t;


class L1PFJetProducer : public edm::EDProducer {
   public:
  explicit L1PFJetProducer(const edm::ParameterSet&);

  ~L1PFJetProducer();

   private:

  /// ///////////////// ///
  /// MANDATORY METHODS ///
  virtual void beginRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void endRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void produce( edm::Event& iEvent, const edm::EventSetup& iSetup ) override;
  /// ///////////////// ///

  float _coneSize;
  int _nJets;
  edm::EDGetTokenT<vector<l1t::PFCandidate>>_l1PFToken;

};


#endif

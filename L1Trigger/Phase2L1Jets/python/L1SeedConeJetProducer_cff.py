import FWCore.ParameterSet.Config as cms

L1SeedConeJetProducer = cms.EDProducer("L1SeedConeJetProducer",
                           L1PFObjects = cms.InputTag("L1PFProducer","l1pfCandidates"),
                           nJets       = cms.uint32(10),
                           coneSize    = cms.double(0.4)
                         )

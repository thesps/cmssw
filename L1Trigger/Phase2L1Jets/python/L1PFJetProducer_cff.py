import FWCore.ParameterSet.Config as cms

L1PFJetProducer = cms.EDProducer("L1PFJetProducer",
                           L1PFObjects = cms.InputTag("L1PFProducer","l1pfCandidates"),
                           nJets       = cms.int32(10),
                           coneSize    = cms.double(0.4)
                         )

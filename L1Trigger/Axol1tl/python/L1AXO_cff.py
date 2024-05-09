import FWCore.ParameterSet.Config as cms

L1AXOProducer = cms.EDAnalyzer("L1AXOProducer",
    muToken    = cms.InputTag("gtStage2Digis:Muon"),
    egToken    = cms.InputTag("gtStage2Digis:EGamma"),
    tauToken   = cms.InputTag("gtStage2Digis:Tau"),
    jetToken   = cms.InputTag("gtStage2Digis:Jet"),
    etSumToken = cms.InputTag("gtStage2Digis:EtSum"),
    nMu = cms.uint32(4),
    nEg = cms.uint32(4),
    nTau = cms.uint32(0),
    nJet = cms.uint32(10),
    model_so_path = cms.string("GTADModel_v3"),
)

def L1AXO(process):
  process.L1AXOProducer = L1AXOProducer.clone()
  process.L1AXO = cms.Path(process.L1AXOProducer)
  process.schedule.append(process.L1AXO)
  return process
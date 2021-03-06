import FWCore.ParameterSet.Config as cms
import sys
#
# Legacy L1 Muon modules still running in 2016 trigger:
#

#  - DT TP emulator
from L1Trigger.DTTrigger.dtTriggerPrimitiveDigis_cfi import *
import L1Trigger.DTTrigger.dtTriggerPrimitiveDigis_cfi
simDtTriggerPrimitiveDigis = L1Trigger.DTTrigger.dtTriggerPrimitiveDigis_cfi.dtTriggerPrimitiveDigis.clone(
    digiTag = 'simMuonDTDigis'
)
#simDtTriggerPrimitiveDigis.debug = cms.untracked.bool(True)

# - CSC TP emulator
from L1Trigger.CSCCommonTrigger.CSCCommonTrigger_cfi import *
import L1Trigger.CSCTriggerPrimitives.cscTriggerPrimitiveDigis_cfi
simCscTriggerPrimitiveDigis = L1Trigger.CSCTriggerPrimitives.cscTriggerPrimitiveDigis_cfi.cscTriggerPrimitiveDigis.clone(
    CSCComparatorDigiProducer = 'simMuonCSCDigis:MuonCSCComparatorDigi',
    CSCWireDigiProducer       = 'simMuonCSCDigis:MuonCSCWireDigi'
)

SimL1TMuonCommon = cms.Sequence(simDtTriggerPrimitiveDigis + simCscTriggerPrimitiveDigis)

#
# Legacy Trigger:
#
#
# - CSC Track Finder emulator
#
import L1Trigger.CSCTrackFinder.csctfTrackDigis_cfi
simCsctfTrackDigis = L1Trigger.CSCTrackFinder.csctfTrackDigis_cfi.csctfTrackDigis.clone(
    SectorReceiverInput = 'simCscTriggerPrimitiveDigis:MPCSORTED',
    DTproducer = 'simDtTriggerPrimitiveDigis'
)
import L1Trigger.CSCTrackFinder.csctfDigis_cfi
simCsctfDigis = L1Trigger.CSCTrackFinder.csctfDigis_cfi.csctfDigis.clone(
    CSCTrackProducer = 'simCsctfTrackDigis'
)
#
# - DT Track Finder emulator
# 
import L1Trigger.DTTrackFinder.dttfDigis_cfi
simDttfDigis = L1Trigger.DTTrackFinder.dttfDigis_cfi.dttfDigis.clone(
    DTDigi_Source  = 'simDtTriggerPrimitiveDigis',
    CSCStub_Source = 'simCsctfTrackDigis'
)
#
# - RPC PAC Trigger emulator
#
from L1Trigger.RPCTrigger.rpcTriggerDigis_cff import *
simRpcTriggerDigis = L1Trigger.RPCTrigger.rpcTriggerDigis_cff.rpcTriggerDigis.clone(
    label = 'simMuonRPCDigis'
)
#
# - Global Muon Trigger emulator
#
import L1Trigger.GlobalMuonTrigger.gmtDigis_cfi
simGmtDigis = L1Trigger.GlobalMuonTrigger.gmtDigis_cfi.gmtDigis.clone(
    DTCandidates   = 'simDttfDigis:DT',
    CSCCandidates  = 'simCsctfDigis:CSC',
    RPCbCandidates = 'simRpcTriggerDigis:RPCb',
    RPCfCandidates = 'simRpcTriggerDigis:RPCf',
#   Note: GMT requires input from calorimeter emulators, namely MipIsoData from GCT
    MipIsoData     = 'simRctDigis'
)
#
#
SimL1TMuon = cms.Sequence(SimL1TMuonCommon + simCsctfTrackDigis + simCsctfDigis + simDttfDigis + simRpcTriggerDigis + simGmtDigis)

#
# Stage-2 Trigger
#
from L1Trigger.L1TTwinMux.simTwinMuxDigis_cfi import *
from L1Trigger.L1TMuonBarrel.simBmtfDigis_cfi import *
from L1Trigger.L1TMuonEndCap.simEmtfDigis_cfi import *
from L1Trigger.L1TMuonOverlap.simOmtfDigis_cfi import *
from L1Trigger.L1TMuon.simGmtCaloSumDigis_cfi import *
from L1Trigger.L1TMuon.simGmtStage2Digis_cfi import *
from L1Trigger.L1TMuonBarrel.simKBmtfStubs_cfi import *
from L1Trigger.L1TMuonBarrel.simKBmtfDigis_cfi import *
from Configuration.Eras.Modifier_stage2L1Trigger_cff import stage2L1Trigger
#
#
stage2L1Trigger.toReplaceWith(SimL1TMuon, cms.Sequence(SimL1TMuonCommon + simTwinMuxDigis + simBmtfDigis + simKBmtfStubs + simKBmtfDigis + simEmtfDigis + simOmtfDigis + simGmtCaloSumDigis + simGmtStage2Digis))

from L1Trigger.ME0Trigger.me0TriggerPseudoDigis_cff import *
from RecoLocalMuon.RPCRecHit.rpcRecHits_cfi import rpcRecHits
rpcRecHits.rpcDigiLabel = 'simMuonRPCDigis'
_phase2_SimL1TMuon = SimL1TMuon.copy()
_phase2_SimL1TMuon.replace(simEmtfDigis, me0TriggerPseudoDigiSequence + rpcRecHits + simEmtfDigis)

from Configuration.Eras.Modifier_phase2_muon_cff import phase2_muon
(stage2L1Trigger & phase2_muon).toReplaceWith( SimL1TMuon, _phase2_SimL1TMuon )

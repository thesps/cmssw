import FWCore.ParameterSet.Config as cms

<<<<<<< HEAD
from L1Trigger.Phase2L1ParticleFlow.l1tSeedConePFJetProducer_cfi import l1tSeedConePFJetProducer, l1tSeedConePFJetEmulatorProducer
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer, l1tDeregionizerProducerExtended as l1tLayer2DeregionizerExtended
l1tSCPFL1PF            = l1tSeedConePFJetProducer.clone(L1PFObjects = 'l1tLayer1:PF')
l1tSCPFL1Puppi         = l1tSeedConePFJetProducer.clone()
l1tSCPFL1PuppiEmulator = l1tSeedConePFJetEmulatorProducer.clone(L1PFObjects = 'l1tLayer2Deregionizer:Puppi')
l1tSCPFL1PuppiCorrectedEmulator = l1tSeedConePFJetEmulatorProducer.clone(L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
                                                                     doCorrections = cms.bool(True),
                                                                     correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
                                                                     correctorDir = cms.string('L1PuppiSC4EmuJets'))
=======
from L1Trigger.Phase2L1ParticleFlow.L1SeedConePFJetProducer_cfi import L1SeedConePFJetProducer, L1SeedConePFJetEmulatorProducer
from L1Trigger.Phase2L1ParticleFlow.DeregionizerProducer_cfi import DeregionizerProducer as l1ctLayer2Deregionizer, DeregionizerProducerExtended as l1ctLayer2DeregionizerExtended
sc4PFL1PF            = L1SeedConePFJetProducer.clone(L1PFObjects = 'l1ctLayer1:PF')
sc4PFL1Puppi         = L1SeedConePFJetProducer.clone()
sc4PFL1PuppiEmulator = L1SeedConePFJetEmulatorProducer.clone(L1PFObject = cms.InputTag('l1ctLayer2Deregionizer', 'Puppi'))
sc4PFL1PuppiCorrectedEmulator = L1SeedConePFJetEmulatorProducer.clone(L1PFObject = cms.InputTag('l1ctLayer2Deregionizer', 'Puppi'),
                                                                      doCorrections = cms.bool(True),
                                                                      correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
                                                                      correctorDir = cms.string('L1PuppiSC4EmuJets'))
sc8PFL1PuppiCorrectedEmulator = L1SeedConePFJetEmulatorProducer.clone(L1PFObject = cms.InputTag('l1ctLayer2Deregionizer', 'Puppi'),
                                                                      coneSize = cms.double(0.8),
                                                                      doCorrections = cms.bool(True),
                                                                      correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
                                                                      correctorDir = cms.string('L1PuppiSC4EmuJets'))
>>>>>>> 5adaf41... Add SC8 jets. Rename SC jets to SC4. Add multiple outputs to jet file writer

_correctedJets = cms.EDProducer("L1TCorrectedPFJetProducer", 
    jets = cms.InputTag("_tag_"),
    correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
    correctorDir = cms.string("_dir_"),
    copyDaughters = cms.bool(False),
    emulate = cms.bool(False)
)

# Using phase2_hgcalV10 to customize the config for all 106X samples, since there's no other modifier for it
from Configuration.Eras.Modifier_phase2_hgcalV10_cff import phase2_hgcalV10
phase2_hgcalV10.toModify(_correctedJets, correctorFile = "L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs.PU200_106X.root")
from Configuration.Eras.Modifier_phase2_hgcalV11_cff import phase2_hgcalV11
phase2_hgcalV11.toModify(_correctedJets, correctorFile = "L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root")

<<<<<<< HEAD
from L1Trigger.Phase2L1ParticleFlow.l1tMHTPFProducer_cfi import l1tMHTPFProducer
l1tSCPFL1PuppiCorrectedEmulatorMHT = l1tMHTPFProducer.clone() 

L1TPFJetsTask = cms.Task(
    l1tLayer2Deregionizer, l1tSCPFL1PF, l1tSCPFL1Puppi, l1tSCPFL1PuppiEmulator, l1tSCPFL1PuppiCorrectedEmulator, l1tSCPFL1PuppiCorrectedEmulatorMHT
)

l1tSCPFL1PuppiExtended         = l1tSeedConePFJetProducer.clone(L1PFObjects = 'l1tLayer1Extended:Puppi')
l1tSCPFL1PuppiExtendedEmulator = l1tSeedConePFJetEmulatorProducer.clone(L1PFObjects = 'l1tLayer2DeregionizerExtended:Puppi')
l1tSCPFL1PuppiExtendedCorrectedEmulator = l1tSeedConePFJetEmulatorProducer.clone(L1PFObjects = 'l1tLayer2DeregionizerExtended:Puppi',
                                                                     doCorrections = cms.bool(True),
                                                                     correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
                                                                     correctorDir = cms.string('L1PuppiSC4EmuJets'))

L1TPFJetsTask = cms.Task(
    l1tLayer2Deregionizer, l1tSCPFL1PF, l1tSCPFL1Puppi, l1tSCPFL1PuppiEmulator, l1tSCPFL1PuppiCorrectedEmulator, l1tSCPFL1PuppiCorrectedEmulatorMHT
)

L1TPFJetsExtendedTask = cms.Task(
    l1tLayer2DeregionizerExtended, l1tSCPFL1PuppiExtended, l1tSCPFL1PuppiExtendedEmulator, l1tSCPFL1PuppiExtendedCorrectedEmulator
=======
from L1Trigger.Phase2L1ParticleFlow.L1MhtPfProducer_cfi import L1MhtPfProducer
sc4PFL1PuppiCorrectedEmulatorMHT = L1MhtPfProducer.clone(jets = cms.InputTag("sc4PFL1PuppiCorrectedEmulator"))
sc8PFL1PuppiCorrectedEmulatorMHT = L1MhtPfProducer.clone(jets = cms.InputTag("sc8PFL1PuppiCorrectedEmulator"))


sc4PFL1PuppiExtended = sc4PFL1Puppi.clone(L1PFObjects = 'l1ctLayer1Extended:Puppi')
sc4PFL1PuppiExtendedEmulator = sc4PFL1PuppiEmulator.clone(L1PFObjects = cms.InputTag('l1ctLayer2DeregionizerExtended', 'Puppi'))
sc4PFL1PuppiExtendedCorrectedEmulator = _correctedJets.clone(jets = 'sc4PFL1PuppiExtendedEmulator', correctorDir = 'L1PuppiSC4EmuDeregJets')

l1PFJetsTask = cms.Task(
    l1ctLayer2Deregionizer, sc4PFL1PF, sc4PFL1Puppi, sc4PFL1PuppiEmulator, sc4PFL1PuppiCorrectedEmulator, sc4PFL1PuppiCorrectedEmulatorMHT,
    sc8PFL1PuppiCorrectedEmulator, sc8PFL1PuppiCorrectedEmulatorMHT
)
l1PFJetsExtendedTask = cms.Task(
    l1ctLayer2DeregionizerExtended, sc4PFL1PuppiExtended, sc4PFL1PuppiExtendedEmulator, sc4PFL1PuppiExtendedCorrectedEmulator
>>>>>>> 5adaf41... Add SC8 jets. Rename SC jets to SC4. Add multiple outputs to jet file writer
)

L1TPFJetsEmulationTask = cms.Task(
    l1tLayer2Deregionizer, l1tSCPFL1PuppiEmulator, l1tSCPFL1PuppiCorrectedEmulator, l1tSCPFL1PuppiCorrectedEmulatorMHT
)


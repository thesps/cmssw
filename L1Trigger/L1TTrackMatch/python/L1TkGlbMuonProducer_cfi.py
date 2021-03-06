import FWCore.ParameterSet.Config as cms

L1TkGlbMuons = cms.EDProducer("L1TkGlbMuonProducer",
    L1MuonInputTag        = cms.InputTag("simGmtStage2Digis",""),
    L1TrackInputTag       = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
    ETAMIN = cms.double(0),
    ETAMAX = cms.double(5.),        # no cut
    ZMAX = cms.double( 25. ),       # in cm
    CHI2MAX = cms.double( 100. ),
    PTMINTRA = cms.double( 2. ),    # in GeV
#    DRmax = cms.double( 0.5 ),
    nStubsmin = cms.int32( 4 ),        # minimum number of stubs
#    closest = cms.bool( True ),
    correctGMTPropForTkZ = cms.bool(True),
    use5ParameterFit = cms.bool(False), #use 4-pars by defaults
    useTPMatchWindows = cms.bool(False),
)



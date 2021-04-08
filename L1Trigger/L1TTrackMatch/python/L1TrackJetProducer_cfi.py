import FWCore.ParameterSet.Config as cms

L1TrackJets = cms.EDProducer('L1TrackJetProducer',
	L1TrackInputTag= cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"),
	trk_zMax = cms.double (15.) ,    # maximum track z
	trk_ptMax = cms.double(200.),    # maximumum track pT before saturation [GeV]
	trk_ptMin = cms.double(2.0),     # minimum track pt [GeV]
   	trk_etaMax = cms.double(2.4),    # maximum track eta
	trk_chi2dofMax=cms.double(10.),	 # maximum track chi2/dof
	trk_bendChi2Max=cms.double(2.2), # maximum track bendchi2
	trk_nPSStubMin=cms.int32(-1),    # minimum PS stubs, -1 means no cut
	minTrkJetpT=cms.double(5.),      # minimum track pt to be considered for track jet
	etaBins=cms.int32(24),
	phiBins=cms.int32(27),
	zBins=cms.int32(60),
	d0_cutNStubs4=cms.double(0.15),
	d0_cutNStubs5=cms.double(0.5),
	lowpTJetMinTrackMultiplicity=cms.int32(2),
	highpTJetMinTrackMultiplicity=cms.int32(3),
	displaced=cms.bool(False), #Flag for displaced tracks
	nStubs4DisplacedChi2_Loose=cms.double(5.0), #Displaced track quality flags for loose/tight
	nStubs4Displacedbend_Loose=cms.double(1.7),
	nStubs5DisplacedChi2_Loose=cms.double(2.75),
	nStubs5Displacedbend_Loose=cms.double(3.5),
	nStubs4DisplacedChi2_Tight=cms.double(12.0),
	nStubs4Displacedbend_Tight=cms.double(1.0),
	nStubs5DisplacedChi2_Tight=cms.double(2.75),
	nStubs5Displacedbend_Tight=cms.double(3.5)
)

L1TrackJetsExtended = cms.EDProducer('L1TrackJetProducer',
	L1TrackInputTag= cms.InputTag("TTTracksFromExtendedTrackletEmulation", "Level1TTTracks"),
	trk_zMax = cms.double (15.) ,    # maximum track z
	trk_ptMax = cms.double(200.),    # maximumum track pT before saturation [GeV]
	trk_ptMin = cms.double(3.0),     # minimum track pt [GeV]
   	trk_etaMax = cms.double(2.4),    # maximum track eta
	trk_chi2dofMax=cms.double(40.),	 # maximum track chi2/dof
	trk_bendChi2Max=cms.double(2.4), # maximum track bendchi2
	trk_nPSStubMin=cms.int32(-1),    # minimum # PS stubs, -1 means no cut
	minTrkJetpT=cms.double(5.),      # minimum track pt to be considered for track jet
	etaBins=cms.int32(24),
	phiBins=cms.int32(27),
	zBins=cms.int32(10),
	d0_cutNStubs4=cms.double(0.15),
	d0_cutNStubs5=cms.double(0.5),
	lowpTJetMinTrackMultiplicity=cms.int32(2),
	highpTJetMinTrackMultiplicity=cms.int32(3),
	displaced=cms.bool(True), #Flag for displaced tracks
	nStubs4DisplacedChi2_Loose=cms.double(5.0), #Displaced track quality flags for loose/tight
	nStubs4Displacedbend_Loose=cms.double(1.7),
	nStubs5DisplacedChi2_Loose=cms.double(2.75),
	nStubs5Displacedbend_Loose=cms.double(3.5),
	nStubs4DisplacedChi2_Tight=cms.double(12.0),
	nStubs4Displacedbend_Tight=cms.double(1.0),
	nStubs5DisplacedChi2_Tight=cms.double(2.75),
	nStubs5Displacedbend_Tight=cms.double(3.5)
)


import FWCore.ParameterSet.Config as cms
rootuple = cms.EDAnalyzer(
    'NtupleMaker',
    dimuons = cms.InputTag("onia2MuMuPAT"),
    #primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    primaryVertices = cms.InputTag("offlinePrimaryVerticesWithBS"),
    offlineBeamSpot = cms.InputTag("offlineBeamSpot"),
    # conversions = cms.InputTag("oniaPhotonCandidates","conversions"), # "conversions"
    muons = cms.InputTag("oniaSelectedMuons"),
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    TriggerSummaryAOD = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    #hltName = cms.untracked.string("HLT"),
    #onia_pdgid = cms.uint32(443),
    #onia_mass_cuts = cms.vdouble(0.1,500.0),
    #isMC = cms.bool(False),
    #OnlyBest = cms.bool(False),
    #run_YY_SAF = cms.bool(True),
    #OnlyGen = cms.bool(False),
    #TriggerMatch = cms.bool(False),
    #VERBOSE = cms.untracked.bool(False),              
    #VERBOSE = cms.untracked.bool(True), 
    triggerList = cms.untracked.vstring(
        'HLT_Dimuon0_Jpsi3p5_Muon2_v',
        # 'HLT_Dimuon0_Jpsi_Muon_v',
        #'HLT_Trimuon5_3p5_2_Upsilon_Muon_v',       #2018  HLT string
        #'HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_v',  #2018 HLT string
        #'HLT_Trimuon2_Upsilon5_Muon_v'            #2017B HLT string
        #'HLT_Trimuon5_3p5_2_Upsilon_Muon_v'       #2017 C,D, E and F HLT string
        #'HLT_Dimuon0_Upsilon_Muon_v',             #2016 HLT string 
        #'HLT_Dimuon0_Phi_Barrel_v',               #2016 HLT string 
        #'HLT_Dimuon13_Upsilon_v',                 #2016 HLT string 
        #'HLT_Dimuon8_Upsilon_Barrel_v',           #2016 HLT string
    ),
    #upsilon_mass = cms.double(3.0969), #9.4603. 3.0969
    #upsilon_sigma = cms.double(0.00004),
    #upsilon_vertex_CL_cut = cms.double(0.005),
    #triggerCuts = cms.uint32(73),  # 36 for Upsilon, 73 for Jpsi
    #TRG_Match_DR = cms.double(0.3),
    #TRG_Match_DP = cms.double(1),
    #TRG_Match_DP_ratio = cms.double(0.1), 
    #best4muonCand = cms.bool(False),
    SecondSource = cms.SecSource(
        "EmbeddedRootSource",
        fileNames = cms.untracked.vstring(
            'SpyFileNameWhichNeedsToBeSet SiStripSpyEventMatcher.SpySource.fileNames'
        ),
        sequential = cms.untracked.bool(True),
    )
)


import FWCore.ParameterSet.Config as cms
rootuple = cms.EDAnalyzer(
    'NtupleMaker',
    # dimuons = cms.InputTag("onia2MuMuPAT"),
    # conversions = cms.InputTag("oniaPhotonCandidates","conversions"),# "conversions"
    # TriggerSummaryAOD = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),# MINAOD #, "", "PAT"
    # primaryVertices = cms.InputTag("offlinePrimaryVertices"),# Monte-Carlo
    offlineBeamSpot = cms.InputTag("offlineBeamSpot"),#, "", "RECO"
    muons = cms.InputTag("slimmedMuons", "", "PAT"),# MINIAOD
    # muons = cms.InputTag("oniaSelectedMuons"),# Monte-Carlo
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    triggerList = cms.untracked.vstring(
        # 'HLT_Dimuon0_Jpsi3p5_Muon2_v',
        'HLT_Dimuon0_Jpsi_Muon_v',
        #'HLT_Trimuon5_3p5_2_Upsilon_Muon_v',       #2018  HLT string
        #'HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_v',  #2018 HLT string
        #'HLT_Trimuon2_Upsilon5_Muon_v'            #2017B HLT string
        #'HLT_Trimuon5_3p5_2_Upsilon_Muon_v'       #2017 C,D, E and F HLT string
        #'HLT_Dimuon0_Upsilon_Muon_v',             #2016 HLT string 
        #'HLT_Dimuon0_Phi_Barrel_v',               #2016 HLT string 
        #'HLT_Dimuon13_Upsilon_v',                 #2016 HLT string 
        #'HLT_Dimuon8_Upsilon_Barrel_v',           #2016 HLT string
    ),
    Candidates = cms.InputTag("packedPFCandidates", "", "PAT"),
    SecondSource = cms.SecSource(
        "EmbeddedRootSource",
        fileNames = cms.untracked.vstring(
            'SpyFileNameWhichNeedsToBeSet SiStripSpyEventMatcher.SpySource.fileNames'
        ),
        sequential = cms.untracked.bool(True),
    )
)

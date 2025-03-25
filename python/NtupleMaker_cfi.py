
import FWCore.ParameterSet.Config as cms
rootuple = cms.EDAnalyzer(
    'NtupleMaker',
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),# MINAOD
    offlineBeamSpot = cms.InputTag("offlineBeamSpot"),
    muons = cms.InputTag("slimmedMuons"),# MINIAOD
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    triggerList = cms.untracked.vstring(
        'HLT_Dimuon0_Jpsi3p5_Muon2_v',
        # 'HLT_Dimuon0_Jpsi_Muon_v',
        # 'HLT_Trimuon5_3p5_2_Upsilon_Muon_v',       #2018  HLT string
        # 'HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_v',  #2018 HLT string
        # 'HLT_Trimuon2_Upsilon5_Muon_v'            #2017B HLT string
        # 'HLT_Trimuon5_3p5_2_Upsilon_Muon_v'       #2017 C,D, E and F HLT string
        # 'HLT_Dimuon0_Upsilon_Muon_v',             #2016 HLT string 
        # 'HLT_Dimuon0_Phi_Barrel_v',               #2016 HLT string 
        # 'HLT_Dimuon13_Upsilon_v',                 #2016 HLT string 
        # 'HLT_Dimuon8_Upsilon_Barrel_v',           #2016 HLT string
    ),
    SecondSource = cms.SecSource(
        "EmbeddedRootSource",
        fileNames = cms.untracked.vstring(
            'SpyFileNameWhichNeedsToBeSet SiStripSpyEventMatcher.SpySource.fileNames'
        ),
        sequential = cms.untracked.bool(True),
    )
)

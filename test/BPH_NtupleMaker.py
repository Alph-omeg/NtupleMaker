# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: skims -s SKIM:BPHSkim --dasquery=file dataset=/Charmonium/Run2016B-PromptReco-v1/AOD -n -1 --data --conditions auto:com10 --python_filename=crab-skim.py --processName=BPHSkim --no_exec

import FWCore.ParameterSet.Config as cms

process = cms.Process('NtupleMaker')

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.Skims_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
# process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # '/store/data/Run2018C/Charmonium/MINIAOD/UL2018_MiniAODv2_BParking-v1/2520000/03BA1A06-DA13-DD4B-B295-345D3D71021D.root'
        # '/store/data/Run2018B/Charmonium/MINIAOD/UL2018_MiniAODv2_BParking-v1/2520000/02D12A94-4055-504E-9ADD-816EF71D2F58.root'
        # '/store/data/Run2018A/Charmonium/MINIAOD/15Feb2022_UL2018-v1/2820000/000CB244-FD5D-AE42-B487-02E61BEED1F1.root'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('Ntuple_2018_B.root'),
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True), SkipEvent = cms.untracked.vstring('ProductNotFound'))

# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v32', '')# MINIAOD

# Schedule definition
process.load('NtupleMaker.NtupleMaker.NtupleMaker_cfi')
process.p = cms.Path(process.rootuple)
process.schedule = cms.Schedule(process.p)

# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: skims -s SKIM:BPHSkim --dasquery=file dataset=/Charmonium/Run2016B-PromptReco-v1/AOD -n -1 --data --conditions auto:com10 --python_filename=crab-skim.py --processName=BPHSkim --no_exec

#2020.10.10 Jinfeng
#Merge the Ntuple part into this file
import FWCore.ParameterSet.Config as cms

process = cms.Process('NtupleMaker')

# import of standard configurations
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
# process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.Skims_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
# '/store/data/Run2016B/Charmonium/AOD/21Feb2020_ver2_UL2016_HIPM-v1/240000/0011355B-1D48-A447-9343-0BEF32F09D9A.root'	
'file:/eos/home-c/chensh/JPsiJPsi/SKIM_tightfilter/SPS/ULPythia2016/CMSSW_10_6_20/src/NtupleMaker/NtupleMaker/test/BPH-MiniAOD-13TeV_1.root'
),
    secondaryFileNames = cms.untracked.vstring()
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('Ntuple_2016_SPS.root'),
)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

# Production Info
# process.configurationMetadata = cms.untracked.PSet(
#     annotation = cms.untracked.string('skims nevts:-1'),
#     name = cms.untracked.string('Applications'),
#     version = cms.untracked.string('$Revision: 1.19 $')
# )

# Output definition

# process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
#     dataset = cms.untracked.PSet(
#         dataTier = cms.untracked.string(''),
#         filterName = cms.untracked.string('')
#     ),
#     eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
#     fileName = cms.untracked.string('skims_SKIM.root'),
#     outputCommands = process.RECOSIMEventContent.outputCommands,
#     splitLevel = cms.untracked.int32(0)
# )

process.BPHSkimSequence = cms.Sequence()#process.oniaPATMuonsWithoutTrigger*process.oniaSelectedMuons*process.onia2MuMuPAT

# Additional output definition
# process.SKIMStreamBPHSkim = cms.OutputModule("PoolOutputModule",
#     SelectEvents = cms.untracked.PSet(
#         SelectEvents = cms.vstring('BPHSkimPath')
#     ),
#     dataset = cms.untracked.PSet(
#         dataTier = cms.untracked.string('USER'),
#         filterName = cms.untracked.string('BPHSkim')
#     ),
#     eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
#     fileName = cms.untracked.string('BPHSkim_UL.root'),
#     outputCommands = cms.untracked.vstring('keep *','drop *_MEtoEDMConverter_*_*',
#     'drop *_*Stage2Digis_*_*',
#     'drop *_ak*_*_*',
#     'drop CTPPS*_*_*_*',
#     'drop *_dedx*_*_*',
#     'drop Hcal*_*_*_*',
#     'drop *_cmsTopTagPFJetsCHS_*_*',
#     'drop recoJetedmRefToBase*_*_*_*',
#     'drop Totem*_*_*_*',
#     'drop recoPFTauDiscriminator_*_*_*',
#     'drop *_selectDigi_*_*',
#     'drop recoBasicJets_*_*_*',
#     'drop recoPFJets*_*_*_*',
#     'drop recoPFMETs*_*_*_*',
#     'drop recoPFTaus*_*_*_*',
#     'drop recPreshowerCluster*_*_*_*',
#  )
# )

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun2_asymptotic_v13', '')		#for UL data 106X_dataRun2_v32

process.oniaSelectedMuons.cut = cms.string(''
#                    'muonID(\"TMOneStationTight\")'
#                    ' && abs(innerTrack.dxy) < 0.3'
#                    ' && abs(innerTrack.dz)  < 20.'
#                    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
#                    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
#                    ' && innerTrack.quality(\"highPurity\")'
#                    ' && (abs(eta) <= 2.5 && pt > 1.8)'
)

#process.onia2MuMuPAT.muons=cms.InputTag('oniaPATMuonsWithoutTrigger')
process.onia2MuMuPAT.higherPuritySelection = cms.string("(isGlobalMuon || isTrackerMuon || (innerTrack.isNonnull && genParticleRef(0).isNonnull)) && abs(innerTrack.dxy)<4 && abs(innerTrack.dz)<35 && muonID('TrackerMuonArbitrated')")
process.onia2MuMuPAT.lowerPuritySelection = cms.string("(isGlobalMuon || isTrackerMuon || (innerTrack.isNonnull && genParticleRef(0).isNonnull)) && abs(innerTrack.dxy)<4 && abs(innerTrack.dz)<35 && muonID('TrackerMuonArbitrated')")

#process.onia2MuMuPAT.higherPuritySelection = cms.string("")
#process.onia2MuMuPAT.lowerPuritySelection = cms.string("")

# Path and EndPath definitions
# process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
# process.SKIMStreamBPHSkimOutPath = cms.EndPath(process.SKIMStreamBPHSkim)

# Schedule definition

process.load('NtupleMaker.NtupleMaker.NtupleMaker_cfi')
process.p = cms.Path(process.rootuple)
process.schedule = cms.Schedule(process.BPHSkimPath, process.p)

# process.rootuple.upsilon_mass = cms.double(3.0969)
# process.rootuple.triggerCuts = cms.uint32(73)
# process.rootuple.isMC = cms.bool(False)                 # is mc?
# process.rootuple.onia_mass_cuts = cms.vdouble(0.1,100)    # you may need to adjust this


# Schedule definition
#process.schedule = cms.Schedule(process.BPHSkimPath,process.RECOSIMoutput_step,process.SKIMStreamBPHSkimOutPath)
#from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
#associatePatAlgosToolsTask(process)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
#from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
#process = customiseEarlyDelete(process)
# End adding early deletion


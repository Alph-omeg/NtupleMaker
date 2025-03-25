from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'Charmonium_2018B_Ntuple'

config.section_('JobType')
config.JobType.psetName = './BPH_NtupleMaker.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['Ntuple_2018_B.root']

config.section_('Data')
config.Data.inputDataset = '/Charmonium/Run2018B-UL2018_MiniAODv2_BParking-v1/MINIAOD'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 16
config.Data.outputDatasetTag = '2018B_Ntuple_chensh_v2'
# NJOBS = 5000  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/chensh/Data2018/B/'#lxplus username here
# config.Data.lumiMask = 'Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt' #2016
# config.Data.lumiMask = 'Cert_294927-306462_13TeV_UL2017_Collisions17_JSON_MuonJSON.txt' #2017
config.Data.lumiMask = 'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON_MuonPhys.txt' #2018

config.section_('User')
config.JobType.allowUndistributedCMSSW = True

config.section_('Site')
config.Site.storageSite = 'T2_CN_Beijing' 

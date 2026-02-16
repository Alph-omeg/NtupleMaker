from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'SPStoJJ_2016_Ntuple'

config.section_('JobType')
config.JobType.psetName = './BPH_NtupleMaker.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['Ntuple_2016_SPS.root']

config.section_('Data')
config.Data.inputDataset = '/Pythia8_SPStoJJ_filter3p5/chensh-MC2016_MiniAOD_chensh_v1_1-c75637bb007d2ca9de9be3457fb23434/USER'

config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1000
config.Data.outputDatasetTag = 'MC2016_Ntuple_chensh_v1'
#NJOBS = 5000  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/chensh/MC2016/SKIM_tightfilter/' #lxplus username here
#config.Data.lumiMask = 'Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt' #2016
#config.Data.lumiMask = 'Cert_294927-306462_13TeV_UL2017_Collisions17_JSON_MuonJSON.txt' #2017
# config.Data.lumiMask = 'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON_MuonPhys.txt' #2018

config.section_('User')
config.JobType.allowUndistributedCMSSW = True

config.section_('Site')
config.Site.storageSite = 'T2_CN_Beijing' 

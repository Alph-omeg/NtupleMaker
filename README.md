# NtupleMaker
This program is dedicated to analyzing `Run2 MINIAOD` charmomium data and runs under `CMSSW_10_6_20`

## 1. Implementation
The program uses PAT level information to reconstruct $pp\to J/\psi+J/\psi$ event and record necessary variables.
The structure of the program is as follows:
```bash
▾ NtupleMaker_repo
│ ▾ python
│ │ __init__.py
│ └ NtupleMaker_cfi.py
│ ▾ src
│ └ NtupleMaker.cc
│ ▾ test
│ │ BPH_NtupleMaker.py
│ │ Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt
│ │ Cert_294927-306462_13TeV_UL2017_Collisions17_MuonJSON.txt
│ │ Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON_MuonPhys.txt
│ └ crab3_Ntuple.py
│ BuildFile.xml
│ Makefile
└ README.md
```
`src/NtupleMaker.cc` is the core of this program, which implements most of analysis utilities.

`python/NtupleMaker_cfi.py` specifies tokens used in C program.

`test/BPH_NtupleMaker.py` allows running on poolsource to check program validity.

`test/crab3_Ntuple.py` submit tasks over whole datasets to CRAB.

## 2. Setup
In `CMSSW_10_6_20/src` directory, create new directory as work directory.
```bash
mkdir NtupleMaker && cd NtupleMaker
```
git clone this repository via
```bash
git clone git@github.com:Alph-omeg/NtupleMaker.git
```
or
```bash
git clone https://github.com/Alph-omeg/NtupleMaker.git
```

## 3. Modify
The following requirement in `src/NtupleMaker.cc` applies for SPS process. Quote this line if dealing with DPS process.
```c++
// line 551
    if(temp_JpsiMom[i] != temp_JpsiMom[j]) continue;//only for SPS
```
Change input(PoolSource) and output file name in `test/BPH_NtupleMaker.py`.
```python
# line 31
    fileNames = cms.untracked.vstring(
        "/path/fileName.root"
    ),
# line 40
        fileName = cms.string('outputFile.root'),
```
Change CRAB configuration in `test/crab3_Ntuple.py`.
```python
# line 6
config.General.requestName = 'requestName'
# line 11
config.JobType.outputFiles = ['outputFile.root']
# line 14
config.Data.inputDataset = '/Charmonium/.../MINIAOD'
# line 19
config.Data.outputDatasetTag = 'outputDsTag'
# line 22
config.Data.outLFNDirBase = '/LFH_Path/'
```

## 4. Run
Compile the program after modification.
```bash
scram b
```
Test the program via
```bash
cmsRun BPH_NtupleMaker.py
```
Submit job to CRAB via
```bash
crab sub crab3_Ntuple.py
```

from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

options = VarParsing('python')

options.register('isMC', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('isPreECALleakage',False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Data taken during ECAL leakage"
)
options.register('globalTag', 'NOTSET',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set global tag"
)
options.register('wantSummary', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('wantFullRECO', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('reportEvery', 100,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "report every N events"
)
options.register('skip', 0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "skip first N events"
)

#options.setDefault('maxEvents', -1)
options.setDefault('maxEvents', 1000)
tag = '2022_SV' if options.isPreECALleakage else '2022EE_metFilter'
options.setDefault('tag', tag)
options.parseArguments()

# global tags:
#    (!) december 2023 - reminiAOD for eras 
# MC pre ECAL leakage  : 130X_mcRun3_2022_realistic_v5           
# MC post ECAL leakage : 130X_mcRun3_2022_realistic_postEE_v6    
# 2022 ABCDE ReReco    : 130X_dataRun3_v2                        # for the time being, not available in 130X for parking
# 2022 FG Prompt       : 130X_dataRun3_PromptAnalysis_v1
if not options.isMC :
    globaltag = '130X_dataRun3_PromptAnalysis_v1'
 
else :
    globaltag = '130X_mcRun3_2022_realistic_v5' if options.isPreECALleakage else '130X_mcRun3_2022_realistic_postEE_v6'


if options._beenSet['globalTag']:
    globaltag = options.globalTag

extension = {False : 'data', True : 'mc'}
outputFileNANO = cms.untracked.string('_'.join(['tau3muNANO', extension[options.isMC], options.tag])+'.root')
outputFileFEVT = cms.untracked.string('_'.join(['xFullEvt', extension[options.isMC], options.tag])+'.root')
if not options.inputFiles :
    if options.isMC :
        # signal channel
        options.inputFiles = ['/store/mc/Run3Summer22MiniAODv4/WtoTauNu_Tauto3Mu_TuneCP5_13p6TeV_pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2540000/956f1823-037d-4a9c-aa2f-50dcf5936f83.root'] if options.isPreECALleakage else \
                            ['/store/mc/Run3Summer22EEMiniAODv4/WtoTauNu_Tauto3Mu_TuneCP5_13p6TeV_pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/2530000/6965991a-4aea-4d25-84b9-82404f0d2b64.root']
        # control channel
        #options.inputFiles = ['/store/mc/Run3Summer22MiniAODv3/DstoPhiPi_Phito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen/MINIAODSIM/124X_mcRun3_2022_realistic_v12-v2/2810000/0da9edba-f8b9-4e0c-8be1-282cdd2b5685.root'] if options.isPreECALleakage else \
        #                     ['/store/mc/Run3Summer22EEMiniAODv3/DstoPhiPi_Phito2Mu_MuFilter_TuneCP5_13p6TeV_pythia8-evtgen/MINIAODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/2810000/00589525-be33-4abd-af78-428bb9ace158.root']


    else :
        options.inputFiles = ['/store/data/Run2022F/ParkingDoubleMuonLowMass0/MINIAOD/22Sep2023-v1/60000/00807cd1-61e0-4a30-917a-f2957e4365b0.root',
                              '/store/data/Run2022F/ParkingDoubleMuonLowMass0/MINIAOD/22Sep2023-v1/60000/0163f2af-82a7-4540-8702-b111818a93d4.root'
        ]
annotation = '%s nevts:%d' % (outputFileNANO, options.maxEvents)

#from Configuration.StandardSequences.Eras import eras
process = cms.Process('Tau3muNANO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.Tau3muNANO.nanoTau3Mu_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents=cms.untracked.uint32(options.skip),
)
process.options = cms.untracked.PSet(
    #Rethrow 
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(options.wantSummary),
)

process.nanoMetadata.strings.tag = annotation
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string(annotation),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = outputFileNANO,
    outputCommands = cms.untracked.vstring(
      'drop *',
      "keep nanoaodFlatTable_*Table_*_*",     # event data
      "keep nanoaodUniqueString_nanoMetadata_*_*",   # basic metadata
    )

)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')

from PhysicsTools.Tau3muNANO.nanoTau3Mu_cff import *
process = nanoAOD_customizeMuonTriggerTau3Mu(process)
#process = nanoAOD_customizeTrackTau3Mu(process)
process = nanoAOD_customizeWnuTau3Mu(process)
process = nanoAOD_customizeTriggerBitsTau3Mu(process)

# Path and EndPath definitions
process.nanoAOD_TauTo3mu_step = cms.Path(process.nanoSequence + process.nanoWnuTau3MuSequence )

# customisation of the process.
if options.isMC:
    from PhysicsTools.Tau3muNANO.nanoTau3Mu_cff import nanoAOD_customizeMC
    nanoAOD_customizeMC(process)

process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
process.schedule = cms.Schedule(
                                process.nanoAOD_TauTo3mu_step,
                                process.endjob_step,
                                process.NANOAODoutput_step
                               )

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
                                   'nanoAOD_TauTo3mu_step'
                                   )
)


process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)    

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

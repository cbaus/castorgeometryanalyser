# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: Configuration/Generator/python/SingleElectronE1000_cfi.py --mc -s GEN,SIM,DIGI,L1,DIGI2RAW,RAW2DIGI,RECO --eventcontent RECOSIM --relval None --datatier GEN-SIM-RECO --conditions auto:mc -n 5 --no_exec
import sys
import FWCore.ParameterSet.Config as cms

from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper

if len(sys.argv) < 3:
        sys.exit('Usage: %s file-name' % sys.argv[0])
print sys.argv
process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.GeometryAll_cff') # added
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5000)
)
#test
print "here it comes"
print XMLIdealGeometryESSource
exit(-1)
#test
# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('Configuration/Generator/python/SingleElectronE1000_cfi.py nevts:5'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string(sys.argv[2]),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

process.generator = cms.EDProducer("FlatRandomEGunProducer",
    PGunParameters = cms.PSet(
        PartID = cms.vint32(211),
        MaxPhi = cms.double(3.14159265359),
        MinPhi = cms.double(-3.14159265359),
        MinEta = cms.double(-6.8),
        MaxEta = cms.double(-5.0),
        MinE = cms.double(29.5),
        MaxE = cms.double(30.5)
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('single electron E 1000'),
    AddAntiParticle = cms.bool(False),
    firstRun = cms.untracked.uint32(1)
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)


# # M E S S A G I N G
# #Module Tracing
# process.Timing = cms.Service("Timing")
# # initialize MessageLogger and output report
# process.load("FWCore.MessageLogger.MessageLogger_cfi")
# process.MessageLogger = cms.Service("MessageLogger",
#                                     destinations=
#                                     cms.untracked.vstring('detailedInfo', 'critical', 'cerr','debuglog'),
#                                     critical=cms.untracked.PSet(threshold = cms.untracked.string('ERROR')),
#                                     detailedInfo=cms.untracked.PSet(threshold = cms.untracked.string('INFO')),
#                                     cerr=cms.untracked.PSet(threshold = cms.untracked.string('WARNING')),
#                                     debuglog=cms.untracked.PSet(threshold = cms.untracked.string('DEBUG'))
#                                     )


# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.raw2digi_step,process.reconstruction_step,process.endjob_step,process.RECOSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq

def customise(process):
    # G E A N T 4   C O N F I G
    process.g4SimHits.Generator.MinEtaCut =-7.0
    process.g4SimHits.Generator.MaxEtaCut = -5.0
    process.g4SimHits.CastorSD.useShowerLibrary = cms.bool(False)
    return(process)

process = customise(process)


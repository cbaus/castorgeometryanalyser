import FWCore.ParameterSet.Config as cms
import os
from subprocess import Popen, PIPE

type='igor'
dir = '/eos/cms/store/temp/user/cbaus/castorgeotest/'+type+'/'
#p = Popen(["/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select", "ls", dir], stdout=PIPE)
#output = ['root://eoscms/' + dir + f[:-1] for f in p.stdout.readlines()]

output = ['root://eoscms//eos/cms/store/group/phys_heavyions/katkov/ttbar_13_710p8_self/step3_RAW2DIGI_L1Reco_RECO_EI01.root']
print output

#files2 = []
#for root, dirs, files in os.walk('/afs/cern.ch/user/c/cbaus/Geometry/CMSSW_7_1_0_pre3/test/' + file):
#for root, dirs, files in os.walk('/afs/cern.ch/user/c/cbaus/Geometry/Old/CMSSW_7_1_0_pre3/test/' + file):
#   for f in files:
#        files2.append('file:' + os.path.join(root,f))

#print files2
process = cms.Process("Demo")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-50000)
)
process.source = cms.Source("PoolSource",
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(output)
                            )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('Histos_' + type + '.root')
                                   )


process.demo = cms.EDAnalyzer('MyPhiAnalyzer')

process.p = cms.Path(process.demo)

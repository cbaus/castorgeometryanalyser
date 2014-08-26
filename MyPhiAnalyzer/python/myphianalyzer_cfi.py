import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use 
    fileNames = cms.untracked.vstring(
#        'file:~/1/GeometryTest/test/NewGeoNoTransNoTilt.root'
        'file:~/1/GeometryTest/src/Geometry/OldGeo.root'
    )
)
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('histodemo2.root')
                                   )


process.demo = cms.EDAnalyzer('MyPhiAnalyzer')

process.p = cms.Path(process.demo)

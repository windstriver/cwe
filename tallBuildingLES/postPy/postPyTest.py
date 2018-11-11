import postPy

boundaryCloudDir = '/lustre/work/wan39502/' + 'tallBuildingLES/testCase/' + \
                   'postProcessing/boundaryCloudPresBuildingExp/'

fvPatchSfFile = './buildingPatchFaceAreaVectors'

probePositionFile = '/lustre/work/wan39502/' + \
                    'tallBuildingLES/ExpTPU/' + \
                    'Location_of_CFD_pressure_probes.csv'
presCoefExpFile = '/lustre/work/wan39502/' + \
                  'tallBuildingLES/ExpTPU/' + \
                  'presCoefExp.csv'

# internalCloudDir = '/lustre/work/wan39502/' + 'tallBuildingLES/testCase/' + \
#                    'postProcessing/internalCloudUUpstream/'

# bc = postPy.BoundaryCloudP(boundaryCloudDir = boundaryCloudDir, \
#                            fvPatchSfFile = fvPatchSfFile)
bcExp = postPy.BoundaryCloudPExp(probePositionFile=probePositionFile, \
                                 presCoefExpFile=presCoefExpFile, \
                                 fvPatchSfFile=fvPatchSfFile)
# baseCoefDF = bc.baseForceMomentCoef()
baseCoefExpDF = bcExp.baseForceMomentCoef()

# ic = postPy.InternalCloudU(Uh, H, B, D, deltaT, internalCloudDir)

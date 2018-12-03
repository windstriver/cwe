import postPy

Uh = 11.11
H = 0.5
B = 0.2
D = 0.1
deltaT = 2e-5

#boundaryCloudDir = '/lustre/work/wan39502/' + 'tallBuildingLES/testCase/' + \
#                   'postProcessing/boundaryCloudPresBuildingExp/'

#fvPatchSfFile = './buildingPatchFaceAreaVectorsExp'

#probePositionFile = '/lustre/work/wan39502/' + \
#                    'tallBuildingLES/ExpTPU/' + \
#                    'Location_of_CFD_pressure_probes.csv'

#presCoefExpFile = '/lustre/work/wan39502/' + \
#                  'tallBuildingLES/ExpTPU/' + \
#                  'presCoefExp.csv'

internalCloudDir = '/lustre/work/wan39502/' + \
                   'emptyDomainDecreaseTimeStep/testCase/' + \
                   'postProcessing/internalCloudUUpstream/'

#bc = postPy.BoundaryCloudP(boundaryCloudDir = boundaryCloudDir, \
#                           fvPatchSfFile = fvPatchSfFile)

#bcExp = postPy.BoundaryCloudPExp(probePositionFile=probePositionFile, \
#                                 presCoefExpFile=presCoefExpFile, \
#                                 fvPatchSfFile=fvPatchSfFile)

#baseCoefDF = bc.baseForceMomentCoef()

#baseCoefExpDF = bcExp.baseForceMomentCoef()

ic = postPy.InternalCloudU(Uh, H, B, D, deltaT, internalCloudDir)

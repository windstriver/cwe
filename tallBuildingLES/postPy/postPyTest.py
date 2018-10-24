import postPy

Uh = 11.11
H = 0.5
B = 0.2
D = 0.1
deltaT = 2e-4

boundaryCloudDir = '/lustre/work/wan39502/' + 'tallBuildingLES/testCase/' + \
                   'postProcessing/boundaryCloudPresBuildingExp/'
fvPatchSfFile = '/lustre/work/wan39502/'+ 'tallBuildingLES/testCase/' + \
                   'constant/polyMesh/writeMesh/buildingPatchFaceAreaVectors'

internalCloudDir = '/lustre/work/wan39502/' + 'tallBuildingLES/testCase/' + \
                   'postProcessing/internalCloudUUpstream/'

#bc = postPy.BoundaryCloud(Uh, H, B, D, deltaT, boundaryCloudDir, fvPatchSfFile)
ic = postPy.InternalCloudU(Uh, H, B, D, deltaT, internalCloudDir)


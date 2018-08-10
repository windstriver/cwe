#!/usr/bin/env python

"""
boundaryCloudBuilding.py

Read OpenFOAM postProcessing/boundaryCloud directory,
and extract pressure time histories at building patch centres.
The pressure time histories are then saved to csv file, with a format of
Pt No,  time.0, time.1, ..., time.N
"""

import os
import numpy as np
import pandas as pd

# Directory of boundaryCloud for building
postDir = '/lustre/work/wan39502/squareCylinderLES/testCase/postProcessing/boundaryCloudBuilding/'
# List of time directories
timeDirList = os.listdir(postDir)
timeNo = np.size(timeDirList)

# Extract the coordinates of building patch centres
ptList = np.genfromtxt(postDir+timeDirList[0]+'/cloud_p.xy')
ptNo = np.shape(ptList)[0]

# Create pandas DataFrame of building patch centres
# Pt.No,    x,    y,    z
ptDataFrame = pd.DataFrame(ptList[:, :-1], \
                           columns = ['x', 'y', 'z'])

# Save the building patch centres(coordinates) to csv file
ptDataFrame.to_csv('buildingPatchCentres.csv')


# Extract pressure time histories at patch centres
# Pt No,  time.0, time.1, ..., time.N
presDataFrame = pd.DataFrame(np.zeros((ptNo, timeNo)), \
                             columns = timeDirList)

presInletDF = pd.read_csv('presInlet.csv', index_col=0)

for timeDir in timeDirList:
    pTemp = np.genfromtxt(postDir+timeDir+'/cloud_p.xy')
    presDataFrame.loc[:, timeDir] = pTemp[:, 3] - presInletDF.loc[0, timeDir]

# Sort the time histories by increasing time
presDataFrame = presDataFrame.sort_index(axis=1)

# Save pressure time histories into csv file
presDataFrame.to_csv('presBuildingPatchCentres.csv')

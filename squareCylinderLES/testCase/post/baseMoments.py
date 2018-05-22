#!/usr/bin/env python

"""
"""

import numpy as np
import pandas as pd

# Read the building patch centres from CSV
ptDataFrame = pd.read_csv('buildingPatchCentres.csv', index_col=0)

# Read the building patch face area vectors
sfDataFrame = pd.read_csv('../buildingPatchFaceAreaVectors', \
                          header=None, names=['Sfx', 'Sfy', 'Sfz'])

# Read the pressure time histories CSV to pandas DataFrame
presDataFrame = pd.read_csv('presBuildingPatchCentres.csv', index_col=0)

# Create base moment DataFrame
timeNo = presDataFrame.shape[1]
momentDataFrame = pd.DataFrame(np.zeros((3, timeNo)), \
                               index=['Mx','My','Mt'], \
                               columns=presDataFrame.columns)

# Building center
H = 0.3658;
B = 0.06096;
xc = 5*H+B/2;
yc = 5*H+B/2;

# Calculate base moments
for tI in presDataFrame.columns:
    Fx = np.multiply(presDataFrame.loc[:,tI], sfDataFrame.loc[:,'Sfx'])
    Mx = np.sum(np.multiply(Fx, ptDataFrame.loc[:,'z']))
    Fy = np.multiply(presDataFrame.loc[:,tI], sfDataFrame.loc[:,'Sfy'])
    My = np.sum(np.multiply(Fy, ptDataFrame.loc[:,'z']))
    Mt = np.sum(np.multiply(Fx, yc-ptDataFrame.loc[:,'y']) + \
                np.multiply(Fy, ptDataFrame.loc[:,'x']-xc))
    momentDataFrame.loc[:,tI] = [Mx, My, Mt]


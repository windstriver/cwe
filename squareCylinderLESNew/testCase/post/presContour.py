"""
presContour.py

Convert the coordinates of pressure probes on building surface to a plane.

Then plot the mean and r.m.s of pressure coefficient contours

"""

import numpy as np
import pandas as pd

# Read the building patch centres from csv
ptDataFrame = pd.read_csv('buildingPatchCentres.csv', index_col=0)

# Read the building patch face area vectors
sfDataFrame = pd.read_csv('../buildingPatchFaceAreaVectors', \
                          header=None, names=['Sfx', 'Sfy', 'Sfz'])

# Read the pressure time histories CSV to pandas DataFrame
presDataFrame = pd.read_csv('presBuildingPatchCentres.csv', index_col=0)

# Create pressure mean and std distribution dataframe
ptNo = ptDataFrame.shape[0]
coefPresContourDataFrame = pd.DataFrame(np.zeros((ptNo, 5)), \
                                    index=ptDataFrame.index, \
                                    columns=['faceID', 'x', 'y', 'mean', 'std'])

# Building size
H = 0.3658
B = 0.06096
rho = 1.225
Vh = 5

for pt in presDataFrame.index:
    if sfDataFrame.loc[pt, 'Sfx'] > 0:    # windward surface
        coefPresContourDataFrame.loc[pt, 'faceID'] = 0
        coefPresContourDataFrame.loc[pt, 'x'] = np.abs(ptDataFrame.loc[pt, 'y'] - (5*H+B))
        coefPresContourDataFrame.loc[pt, 'y'] = ptDataFrame.loc[pt, 'z']
    elif sfDataFrame.loc[pt, 'Sfy'] < 0:    # side surface 1
        coefPresContourDataFrame.loc[pt, 'faceID'] = 1
        coefPresContourDataFrame.loc[pt, 'x'] = np.abs(ptDataFrame.loc[pt, 'x'] - 5*H -1) + B
        coefPresContourDataFrame.loc[pt, 'y'] = ptDataFrame.loc[pt, 'z']
    elif sfDataFrame.loc[pt, 'Sfx'] < 0:    # leeward surface
        coefPresContourDataFrame.loc[pt, 'faceID'] = 2
        coefPresContourDataFrame.loc[pt, 'x'] = np.abs(ptDataFrame.loc[pt, 'y'] - 5*H) + 2*B
        coefPresContourDataFrame.loc[pt, 'y'] = ptDataFrame.loc[pt, 'z']
    elif sfDataFrame.loc[pt, 'Sfy'] > 0:    # side surface 2
        coefPresContourDataFrame.loc[pt, 'faceID'] = 3
        coefPresContourDataFrame.loc[pt, 'x'] = np.abs(ptDataFrame.loc[pt, 'y'] - 5*H - B - 1) + 3*B
        coefPresContourDataFrame.loc[pt, 'y'] = ptDataFrame.loc[pt, 'z']
    else:    # top surface
        coefPresContourDataFrame.loc[pt, 'faceID'] = 4
        coefPresContourDataFrame.loc[pt, 'x'] = np.abs(ptDataFrame.loc[pt,'x'] - 5*H)
        coefPresContourDataFrame.loc[pt, 'y'] = np.abs(ptDataFrame.loc[pt,'y'] - 5*H)

    coefPresContourDataFrame.loc[pt, 'mean'] = np.mean(presDataFrame.iloc[pt, -3000:-10] / (1/2*rho*Vh**2) )
    coefPresContourDataFrame.loc[pt, 'std'] = np.std(presDataFrame.iloc[pt, -3000:-10] / (1/2*rho*Vh**2) )

    
# Save windward, leeward, and side faces pressure coefficients to csv
coefPresContourDataFrame[coefPresContourDataFrame.faceID != 4].to_csv('coefPresContour.csv')

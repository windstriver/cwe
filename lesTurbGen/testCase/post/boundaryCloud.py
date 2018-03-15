import os
import numpy as np
import pandas as pd

udf = pd.DataFrame(columns=['u', 'v', 'w'])

postDir = '../postProcessing/boundaryCloud/'
timeDirList = os.listdir(postDir)
for timeDir in timeDirList:
    uTemp = np.genfromtxt(postDir+timeDir+'/cloud_U.xy')
    udf.loc[float(timeDir)] = uTemp[3:]

udf = udf.sort_index(axis=0)

udf.to_csv('inlet_velocity.csv')

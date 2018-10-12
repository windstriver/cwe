import os
import numpy as np
import pandas as pd

class BoundaryCloud:
    def __init__(self, Uh, boundaryCloudDir, fvPatchSfFile):
        '''
        Constructor by the OpenFOAM postProcessing/boundaryCloud dir and
        boundary patch face area vectors.
        '''
        # Reference velocity at the building height
        self.Uh = Uh
        
        # List of time directories
        timeDirList = os.listdir(boundaryCloudDir)

        # Extract the coordinates of boundary patch face centers(probes)
        try:
            self.probePosition = np.genfromtxt(boundaryCloudDir + \
                                               timeDirList[0] + \
                                               '/cloud_p.xy', \
                                               usecols = (0, 1, 2))
        except:
            self.probePosition = np.genfromtxt(boundaryCloudDir + \
                                               timeDirList[0] + \
                                               '/points_p.xy', \
                                               usecols = (0, 1, 2))

        # Extract pressure coefficient time histories
        # number of probes
        self.noProbe = np.shape(self.probePosition)[0]
        # number of time steps
        self.noTimeStep = np.size(timeDirList)
        # initialization of the pandas dataFrame to store the
        # pressure coefficient time series
        self.presCoefDF = pd.DataFrame(np.zeros((self.noProbe, self.noTimeStep)), \
                                       columns = timeDirList)
        # extractin
        for timeDir in timeDirList:
            try:
                presCoef = np.genfromtxt(boundaryCloudDir + \
                                         timeDir + \
                                         '/cloud_p.xy', \
                                         usecols = (3))
            except:
                presCoef = np.genfromtxt(boundaryCloudDir + \
                                         timeDir + \
                                         '/points_p.xy', \
                                         usecols = (3))

            presCoef = presCoef / (0.5 * self.Uh**2)
            self.presCoefDF.loc[:, timeDir] = presCoef

        # extraction end
        # sort the time histories by increasing time
        self.presCoefDF = self.presCoefDF.sort_index(axis=1)

        # Boundary patch face area vectors
        self.fvPatchSf = np.genfromtxt(fvPatchSfFile, delimiter=',')


def integralTimeScale(x, deltaT):
    '''
    Estimate the integral time scale by autocorrelation function
    Input:
        x         signal
        deltaT    time step of the signal
    Output:
        integral time scale
    '''
    N = len(x)
    autoCorrFun = np.correlate(x, x, mode='full') / np.correlate(x, x)
    mask = np.ones(len(autoCorrFun), dtype=bool)
    mask[:N-1] = False
    autoCorrFun = autoCorrFun[mask]
    mask = np.ones(len(autoCorrFun), dtype=bool)
    mask[np.argmax(autoCorrFun),:] = False
    autoCorrFun = autoCorrFun[mask]
    
    return np.sum(autoCorrFun) * deltaT



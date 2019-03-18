import os
import numpy as np

class Cloud:
    def __init__(self, cloudDir, varName, deltaT, cutinTime):
        '''
        Constructor by the OpenFOAM postProcessing/*Cloud dir.
        '''
        
        # List of time directories
        timeNameList = np.array(os.listdir(cloudDir))
        # number of time steps
        numTimeStep = np.size(timeNameList)
        # convert time name to value
        timeValueList = timeNameList.astype(np.float)
        # combine into the time array
        dtype = [('timeName','<U10'), ('timeValue',float)]
        self.time = np.array([(timeNameList[i], timeValueList[i]) \
                         for i in range(numTimeStep)], dtype = dtype)
        self.time = np.sort(self.time, order='timeValue')

        # cut-in time index
        cutinTimeIndex = np.floor(cutinTime/deltaT).astype(int) - 1
        mask = np.ones(len(self.time), dtype=bool)
        mask[:cutinTimeIndex] = False
        self.time = self.time[mask]
        # Extract the coordinates of boundary patch face centers(probes)
        
        cloudTemp = np.genfromtxt(cloudDir + \
                                  self.time[-1][0] + \
                                  '/cloud_{:s}.xy'.format(varName))
        self.probeLoc = cloudTemp[:,0:3]

        # number of probes
        self.numProbe = np.shape(self.probeLoc)[0]
        # number of variable components
        numVarComp = np.shape(cloudTemp)[1] - 3

        # initialization of the numpy array to store the sampled data
        self.cloudData = np.zeros((numTimeStep-cutinTimeIndex,self.numProbe,numVarComp))

        # extract the cloud data
        for tI in range(len(self.time)):
            try:
                self.cloudData[tI,:,:] = np.genfromtxt(cloudDir + \
                                            self.time[tI][0] + \
                                            '/cloud_{:s}.xy'.format(varName), \
                                            usecols = tuple(range(3,numVarComp+3)))
                
            except ValueError:
                self.cloudData[tI,:,:] = np.nan

        # extraction end



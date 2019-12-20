#!/usr/bin/env python

"""
massFluxCorrection.py

Inlet mass flux correction (Kim et al., 2013)

"""

import numpy as np
import h5py

inputHDF5File = '../turbGen/inflowTurb.h5'
outputHDF5File = '../inflowTurbMFC.h5'
ofCase = '../testCase/'

# building height
H = 0.5
# computational domain height
Z = 10*H
# computational domain width
Y = 4*H

# delta = 1.0;

##############################################################################
# Calculate the prescribed bulk velocity from mean wind velocity profile
##############################################################################
# reference height for the mean velocity profile
Href = 0.5
# mean velocity at the reference height
Uref = 11.1438
# power law exponent of the mean velocity profile
alphau = 1/4;
# bulk velocity
Ub = 1/(alphau+1) * 4**(alphau) * Uref
# y0 = 0.02/300;
# Ub = (np.log(4*H/y0)-1)/np.log(Href/y0)*Uref
print(f'Prescribed bulk velocity is {Ub:.3f} m/s')

##############################################################################
# Calculate the instantaneous bulk velocity from synthetic velocities
##############################################################################
# read the inlet patch face area vectors
sfInlet = np.genfromtxt(ofCase+\
                        'constant/meshInfo/faceAreaVectorsInlet.csv',\
                        delimiter=',')
sMagInlet = np.abs(sfInlet[:, 0])
S = np.sum(sMagInlet)
print('Area of inlet patch is {:.5f}({:.5f}) m^2'.\
      format(S, Y*Z))
# read synthetic velocities from HDF5 file
h5FileIn = h5py.File(inputHDF5File, 'r')
uIn = h5FileIn['U']
vIn = h5FileIn['V']
wIn = h5FileIn['W']
uMeanIn = h5FileIn['UMEAN']
# shape of the velocity matrix
nPt, nTime = uIn.shape
print(f'Shape of input velocities, nPt = {nPt}, nTime = {nTime}')
# total velocity = mean + turb
uTot = uIn[:,:] + uMeanIn[:][:,None]
# instantaneous bulk velocity
UbT = np.dot(uTot.T, sMagInlet) / S
print('instantaneous bulk velocity:')
print(UbT)
np.savetxt('UbT.csv', UbT)

##############################################################################
# Correct the instantaneous velocity
##############################################################################
uOut = Ub * ((uIn[:,:]+uMeanIn[:][:,None])/UbT) - uMeanIn[:][:,None]
vOut = Ub * (vIn[:,:]/UbT)
wOut = Ub * (wIn[:,:]/UbT)
# corrected instantaneous bulk velocity
UbTCorrected = np.dot((uOut[:,:]+uMeanIn[:][:,None]).T, sMagInlet) / S
print('Corrected instantaneous bulk velocity:')
print(UbTCorrected)

##############################################################################
# Write the corrected velocity to HDF5 file
##############################################################################
h5FileOut = h5py.File(outputHDF5File, 'w')
dt = np.dtype('>d')

dset = h5FileOut.create_dataset('U', shape=(nPt, nTime), dtype=dt)
dset.write_direct(uOut)

dset = h5FileOut.create_dataset('V', shape=(nPt, nTime), dtype=dt)
dset.write_direct(vOut)

dset = h5FileOut.create_dataset('W', shape=(nPt, nTime), dtype=dt)
dset.write_direct(wOut)

dset = h5FileOut.create_dataset('UMEAN', shape=(uMeanIn.shape[0],), dtype=dt)
dset.write_direct(uMeanIn[:])

h5FileIn.close()
h5FileOut.close()


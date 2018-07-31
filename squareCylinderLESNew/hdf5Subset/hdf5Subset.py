import numpy as np
import h5py

fileIn = h5py.File('../lesinletfull.h5', 'r')
fileOut = h5py.File('lesinlet.h5', 'w')

u = fileIn['U']
v = fileIn['V']
w = fileIn['W']
Umean = fileIn['UMEAN']

N = 6001

dt = np.dtype('>d')

dset = fileOut.create_dataset('U', shape=(u.shape[0], N), dtype=dt)
dset.write_direct(u[:, :N])

dset = fileOut.create_dataset('V', shape=(u.shape[0], N), dtype=dt)
dset.write_direct(v[:, :N])

dset = fileOut.create_dataset('W', shape=(u.shape[0], N), dtype=dt)
dset.write_direct(w[:, :N])

dset = fileOut.create_dataset('UMEAN', shape=(Umean.shape[0],), dtype=dt)
dset.write_direct(Umean[:])

fileOut.close()

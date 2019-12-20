import numpy as np
import h5py

fileIn = h5py.File('../inflowTurbMFC.h5', 'r')

u = fileIn['U']
v = fileIn['V']
w = fileIn['W']
Umean = fileIn['UMEAN']

N = 5000

dt = np.dtype('>d')

# 0 to 1
for i in range(6):
    fileOut = h5py.File('../inflowTurbMFC{}to{}.h5'\
        .format(str(i), str(i+1)), 'w')

    dset = fileOut.create_dataset('U', shape=(u.shape[0], N+1), dtype=dt)
    dset.write_direct(u[:, i*N:N+1+i*N])

    dset = fileOut.create_dataset('V', shape=(u.shape[0], N+1), dtype=dt)
    dset.write_direct(v[:, i*N:N+1+i*N])

    dset = fileOut.create_dataset('W', shape=(u.shape[0], N+1), dtype=dt)
    dset.write_direct(w[:, i*N:N+1+i*N])

    dset = fileOut.create_dataset('UMEAN', shape=(Umean.shape[0],), dtype=dt)
    dset.write_direct(Umean[:])
    
    fileOut.close()



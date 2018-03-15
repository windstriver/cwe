import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import h5py
import wind as wp
from matplotlib import rc

rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)
rc('lines', linewidth=0.5)

case = "testCase/"
fname = case+"postProcessing/probes/0/U"
foname = fname+'.out'
with open(fname, 'rt') as fin:
    with open(foname, 'wt') as fout:
        for line in fin:
            fout.write(line.replace('(', ' ').replace(')', ' '))

timeHistDat = np.genfromtxt(foname, comments='#')
time = timeHistDat[:, 0]
# U time history at inlet top point
UTopInlet = timeHistDat[:, 1]
UTopBuild = timeHistDat[:, 4]

fh5 = h5py.File('lesinlet.h5', 'r')
UTopHDF5 = fh5['U'][6, :]

UBCDat = np.genfromtxt('boundaryCloud.csv', delimiter=',', skip_header=1)
uInlet = UBCDat[:, 1]
uInletTurb = uInlet[-2048:] - np.mean(uInlet[-2048:])

# Plot the U time history at inlet and building top
fig, ax = plt.subplots()
ax.plot(time[-2048:], uInletTurb[-2048:], 'r', label='Inlet')
# ax.plot(time, UTopBuild-np.mean(UTopInlet), 'b', label='Building Top')
ax.plot(fh5['TIME'][-2048:], UTopHDF5[-2048:]-np.mean(UTopHDF5[-2048:]), 'k', label='Input')
legend = ax.legend()
ax.set_title(r'$U$ Time History')
# ax.set_xlim([0, 6])
ax.set_xlabel(r'$\mathrm{Time} (\mathrm{s})$')
ax.set_ylabel(r'$U (\mathrm{m/s})$')
plt.savefig(case+'postProcessing/u_ts.pdf')

# Spectrum of the velocity time history
dt = 2e-3
fs = 1/dt
H = 0.75

fInlet, SuInlet = signal.periodogram(uInletTurb, fs)
fBuild, SuBuild = signal.periodogram(UTopBuild-np.mean(UTopBuild), fs)
# fHDF5, SuHDF5 = signal.periodogram(UTopHDF5[-4096:] -
                                   # np.mean(UTopHDF5[-4096:]), fs)
SuKarmon = wp.Su(H, fInlet)

fig, ax = plt.subplots()
ax.loglog(fInlet, SuKarmon, 'k', label='von Karmon')
ax.loglog(fInlet, SuInlet, 'r', label='Inlet')
# ax.loglog(fBuild, SuBuild, 'b', label='Building')
# ax.loglog(fHDF5, SuHDF5, 'g', label='Input')
ax.set_xlim([1, 100])
ax.set_ylim([1e-6, 10])
# ax.set_title(r'$U$ Spectrum')
ax.set_xlabel(r'$f(\mathrm{Hz})$')
ax.set_ylabel(r'$S_u(\mathrm{m^2 s^2 / Hz})$')
ax.legend()

plt.savefig(case+'postProcessing/u_spectrum.pdf')

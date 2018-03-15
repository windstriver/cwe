import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import signal
import wind

# parameters
H = 0.5
deltaT = 2e-3

# Data conversion
# import probes
# import boundaryCloud

# matplotlib plot setting
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)
rc('lines', linewidth=0.5)

# extract data
probesDat = np.genfromtxt('internal_cell_velocity.csv', comments='#')
boundaryCloudDat = np.genfromtxt('inlet_velocity.csv',
                                 delimiter=',', skip_header=1)
time = probesDat[:, 0]
uInlet = boundaryCloudDat[:, 1]
uInletCellCen = probesDat[:, 1]
uBuild = probesDat[:, 4]

# velocity spectrum
time = time[-2048:]
uInlet = uInlet[-2048:] - np.mean(uInlet[-2048:])
uInletCellCen = uInletCellCen[-2048:] - np.mean(uInletCellCen[-2048:])
uBuild = uBuild[-2048:] - np.mean(uBuild[-2048:])

fs = 1/deltaT
fuInlet, SuInlet = signal.periodogram(uInlet, fs)
fuInletCellCen, SuInletCellCen = signal.periodogram(uInletCellCen, fs)
fuBuild, SuBuild = signal.periodogram(uBuild, fs)

SuKarmon = wind.Su(H, fuInlet)

fig, ax = plt.subplots()
ax.loglog(fuInlet, SuKarmon, 'k', label='von Karmon')
ax.loglog(fuInlet, SuInlet, 'r', label='Inlet')
ax.loglog(fuInletCellCen, SuInletCellCen, 'b', label='Inlet Cell Center')
# ax.loglog(fuBuild, SuBuild, 'g', label='Building')
ax.set_xlim([1, 100])
ax.set_ylim([1e-6, 10])
ax.set_xlabel(r'$f(\mathrm{Hz})$')
ax.set_ylabel(r'$S_u(\mathrm{m^2 s^2 / Hz})$')
ax.legend()
plt.savefig('u_spectrum.pdf')

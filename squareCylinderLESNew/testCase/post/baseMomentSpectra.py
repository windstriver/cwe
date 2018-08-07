"""
baseMomentSpectra.py

Calculate the spectra of base moment time history.

"""

import numpy as np
from scipy import signal
import matplotlib as mpl
mpl.use('SVG')
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)

# Read the base moment time history from csv
baseMomentData = np.genfromtxt('baseMomentTimeHistoryNew.csv', \
                               delimiter=',', skip_header=1)

time = baseMomentData[-5000:,0]
My = baseMomentData[-5000:,1]
Mx = baseMomentData[-5000:,2]
Mt = baseMomentData[-5000:,3]

# Estimate power spectral density using a periodogram
deltaT = 1e-3
fs = 1/deltaT    # sampling frequency
B = 0.06096      # building width
Vh = 5           # reference velocity at building height

# Read experimental PSD data
PSD_My = np.genfromtxt('PSD_along_wind_base_moment.csv', \
                       delimiter=',')
PSD_Mx = np.genfromtxt('PSD_cross_wind_base_moment.csv', \
                       delimiter=',')
PSD_Mt = np.genfromtxt('PSD_torsional_base_moment.csv', \
                       delimiter=',')

# Along wind base moment power spectral density
freq, SpectraMy = signal.welch(My, fs, window='hann', \
                               nperseg=1024, noverlap=512)
fig = plt.figure(num=None, figsize=(8,3.5))
fig.subplots_adjust(wspace=0.5)

ax = fig.add_subplot(131)
ax.loglog(freq*B/Vh, np.multiply(SpectraMy, freq)/np.var(My), \
          'b-', label='LES')
ax.loglog(PSD_My[:,0], PSD_My[:,1], 'r--', label='Exp.')
ax.set_xlim([np.min(PSD_My[:,0]), np.max(PSD_My[:,0])])
ax.set_ylim([1e-3, 1])
ax.set_aspect(0.8)
ax.set_xlabel(r'$fB/V_h$')
ax.set_ylabel(r'$S_{My} f/\sigma_{My}^2$')
ax.set_title(r'$M_y$')

# Cross wind base moment power spectral density
freq, SpectraMx = signal.welch(Mx, fs, window='hann', \
                               nperseg=1024, noverlap=512)
ax = fig.add_subplot(132)
ax.loglog(freq*B/Vh, np.multiply(SpectraMx, freq)/np.var(Mx), 'b-')
ax.loglog(PSD_Mx[:,0], PSD_Mx[:,1], 'r--')
ax.set_xlim([np.min(PSD_My[:,0]), np.max(PSD_My[:,0])])
ax.set_ylim([1e-3, 1])
ax.set_aspect(0.8)
ax.set_xlabel(r'$fB/V_h$')
ax.set_ylabel(r'$S_{Mx} f/\sigma_{Mx}^2$')
ax.set_title(r'$M_x$')

# Torsional moment power spectral density
freq, SpectraMt = signal.welch(Mt, fs, window='hann', \
                               nperseg=1024, noverlap=512)
ax = fig.add_subplot(133)
ax.loglog(freq*B/Vh, np.multiply(SpectraMt, freq)/np.var(Mt), \
          'b-')
ax.loglog(PSD_Mt[:,0], PSD_Mt[:,1], 'r--')
ax.set_xlim([np.min(PSD_My[:,0]), np.max(PSD_My[:,0])])
ax.set_ylim([1e-3, 1])
ax.set_aspect(0.8)
ax.set_xlabel(r'$fB/V_h$')
ax.set_ylabel(r'$S_{Mt} f/\sigma_{Mt}^2$')
ax.set_title(r'$M_t$')

fig.legend(bbox_to_anchor=(0.5,0), loc='lower center', ncol=2)

fig.savefig('PSD_base_moments.svg', bbox_inches='tight')


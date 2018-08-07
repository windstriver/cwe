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
baseMomentData = np.genfromtxt('baseMomentTimeHistory.csv', \
                               delimiter=',', skip_header=1)

time = baseMomentData[-5000:,0]
My = baseMomentData[-5000:,1]
Mx = baseMomentData[-5000:,2]
Mt = baseMomentData[-5000:,3]

# Estimate power spectral density using a periodogram
deltaT = 1e-3
fs = 1/deltaT    # sampling frequency
B = 0.06096      # building width
D = B;
H = 0.3658
Vh = 5           # reference velocity at building height
rho = 1.225;

# Plot base moments time history
Myref = 1/2*rho*Vh**2*B*H**2
Mxref = 1/2*rho*Vh**2*D*H**2
Mtref = 1/2*rho*Vh**2*D*B*H

fig = plt.figure()
ax = fig.add_subplot(311)
ax.plot(time, My/Myref)
ax.set_xlabel('Time (s)')
ax.set_ylabel(r'$M_y/M_{yref}$')
ax.set_title('Along-wind base moment time history')

ax = fig.add_subplot(312)
ax.plot(time, Mx/Mxref)
ax.set_xlabel('Time (s)')
ax.set_ylabel(r'$M_x/M_{xref}$')
ax.set_title('Cross-wind base moment time history')

ax = fig.add_subplot(313)
ax.plot(time, Mt/Mtref)
ax.set_xlabel('Time (s)')
ax.set_ylabel(r'$M_t/M_{tref}$')
ax.set_title('Torsional base moment time history')

fig.savefig('base_moments_time_history.svg')


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


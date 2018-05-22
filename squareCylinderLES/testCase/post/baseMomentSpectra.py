#!/usr/bin/env python

"""
baseMomentSpectra.py

Calculate the spectra of base moment time history.

"""

import numpy as np
from scipy import signal
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# Read the base moment time history from csv
baseMomentData = np.genfromtxt('baseMomentTimeHistory.csv', \
                               delimiter=',', skip_header=1)

time = baseMomentData[-3000:,0]
Mx = baseMomentData[-3000:,1]
My = baseMomentData[-3000:,2]
Mt = baseMomentData[-3000:,3]

# Estimate power spectral density using a periodogram
deltaT = 1e-3
fs = 1/deltaT    # sampling frequency
B = 0.06096      # building width
Vh = 5           # reference velocity at building height

# Along wind base moment power spectral density
freq, SpectraMx = signal.periodogram(Mx, fs)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(freq*B/Vh, np.multiply(SpectraMx, freq)/np.var(Mx))
ax.set_xlim([1e-3, 1])
ax.set_xlabel(r'$fB/V_h$')
ax.set_ylabel(r'$S_{Mx} f/\sigma_{Mx}^2$')
ax.set_title('Along Wind Base Moment PSD')
fig.savefig('spectraMx.pdf')

# Cross wind base moment power spectral density
freq, SpectraMy = signal.periodogram(My, fs)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(freq*B/Vh, np.multiply(SpectraMy, freq)/np.var(My))
ax.set_xlim([1e-3, 1])
ax.set_xlabel(r'$fB/V_h$')
ax.set_ylabel(r'$S_{My} f/\sigma_{My}^2$')
ax.set_title('Cross Wind Base Moment PSD')
fig.savefig('spectraMy.pdf')

# Torsional moment power spectral density
freq, SpectraMt = signal.periodogram(Mt, fs)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(freq*B/Vh, np.multiply(SpectraMt, freq)/np.var(Mt))
ax.set_xlim([1e-3, 1])
ax.set_xlabel(r'$fB/V_h$')
ax.set_ylabel(r'$S_{Mt} f/\sigma_{Mt}^2$')
ax.set_title('Torsional Moment PSD')
fig.savefig('spectraMt.pdf')

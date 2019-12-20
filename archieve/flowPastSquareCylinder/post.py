import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

data = np.genfromtxt('forceCoeffs.dat', skip_header=9)
time = data[:, 0]
Cd = data[:, 2]
Cl = data[:, 3]

# Cd_avg = signal.medfilt(Cd, 1001)
Cd_avg = np.convolve(Cd, np.ones((2000,))/2000, mode='valid')

figCd = plt.figure()
ax = figCd.add_subplot(111)
line1, = ax.plot(time, Cd, label='intantaneous')
line2, = ax.plot(time[1999:], Cd_avg, 'r', label='average - 2000 iterations')
ax.set_ylim([1.5, 2.5])
ax.set_xlim([2.2, 20])
ax.set_xlabel('Time')
ax.set_ylabel(r'$C_D$')
ax.legend()
figCd.savefig('Cd-hist.pdf')

figCl = plt.figure()
ax = figCl.add_subplot(111)
ax.plot(time, Cl)
ax.set_xlabel('Time')
ax.set_ylabel(r'$C_L$')
figCl.savefig('Cl-hist.pdf')

N = 22829
dt = np.sum(time[1:] - time[:N-1])/(N-1)
fs = 1/dt

f, Pxx_den = signal.welch(Cl, fs, nperseg=N)

D = 0.04
U = 0.54
f_norm = f*D/U

figSPD = plt.figure()
ax = figSPD.add_subplot(111)
ax.plot(f_norm, Pxx_den)
ax.set_xlabel(r'$\frac{nD}{U}$')
ax.set_ylabel('SPD')
ax.set_xlim([0, 0.01*f_norm[-1]])
figSPD.savefig('SPD_CL.pdf')

St = f_norm[np.argmax(Pxx_den)]
print('Strouhal value is: ', St)

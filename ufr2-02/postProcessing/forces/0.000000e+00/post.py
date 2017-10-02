import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

data = np.genfromtxt('forceCoeffs.dat', skip_header=9)
time = data[:, 0]
Cl = data[:, 3]

# plt.plot(time, Cl)
# plt.xlabel('Time')
# plt.ylabel('Cl')
# plt.show()

N = 22829
dt = np.sum(time[1:] - time[:N-1])/(N-1)
fs = 1/dt

f, Pxx_den = signal.welch(Cl, fs, nperseg=N)

D = 0.04
U = 0.54
f_norm = f*D/U

plt.plot(f_norm, Pxx_den)
plt.xlabel(r'$\frac{nD}{U}$')
plt.ylabel('SPD')
plt.xlim([0, 0.01*f_norm[-1]])
plt.show()

St = f_norm[np.argmax(Pxx_den)]
print('Strouhal value is: ', St)

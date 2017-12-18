import numpy as np
import matplotlib.pyplot as plt
import h5py

file = "timeHist.csv"
data = np.genfromtxt(file, skip_header=1, delimiter=',')
time = data[:, 3]
Ux = data[:, 5]
Uy = data[:, 6]
Uz = data[:, 7]

f = h5py.File('lesinlet.h5', 'r')
# dsetUx = f.require_dataset('U', (100, 2001), '>f8')

# fig = plt.figure()
# ax = fig.add_subplot(111)
# line1, = ax.plot(time, Ux, 'r', label='Ux')
# # line2, = ax.plot(time, Uy, 'b', label='Uy')
# # line3, = ax.plot(time, Uz, 'g', label='Uz')
# ax.set_xlabel('Time (s)')
# ax.set_ylabel('Velocity (m/s)')
# ax.legend()
# fig.savefig('timeHist.pdf')



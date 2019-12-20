import numpy as np
import matplotlib.pyplot as plt

fileInput = ["U_epsilon_k_over_line_x0.csv",
             "U_epsilon_k_over_line_x2500.csv",
             "U_epsilon_k_over_line_x4000.csv"]

flowVariable = ['U', 'k', 'epsilon']

fileNum = len(fileInput)

data = np.genfromtxt(fileInput[0], delimiter=',', skip_header=1)
pointNum = np.shape(data)[0]

z = np.zeros((pointNum, fileNum))
U = np.zeros((pointNum, fileNum))
epsilon = np.zeros((pointNum, fileNum))
k = np.zeros((pointNum, fileNum))

for i in range(fileNum):
    data = np.genfromtxt(fileInput[i], delimiter=',', skip_header=1)
    z[:, i] = data[:, 3]
    U[:, i] = data[:, 4]
    epsilon[:, i] = data[:, 9]
    k[:, i] = data[:, 10]

# Plot of U
fig, ax = plt.subplots(1, 1)
line0, = ax.plot(U[:,0], z[:,0], 'r-', label='RH')
line2500, = ax.plot(U[:,1], z[:,1], 'b--', linewidth=1.0, label="x=2500m")
line4000, = ax.plot(U[:,2], z[:,2], 'g-.', linewidth=1.0, label="x=4000m")
ax.set_xlim([4,18])
ax.set_ylim([0, 300])
ax.set_xlabel(r'$U$ (m/s)')
ax.set_ylabel(r'$z$ (m)')
ax.legend(loc='best')
fig.savefig("U.pdf")
plt.close()

# Plot of k
fig, ax = plt.subplots(1, 1)
line0, = ax.plot(k[:,0], z[:,0], 'r-', label='RH')
line2500, = ax.plot(k[:,1], z[:,1], 'b--', linewidth=1.0, label="x=2500m")
line4000, = ax.plot(k[:,2], z[:,2], 'g-.', linewidth=1.0, label="x=4000m")
ax.set_xlim([1, 2])
ax.set_ylim([0, 300])
ax.set_xlabel(r'$k$ ($\mathrm{m}^2/\mathrm{s}^2$)')
ax.set_ylabel(r'$z$ (m)')
ax.legend(loc='best')
fig.savefig("k.pdf")
plt.close()


# Plot of epsilon
fig, ax = plt.subplots(1, 1)
line0, = ax.semilogx(epsilon[:,0], z[:,0], 'r-', label='RH')
line2500, = ax.semilogx(epsilon[:,1], z[:,1], 'b--', linewidth=1.0, label="x=2500m")
line4000, = ax.semilogx(epsilon[:,2], z[:,2], 'g-.', linewidth=1.0, label="x=4000m")
# ax.set_xlim([1, 2])
ax.set_ylim([0, 300])
ax.set_xlabel(r'$\varepsilon$ ($\mathrm{m}^2/\mathrm{s}^3$)')
ax.set_ylabel(r'$z$ (m)')
ax.legend(loc='best')
fig.savefig("epsilon.pdf")
plt.close()


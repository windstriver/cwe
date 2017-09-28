import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('p', skip_header=8)

time = data[-100:, 0]
p_0 = data[-100:, 1]

fig = plt.figure()
axes = fig.add_subplot(111)
line, = axes.plot(time, p_0)
axes.set_xlabel('Time')
axes.set_ylabel('p')
fig.savefig('p.png')

# Plot the spectral dencity of wind velocity

x = data[:, 1]

Fs = 1e5
delta_t = 1.0/Fs
N = len(x)
Td = N*delta_t
delta_f = 1.0/Td

t = np.arange(0, N)*delta_t    # time array
f = np.arange(0, N/2)*delta_f

x_fft = np.fft.fft(x)
# print(np.size(x_fft))
P2 = abs(x_fft) / N
P1 = P2[0:int(N/2)]
P1[1:int(N/2)] = (2*P1[1:int(N/2)])**2 / (2*delta_f)

f_norm = f / np.mean(x)
P1_norm = P1*f/np.var(x)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.semilogx(f_norm, P1_norm)
ax2.set_xlabel(r'$f$')
ax2.set_ylabel(r'$S_x(f)$')
ax2.set_title('Spectral density of pressure')
fig2.savefig('p_spetral.png')

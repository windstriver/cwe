import numpy as np
import matplotlib.pyplot as plt

xyFileDir = './postProcessing/singleGraph/100/line_sigmaxx.xy'

xyData = np.genfromtxt(xyFileDir)
y = xyData[:, 0]
Sxx = xyData[:, 1]

# Analytical solution for the stress at x = 0
def Sxx_theory(y):
    R = 0.5
    sigma0 = 10e3
    return sigma0*(1+R**2/(2*y**2)+3*R**4/(2*y**4))

Sxx2 = Sxx_theory(y)

plt.plot(y, Sxx, 'bo', y, Sxx2, 'g-')
plt.xlabel(r'Distance, $y$ (m)')
plt.ylabel(r'Stress $(\sigma_{xx})_{x=0}\ (\mathrm{kPa})$')
plt.xlim([0.5, 2.0])
plt.savefig("sigmaxx_x_0.pdf")

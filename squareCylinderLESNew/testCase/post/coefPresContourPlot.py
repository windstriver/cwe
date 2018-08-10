"""
Plot the pressure coefficient contours

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.patches import Rectangle

coefPresContourDataFrame = pd.read_csv('coefPresContour.csv', index_col=0)
H = 0.3658
B = 0.06096
Vh = 5

# Plot of windward pressure coefficient
coefPresWindward = coefPresContourDataFrame[coefPresContourDataFrame.faceID == 2]
x = coefPresWindward.loc[:, 'x'] - 2*B
y = coefPresWindward.loc[:, 'y']

coefPresMean = coefPresWindward.loc[:, 'mean']
coefPresStd = coefPresWindward.loc[:, 'std']

fig, (ax1, ax2) = plt.subplots(ncols = 2)
#fig.subplots_adjust(wspace=0.1)

#Rectangle((0,0), B, H)
xi = np.linspace(0, B, 20)
yi = np.linspace(0, H, 120)

triang = tri.Triangulation(x, y)
interpolatorMean = tri.LinearTriInterpolator(triang, coefPresMean)
Xi, Yi = np.meshgrid(xi, yi)
zi = interpolatorMean(Xi, Yi)
cs = ax1.contour(xi, yi, zi, 14, linewidths=0.5, colors='k')
#plt.clabel(cs, inline=1)
cntr1 = ax1.contourf(xi, yi, zi, 14, cmap="RdBu_r")
fig.colorbar(cntr1, ax=ax1)
ax1.plot(x, y, 'ko', ms=1)
ax1.axis((0, B, 0, H))
ax1.set_aspect(0.6)
ax1.set_title(r'$C_p$ mean')

interpolatorStd = tri.LinearTriInterpolator(triang, coefPresStd)
zi = interpolatorStd(Xi, Yi)
ax2.contour(xi, yi, zi, 14, linewidths=0.5, colors='k')
cntr2 = ax2.contourf(xi, yi, zi, 14, cmap="RdBu_r")
fig.colorbar(cntr2, ax=ax2)
ax2.plot(x, y, 'ko', ms=1)
ax2.axis((0, B, 0, H))
ax2.set_aspect(0.6)
ax2.set_title(r'$C_p$ std')

fig.savefig('pres_coef_leeward.svg')

# Plot the pressure coefficient time history
# err = B/8
# presDataFrame = pd.read_csv('presBuildingPatchCentres.csv', index_col=0)
# ptIndex = coefPresWindward.loc[(np.abs(coefPresWindward.x-B/2)<err) & (np.abs(coefPresWindward.y-2*H/3)<err)].index.values[0]
# time = presDataFrame.columns.values.astype(float)
# Cp = presDataFrame.loc[ptIndex, :].values / (1/2*Vh**2)

# fig, ax = plt.subplots()
# ax.plot(time[-3000:-10], Cp[-3000:-10])
# ax.set_xlabel('Time (s)')
# ax.set_ylabel(r'$C_p$')
# fig.savefig('cp_side.svg')

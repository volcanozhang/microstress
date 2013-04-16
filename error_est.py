import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pylab

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

fig = pylab.figure()
ax = Axes3D(fig)

name = '101x101'
ini = np.loadtxt(name+'.pts', skiprows = 1)[:, 2:4]
res = np.loadtxt(name+'.stp', skiprows = 1)[:, 2:4]
dd = res - ini
dd[:, 1] = -dd[:, 1]
ddx = dd[:, 0].reshape(101, 101)
ddy = dd[:, 1].reshape(101, 101)

ddx_, ddy_ = np.mgrid[0:1:complex(0,101), 0:1:complex(0,101)]
errx, erry = ddx-ddx_, ddy-ddy_
surf = ax.plot_surface(ddx_, ddy_, errx, rstride=1, cstride=1, cmap=cm.jet)

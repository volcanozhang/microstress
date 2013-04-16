import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pylab

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

from common import stringint

fig = pylab.figure()
ax = Axes3D(fig)
ddx_, ddy_ = np.mgrid[0:1:complex(0,101), 0:1:complex(0,101)]

dd = np.zeros((10, 10201, 2))
mean = np.zeros((10201, 2))
for i in range(10):
    name = stringint(i+1, 3)
    ini = np.loadtxt('101x101.pts', skiprows = 1)[:, 2:4]
    res = np.loadtxt(name, skiprows = 1)[:, 2:4]
    dd[i] = res - ini
    dd[i, :, 1] = -dd[i, :, 1]
for i in range(10201):
    for j in range(2):
        mean[i, j] = dd[:, i, j].mean()

meanx, meany = mean[:, 0].reshape((101, 101)), mean[:, 1].reshape((101, 101))
errx, erry = meanx-ddx_, meany-ddy_
surf = ax.plot_surface(ddx_, ddy_, errx, rstride=1, cstride=1, cmap=cm.jet)

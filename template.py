import numpy as np
import scipy as sp
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pylab

from math import pow, sqrt, exp, log, sin, cos, pi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

fig = pylab.figure()
ax = Axes3D(fig)

m = np.zeros((2594, 2774))
m[2422, 477] = 1
mask = ma.make_mask(m)
"""
# for displaying
numx, numy = 50, 100
scalex, scaley = full2theta, fullgamma
func = gaussian
x, y = np.mgrid[0:scalex:complex(0,numx), 0:scaley:complex(0,numy)]

z = np.zeros((numx, numy))
for i in range(numx):
    for j in range(numy):
        z[i, j] = func(x[i,j], y[i,j])

surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet)
"""

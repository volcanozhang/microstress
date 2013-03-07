import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pylab

from math import pow, sqrt, exp, log, sin, cos, pi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

fig = pylab.figure()
ax = Axes3D(fig)

# for displaying
numx, numy = 50, 100
scalex, scaley = full2theta, fullgamma
func = gaussian

ar = np.arange(0, scalex+scalex/numx, scalex/numx)
x = np. zeros((numx+1, numy+1))
for i in range(numy+1):
    x[:, i] = ar

ar = np.arange(0, scaley+scaley/numy, scaley/numy)
y = np. zeros((numx+1, numy+1))
for i in range(numx+1):
    y[i] = ar

z = np.zeros((numx+1, numy+1))
for i in range(numx+1):
    for j in range(numy+1):
        z[i, j] = func(x[i,j], y[i,j])

surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet)

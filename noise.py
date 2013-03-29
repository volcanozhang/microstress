import numpy as np
import numpy.ma as ma
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pylab

from math import pow, sqrt, exp, log, sin, cos, pi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

from common import stringint

fig = pylab.figure()
ax = Axes3D(fig)

offset = 4096
nimage = 100
framedim = (2594, 2774)
images = np.zeros((100, framedim[0], framedim[1]))
nb_elem = framedim[0] * framedim[1]
formatdata = np.uint16
dirpath = '/home/fengguo/Data/Si1g_5N/nomvt/'

bkg = np.load("background.npy")
for i in range(100):
    path = dirpath + 'S1gnomvt_%s_mar.tiff'%stringint(i, 4)
    f = open(path, 'rb')
    f.seek(offset)
    images[i] = np.fromfile(f, dtype = formatdata, count = nb_elem).reshape(framedim)+bkg-100
    f.close()

#np.save('21_Si1g_5N_nomvt', images)
mean, std = np.zeros(framedim), np.zeros(framedim)
for i in range(framedim[0]):
    for j in range(framedim[1]):
        mean[i, j], std[i, j] = images[:, i, j].mean(), images[:, i, j].std()

np.save('21_Si1g_5N_nomvt_mean', mean)
np.save('21_Si1g_5N_nomvt_std', std)
#mean = images[:]
"""
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
"""

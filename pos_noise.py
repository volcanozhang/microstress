import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pylab

from scipy.stats import linregress
from PeakSearch import PeakSearch
from common import stringint
from math import pow, sqrt, exp, log, sin, cos, pi
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

fig = pylab.figure()
ax = Axes3D(fig)

bkg = np.load('background.npy')

def mean_var_xy(xy, rx = 10, ry = 10):
    int_x, int_y = int(round(xy[0])), int(round(xy[1]))

    xys = []
    for i in range(int_x-rx, int_x+rx+1):
        for j in range(int_y-ry, int_y+ry+1):
            xys.append((i, j))
    xys = tuple(xys)

    subimages = np.zeros((100, len(xys)))
    dirpath = '/home/fengguo/Data/Si2g_200N/nomvt/'
    for i in range(100):
        path = dirpath + 'S2gnomvt_%s_mar.tiff'%stringint(i,4)
        f = open(path, 'rb')
        f.seek(offset)
        image = np.fromfile(f, dtype = formatdata, count = nb_elem).reshape(framedim)
        for j in range(len(xys)):
            subimages[i,j] = image[xys[j]]+bkg[xys[j]]-100
        f.close()
    mean, std = np.zeros(len(xys)), np.zeros(len(xys))
    for i in range(len(xys)):
        mean[i], std[i] = subimages[:,i].mean(), subimages[:,i].std()
    var = std**2
    return mean, var

offset = 4096
nimage = 100
framedim = (2594, 2774)
images = np.zeros((nimage, framedim[0], framedim[1]))
nb_elem = framedim[0] * framedim[1]
formatdata = np.uint16
dirpath = '/home/fengguo/Data/Si2g_200N/nomvt/'

ref = dirpath + 'S2gnomvt_0000_mar.tiff'
peaks = PeakSearch(ref,IntensityThreshold=150)[0]

rx, ry = 10, 10

xys = peaks[:, 0:2]
xykb = np.zeros((peaks.shape[0], 4))
for n in range(xys.shape[0]):#xys.shape[0]):
    pixels = []
    int_x, int_y = int(round(xys[n,1])), int(round(xys[n,0]))

    for i in range(int_x-rx, int_x+rx+1):
        for j in range(int_y-ry, int_y+ry+1):
            pixels.append((i, j))
    pixels = tuple(pixels)

    subimages = np.zeros((100, len(pixels)))
    #dirpath = '/home/fengguo/Data/Si2g_200N/nomvt/'
    for i in range(100):
        path = dirpath + 'S2gnomvt_%s_mar.tiff'%stringint(i,4)
        f = open(path, 'rb')
        f.seek(offset)
        image = np.fromfile(f, dtype = formatdata, count = nb_elem).reshape(framedim)
        for j in range(len(pixels)):
            subimages[i,j] = image[pixels[j]]+bkg[pixels[j]]-100
        f.close()
    mean, std = np.zeros(subimages.shape[1]), np.zeros(subimages.shape[1])
    for i in range(subimages.shape[1]):
        mean[i], std[i] = subimages[:,i].mean(), subimages[:,i].std()
    var = std**2
    k, b = linregress(mean, var)[0:2]
    xykb[n] = xys[n,0], xys[n,1], k, b
np.save("xykb",xykb)

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

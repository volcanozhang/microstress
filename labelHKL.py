import numpy as np
import scipy as sp
import os
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
#import pylab

from math import pow, sqrt, exp, log, sin, cos, pi
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm

#fig = pylab.figure()
#ax = Axes3D(fig)

m = np.zeros((2594, 2774))
m[2422, 477] = 1
mask = ma.make_mask(m)

name = 'peaks'

cor_apd, fit_apd = '.cor', '.fit'
path = '/home/fengguo/microstress/'
cor_path = path+'dat_'+name+cor_apd
fit_path = path+'dat_'+name+fit_apd

#ids = np.loadtxt(fname=fit_path,skiprows=5,usecols=(0))
hkls = np.loadtxt(fname=fit_path,skiprows=5,usecols=(0,2,3,4)).astype(int)
xys = np.loadtxt(fname=cor_path,skiprows=1,usecols=(2,3))

hklxy = dict()
for i in range(hkls.shape[0]):
    x, y = xys[hkls[i,0],0],xys[hkls[i,0],1]
    h, k, l = hkls[i,1],hkls[i,2],hkls[i,3]
    hklxy[(h, k, l)] = (x, y)
    plt.text(x, y, '(%i%i%i)'%(h,k,l))
offset = 4096
framedim = (2594, 2774)
nb_elem = framedim[0]*framedim[1]
formatdata = np.uint16
path = os.path.join('/', 'home', 'fengguo', 'Data', '21Feb13', 'Si1g_5N', 'nomvt', 'S1gnomvt_0000_mar.tiff')

f = open(path, 'rb')
f.seek(offset)
raw_image = np.fromfile(f, dtype = formatdata, count = nb_elem).reshape(framedim)
plt.imshow(np.log(raw_image))
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

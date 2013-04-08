import numpy as np
from numpy import random
import scipy as sp
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pylab

import libtiff
from libtiff import TIFFfile, TIFFimage

from math import pow, sqrt, exp, log, sin, cos, pi
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm

from common import xy2angle, angle2xy

D2R = pi/180.
R2D = 180./pi

#fig = pylab.figure()
#ax = Axes3D(fig)
hklfit = np.load('hklfit.npy').item()
hklxy = np.load('hkl_xy.npy').item()

hkl = hklfit.keys()[0]
alpha = hklfit[hkl][0]
x0, y0 = hklxy[hkl]

sizex, sizey = 31, 31

himage, limage = np.zeros((sizex, sizey)), np.zeros((sizex, sizey))
twoTheta0, Chi0 = xy2angle(x0, y0)
twotheta0, chi0 = twoTheta0 * D2R, Chi0 * D2R

FWHMx, FWHMy, xi, Imax = 0.003, 0.0028, 0., 1700.
full2theta, minchi, maxchi = pi, -pi, pi
sqcosxi, sin2xi, coef = pow(cos(xi), 2), sin(2*xi), 8*log(2)
sqsinxi = 1-sqcosxi
coefx, coefy = coef/pow(FWHMx, 2), coef/pow(FWHMy, 2)
A, B, C = (sqcosxi*coefx+sqsinxi*coefy)/2, sin2xi*(coefy-coefx)/4, (sqcosxi*coefy+sqsinxi*coefx)/2

gaussian = lambda x, y: Imax*exp(-(A*(x-twotheta0)**2+2*B*(x-twotheta0)*(y-chi0)+C*(y-chi0)**2))

dd = 59.836

expe, var = np.load('exp_var.npy')[0: 2]
def gen_prototype(sizex, sizey, ddh, ddl):
    def gaussian_xy_low(x, y, dd = ddl):
        cenx, ceny = angle2xy(twoTheta0, Chi0, dd = dd)
        two_Theta, Chi = xy2angle(x - sizex/2 + cenx, y - sizey/2 + ceny, dd = dd)
        two_theta, _chi = two_Theta * D2R, Chi * D2R
        return gaussian(two_theta, _chi)
    def gaussian_xy_high(x, y, dd = ddh):
        cenx, ceny = angle2xy(twoTheta0, Chi0, dd = dd)
        two_Theta, Chi = xy2angle(x - sizex/2 + cenx, y - sizey/2 + ceny, dd = dd)
        two_theta, _chi = two_Theta * D2R, Chi * D2R
        return gaussian(two_theta, _chi)
    for i in range(sizex):
        for j in range(sizey):
            lower, upper = lambda x: j, lambda x: j+1
            es_limage[i, j] = integrate.dblquad(gaussian_xy_low, i, i+1, lower, upper)[0]
            es_himage[i, j] = integrate.dblquad(gaussian_xy_high, i, i+1, lower, upper)[0]
    es_himage, es_limage = np.zeros((sizex, sizey)), np.zeros((sizex, sizey))
    np.save('es_limage', es_limage)
    np.save('es_himage', es_himage)

def SpotArray(es_l, es_h, num = 100):
    es_limage, es_himage = np.load(es_l), np.load(es_h)
    sizex, sizey = es_limage.shape
    limage, himage = np.zeros((sizex * num, sizey), np.uint16), np.zeros((sizex * num, sizey), np.uint16)
    for k in range(num):
        for i in range(sizex):
            for j in range(sizey):
                lesti, hesti = es_limage[i, j], es_himage[i, j]
                limage[sizex * k + i, j] = int(round(random.normal(lesti+expe, sqrt(alpha*lesti + var))))
                himage[sizex * k + i, j] = int(round(random.normal(hesti+expe, sqrt(alpha*hesti + var))))
    tiff = TIFFimage(limage, description = '')
    tiff.write_file('limage', compression = 'none')
    del tiff
    tiff = TIFFimage(himage, description = '')
    tiff.write_file('himage', compression = 'none')
    del tiff
    return limage, himage


# for displaying the spot
"""
es_image = np.load('es_image.npy')
images = np.zeros((100, 21, 21))
for k in range(100):
    for i in range(21):
        for j in range(21):
            esti = es_image[i, j]
            images[k, i, j] = round(random.normal(esti+expe, sqrt(alpha*esti + var)))
mean, std = np.zeros((21,21)), np.zeros((21,21))
for i in range(21):
    for j in range(21):
        mean[i, j] = images[:, i, j].mean()
        std[i, j] = images[:, i, j].std()
"""
"""
numx, numy = 50, 100
minx, maxx, miny, maxy = 0, full2theta, minchi, maxchi
func = gaussian
x, y = np.mgrid[minx:maxx:complex(0,numx),miny:maxy:complex(0,numy)]

z = np.zeros((numx, numy))
for i in range(numx):
    for j in range(numy):
        z[i, j] = func(x[i,j], y[i,j])

surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet)
"""

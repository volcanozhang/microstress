import numpy as np
from numpy import sin, cos, pi
import scipy as sp
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pylab

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

from common import stringint

FIT = True

fig = pylab.figure()
ax = Axes3D(fig)
ddx_, ddy_ = np.mgrid[0:1:complex(0,101), 0:1:complex(0,101)]

num = 100
dd = np.zeros((num, 10201, 2))
mean = np.zeros((10201, 2))
std = np.zeros((10201, 2))
for i in range(num):
    name = 'XYfin%i.pts'%i
    ini = np.loadtxt('101x101.pts', skiprows = 1)[:, 2:4]
    res = np.loadtxt(name, skiprows = 1)[:, 2:4]
    dd[i] = res - ini
    dd[i, :, 1] = -dd[i, :, 1]

for i in range(10201):
    for j in range(2):
        mean[i, j] = dd[:, i, j].mean()
        std[i,j] = dd[:, i, j].std()

meanx, meany, stdx, stdy = mean[:, 0].reshape((101, 101)), mean[:, 1].reshape((101, 101)), std[:, 0].reshape((101, 101)), std[:, 1].reshape((101, 101))
errx, erry = meanx-ddx_, meany-ddy_
varx, vary = stdx**2, stdy **2

#surf = ax.plot_surface(ddx_, ddy_, errx, rstride=1, cstride=1, cmap=cm.jet)
ax.set_xlabel('$x$ displacement',size='xx-large')
ax.set_ylabel('$y$ displacement',size='xx-large')
ax.set_zlabel('$x$ displacement error',size='xx-large')

#pylab.savefig('error_pm.pdf',format='pdf', bbox_inches='tight', pad_inches=0)
def residuals(AB, err):
    A, B = AB
    surf = A*sin(2*pi*ddx_)-B*sin(2*pi*ddy_)
    return (err-surf).ravel()
def residuals1(A, err):
    surf = A*sin(4*pi*ddx_)
    return (err-surf).ravel()


err = errx
if FIT:
    #x, y = np.mgrid[0:1:complex(0,101), 0:1:complex(0,101)]
    #sur = sin(2*pi*x)-0.1*sin(2*pi*y)
    #surf = ax.plot_surface(ddx_, ddy_, sur, rstride=1, cstride=1, cmap=cm.jet)
    A0, B0 = err[:,50].max(), err[0,:].max()
    A, B = leastsq(residuals, (A0, B0), args = err)[0]
    s = A*sin(2*pi*ddx_)-B*sin(2*pi*ddy_)
    #surf1 = ax.plot_surface(ddx_, ddy_, s, rstride=1, cstride=1, cmap=cm.jet)
    #r = errx-s
    #A0 = r.max()
    #A = leastsq(residuals1,A0,args=r)[0]
    #s = A*sin(4*pi*ddx_)

surf1 = ax.plot_surface(ddx_, ddy_, errx-s, rstride=1, cstride=1, cmap=cm.jet)
pylab.show()
"""
dd = np.zeros((10201, 2))
name = '38000id'
ini = np.loadtxt('101x101.pts', skiprows = 1)[:, 2:4]
res = np.loadtxt(name, skiprows = 1)[:, 2:4]
dd = res - ini
dd[:, 1] = -dd[:, 1]
ddx, ddy = dd[:, 0].reshape((101, 101)), dd[:, 1].reshape((101, 101))
errx, erry = ddx - ddx_, ddy - ddy_
surf = ax.plot_surface(ddx_, ddy_, errx, rstride=1, cstride=1, cmap=cm.jet)
ax.set_xlabel('$x$ displacement',size='xx-large')
ax.set_xmargin(0.5)
ax.set_ylabel('$y$ displacement',size='xx-large')
ax.set_zlabel('$x$ displacement error',size='xx-large')
pylab.show()
"""

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

FWHMx, FWHMy, chi, Imax = 2., 1., 0., 10.
twotheta, gamma, full2theta, fullgamma = pi/2, pi, pi, 2*pi
sqcoschi, sin2chi, coef = pow(cos(chi), 2), sin(2*chi), 8*log(2)
sqsinchi = 1-sqcoschi
coefx, coefy = coef/pow(FWHMx, 2), coef/pow(FWHMy, 2)
A, B, C = (sqcoschi*coefx+sqsinchi*coefy)/2, sin2chi*(coefy-coefx)/4, (sqcoschi*coefy+sqsinchi*coefx)/2

gaussian = lambda x, y: Imax*exp(-(A*pow(x-twotheta, 2)+2*B*(x-twotheta)*(y-gamma)+C*pow(y-gamma, 2)))

# for displaying the spot
numx, numy = 50, 100
scalex, scaley = full2theta, fullgamma
func = gaussian
x, y = np.mgrid[0:scalex:complex(0,numx),0:scaley:complex(0,numy)]

z = np.zeros((numx, numy))
for i in range(numx):
    for j in range(numy):
        z[i, j] = func(x[i,j], y[i,j])

surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet)

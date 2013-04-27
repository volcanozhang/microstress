import numpy as np
import matplotlib.pyplot as plt

from common import stringint
from scipy.stats import linregress
from numpy import linalg

dirpath = '/home/fengguo/Data/Si1g_5N/tdet/'
ptspath = dirpath + 'XY.pts'
stppath = dirpath + '150.stp'

pts = np.loadtxt(ptspath, skiprows = 1)[:, 2:4]
stp = np.loadtxt(stppath, skiprows = 1)[:, 2:4]

num = pts.shape[0]
#stps = stp.shape[0]/num

#stp = stp.reshape((stps, num, 2))

#plt.plot(stp[:,0,1], stp[:,0,0])
#plt.show()
#sx = linregress(stp[0,:,0], stp[1,:,0])[0]
#xc = res/(1-s)
#sy = linregress(stp[0,:,1], stp[1,:,1])[0]
#sx1 = linregress(stp[1,:,0], stp[2,:,0])[0]
x0, y0, x1, y1 = pts[:,0], pts[:,1], stp[:,0], stp[:,1]

my, mx = np.zeros((y0.shape[0])), np.zeros((x0.shape[0]))
for i in range(my.shape[0]):
    my[i] = (x0[i] * y1[i] - x1[i] * y0[i]) / (x0[i] - x1[i])
    mx[i] = -(y0[i] - y1[i]) / (x0[i] - x1[i])
xc, yc = linregress(mx, my)[0:2]

plt.axis([0,2800,0,2800])
for i in range(stp.shape[0]):
    k, b = linregress(np.array([x0[i],x1[i]]), np.array([-y0[i],-y1[i]]))[0:2]
    t = np.arange(0,2800,200)
    s = k*t+b
    plt.plot(t,s)
    plt.plot(np.array([x0[i],x1[i]]), np.array([-y0[i],-y1[i]]),'r.')
plt.plot(xc, -yc, 'ro')
plt.show()
"""
k, b = linregress(np.array([x0[0],x1[0]]), np.array([-y0[0],-y1[0]]))[0:2]
t = np.arange(0,2800,200)
s = k*t+b
plt.plot(t,s)
plt.plot(np.array([x0[0],x1[0]]), np.array([-y0[0],-y1[0]]),'r.')
"""
"""
stps = np.zeros((13, 39, 2))
for i in range(1,14):
    stppath = dirpath + stringint(i*10, 4) + '.stp'
    stps[i-1] = np.loadtxt(stppath, skiprows = 1)[:, 2:4]
"""

import numpy as np
import scipy as sp
from scipy.stats import linregress
import os
#from scipy import optimize, stats
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

hklxy, path = np.load('hklxy_peak_20.npy')

offset = 4096
framedim = (2594, 2774)
nb_elem = framedim[0]*framedim[1]
formatdata = np.uint16

bkg = np.load("background.npy")
print hklxy.keys()
def stringint(k, n):
    strint = str(k)
    res = '0' * (n - len(strint)) + strint
    return res

def show(hkl, rx = 10, ry = 10):
    f = open(path, 'rb')
    f.seek(offset)
    image = np.fromfile(f, dtype = formatdata, count = nb_elem).reshape(framedim)+bkg-100
    f.close()
    cen = hklxy[hkl]
    int_x, int_y = int(round(cen[1])), int(round(cen[0]))
    subimage = image[int_x-rx:int_x+rx+1, int_y-ry:int_y+ry+1]
    plt.imshow(np.log(subimage))
    plt.show()
    return subimage

def mean_var_hkl(hkl, rx = 10, ry = 10):
    cen = hklxy[hkl]
    int_x, int_y = int(round(cen[1])), int(round(cen[0]))

    xys = []
    for i in range(int_x-rx, int_x+rx+1):
        for j in range(int_y-ry, int_y+ry+1):
            xys.append((i, j))
    xys = tuple(xys)

    subimages = np.zeros((100, len(xys)))
    dirpath = '/home/fengguo/Data/Si1g_5N/nomvt/'
    for i in range(100):
        path = dirpath + 'S1gnomvt_%s_mar.tiff'%stringint(i,4)
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
def mean_var_xy(xy, rx = 10, ry = 10):
    int_x, int_y = int(round(xy[0])), int(round(xy[1]))

    xys = []
    for i in range(int_x-rx, int_x+rx+1):
        for j in range(int_y-ry, int_y+ry+1):
            xys.append((i, j))
    xys = tuple(xys)

    subimages = np.zeros((100, len(xys)))
    dirpath = '/home/fengguo/Data/Si1g_5N/nomvt/'
    for i in range(100):
        path = dirpath + 'S1gnomvt_%s_mar.tiff'%stringint(i,4)
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
"""
mean, var = mean_var_hkl((-5,-5,-7),rx=7,ry=7)
k, b = linregress(mean, var)[0:2]
fig, ax = plt.subplots(1)
plt.axis([0, mean.max()+1, 0, var.max()+1])
plt.plot(mean, var, '.')
t = np.arange(0, mean.max(), mean.max()/10)
s = k * t + b
ax.set_xlabel('mean')
ax.set_ylabel('variance')
plt.plot(t, s, 'r')
textstr = 'slope$=\\alpha$\ninterception with X axis$=\mu-\sigma^2/\\alpha$'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
"""
"""
res = []
for key in hklxy.keys():
    mean, var = mean_var_hkl(key, rx=7, ry=7)
    k, b = linregress(mean, var)[0:2]
    res.append([1/k, -b/k])
res = np.array(res)
fig, ax = plt.subplots(1)
ax.set_xlabel('$\\frac{1}{\\alpha}$')
ax.set_ylabel('$\mu-\\frac{\sigma^2}{\\alpha}$')
plt.plot(res[:,0], res[:,1], '.')
"""
res = []
peaks = np.load('peaks.npy')[2]
for i in range(peaks.shape[0]):
    mean, var = mean_var_xy((peaks[i, 1], peaks[i, 0]), rx=7, ry=7)
    k, b = linregress(mean, var)[0:2]
    res.append([1/k, -b/k])
res = np.array(res)
fig, ax = plt.subplots(1)
ax.set_xlabel('$\\frac{1}{\\alpha}$')
ax.set_ylabel('$\mu-\\frac{\sigma^2}{\\alpha}$')
plt.plot(res[:,0], res[:,1], '.')
k, b = linregress(res[:,0], res[:,1])[0:2]
t = np.arange(res[:,0].min()-1, res[:,0].max()+1, (res[:,0].max()-res[:,0].min()+2)/10)
s = k * t + b
plt.plot(t, s, 'r')
plt.show()

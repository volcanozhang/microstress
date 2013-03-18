import numpy as np
import os
#from scipy import optimize, stats
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

hklxy, path = np.load('hklxy_peak_20.npy')

offset = 4096
framedim = (2594, 2774)
nb_elem = framedim[0]*framedim[1]
formatdata = np.uint16


print hklxy.keys()
def stringint(k, n):
    strint = str(k)
    res = '0' * (n - len(strint)) + strint
    return res

def show(hkl, rx = 10, ry = 10):
    f = open(path, 'rb')
    f.seek(offset)
    image = np.fromfile(f, dtype = formatdata, count = nb_elem).reshape(framedim)
    f.close()
    cen = hklxy[hkl]
    int_x, int_y = int(round(cen[1])), int(round(cen[0]))
    subimage = image[int_x-rx:int_x+rx+1, int_y-ry:int_y+ry+1]
    plt.imshow(np.log(subimage))
    plt.show()
    return subimage

def mean_std(hkl, rx = 10, ry = 10):
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
            subimages[i,j] = image[xys[j]]
        f.close()
    mean, std = np.zeros(len(xys)), np.zeros(len(xys))
    for i in range(len(xys)):
        mean[i], std[i] = subimages[:,i].mean(), subimages[:,i].std()
    return mean, std

mean, std = mean_std((-2,-2,-4),rx=7,ry=7)

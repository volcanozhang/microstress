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

f = open(path, 'rb')
f.seek(offset)
image = np.fromfile(f, dtype = formatdata, count = nb_elem).reshape(framedim)
f.close()

print hklxy.keys()
def show(hkl, rx = 10, ry = 10):
    cen = hklxy[hkl]
    int_x, int_y = int(round(cen[1])), int(round(cen[0]))
    subimage = image[int_x-rx:int_x+rx+1, int_y-ry:int_y+ry+1]
    plt.imshow(np.log(subimage))
    plt.show()
    return subimage

def xys(hkl, rx = 10, ry = 10):
    cen = hklxy[hkl]
    int_x, int_y = int(round(cen[1])), int(round(cen[0]))
    subimage = image[int_x-rx:int_x+rx+1, int_y-ry:int_y+ry+1]
    xys = []
    for i in range(int_x-rx, int_x+rx+1):
        for j in range(int_y-ry, int_y+ry+1):
            xys.append((i, j))
    return tuple(xys)

import numpy as np
#from scipy import optimize, stats
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import scipy as sp
from numpy import random
from scipy import interpolate

import libtiff
from libtiff import TIFFfile, TIFFimage

def image_shifted(mean, std, sx = 21, sy = 21, nx = 101, ny = 101, kx=1, ky=1, initial = False, noisy = True, rep = 100):
    mx, my = mean.shape
    if mx<sx+1 or my<sy+1:
        return 0
    image = np.zeros((sx*nx,sy*ny),dtype=np.uint16)
    fimage = np.zeros((sx*nx,sy*ny))
    fvar = np.zeros((sx*nx,sy*ny))
    if initial:
       for m in range(nx):
           for n in range(ny):
               image[m*sx:(m+1)*sx, n*sy:(n+1)*sy] = np.int_(np.round_(mean[0:sx,0:sy]))
       tiff = TIFFimage(image, description = '')
       tiff.write_file('image_ini', compression = 'none')
       del tiff
       return 0

    accum=np.zeros((sx+2, sy+2))
    #sub = np.zeros((sx+1, sy+1))

    for m in range(nx):
        for n in range(ny):
            for i in range(sx+2):
                for j in range(sy+2):
                    accum[i,j] = mean[0:i,0:j].sum()

            spline = interpolate.RectBivariateSpline(np.arange(0,sx+2),np.arange(0,sy+2),accum,kx=kx,ky=ky,s=0)
            subimage = fimage[m*sx:(m+1)*sx,n*sy:(n+1)*sy]
            for i in range(sx):
                for j in range(sy):
                    subimage[i,j] = spline.ev(m/(nx-1.)+i+1,n/(ny-1.)+j+1)-spline.ev(m/(nx-1.)+i,n/(ny-1.)+j+1)-spline.ev(m/(nx-1.)+i+1,n/(ny-1.)+j)+spline.ev(m/(nx-1.)+i,n/(ny-1.)+j)
    var = std**2
    if noisy:
        for m in range(nx):
            for n in range(ny):
                for i in range(sx+2):
                    for j in range(sy+2):
                        accum[i,j] = var[0:i,0:j].sum()
                spline = interpolate.RectBivariateSpline(np.arange(0,sx+2),np.arange(0,sy+2),accum,kx=kx,ky=ky,s=0)
                subvar = fvar[m*sx:(m+1)*sx,n*sy:(n+1)*sy]
                for i in range(sx):
                    for j in range(sy):
                        subvar[i,j] = spline.ev(m/(nx-1.)+i+1,n/(ny-1.)+j+1)-spline.ev(m/(nx-1.)+i,n/(ny-1.)+j+1)-spline.ev(m/(nx-1.)+i+1,n/(ny-1.)+j)+spline.ev(m/(nx-1.)+i,n/(ny-1.)+j)
        for k in range(100):
            for i in range(sx*nx):
                for j in range(sy*ny):
                    image[i,j] = int(round(random.normal(fimage[i,j],np.sqrt(fvar[i,j]))))
            tiff = TIFFimage(image, description = '')
            tiff.write_file('image%i'%k, compression = 'none')
            del tiff
    else:
        for i in range(sx*nx):
            for j in range(sy*ny):
                image[i,j] = int(round(fimage[i,j]))
        tiff = TIFFimage(image, description = '')
        tiff.write_file('image', compression = 'none')
        del tiff
    return 0
i=0
mean, std = np.load('mean%i.npy'%i)[5:,5:], np.load('std%i.npy'%i)[5:,5:]
image = image_shifted(mean, std,noisy=True)

import numpy as np
import scipy as sp
import os
#from scipy import optimize, stats
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.cm as cm
from numpy import pi
from scipy import cos, sin

from scipy.fftpack import fft2, ifft2, dct

ind = np.load('individual.npy')

def fft_shift(image, ux, uy):
    fimage = fft2(ind)
    nx, ny = image.shape
    shift = np.zeros((nx,ny), dtype=np.complex)
    for i in range(nx):
        for j in range(ny):
            pha = 2*pi*(i*ux+j*uy)/nx
            shift[i, j] = complex(cos(pha), sin(pha))
    shifted = ifft2(fimage * shift)
    shifted = np.int_(np.real(shifted))
    return shifted
shifted = fft_shift(ind,-10,0)
#plt.imshow(ind)
#plt.show()

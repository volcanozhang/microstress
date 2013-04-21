planck = 4.135667516e-18
lightspeed = 299792458
C = planck * lightspeed / 1e-10
import numpy as np
from numpy import sqrt, cos, sin, arccos, arcsin

from common import Euler, multi_hkl, PI, D2R, R2D

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

def hkl_diamond(emin = 8, emax = 22, a = 3.5):
    #_wlmax = C / emin
    #wlmax = min(_wlmax, 2 * a)
    wlmin, wlmax = C / emax, C / emin

    hkllim = int(round(2 * a / wlmin))+1
    hklkey = set()
    for h in range(hkllim):
        for k in range(hkllim):
            for l in range(hkllim):
                if (h%2==0 and k%2==0 and l%2==0) or (h%2 and k%2 and l%2):
                    if (h+k+l)%4 != 2:
                        nh, nk, nl = sorted([h,k,l])
                        hklkey.add(nh*hkllim*hkllim+nk*hkllim+nl)
    hklkey.discard(0)
    hklhar = []
    for n in hklkey:
        h = n / (hkllim * hkllim)
        k = (n-h*hkllim*hkllim) / hkllim
        l = n % hkllim
        hklhar.append([h,k,l])
    #hkl = np.array(hkl)
#   remove harmonic parts.
    #hkl_har = hkl.tolist()
    hklnohar = []
    for hkl0 in hklhar:
        hklnohar.append(hkl0)
        for hkl1 in hklhar:
            if hkl0[0]*hkl1[1]-hkl0[1]*hkl1[0]==0 and hkl0[0]*hkl1[2]-hkl0[2]*hkl1[0]==0 and hkl0[2]*hkl1[1]-hkl0[1]*hkl1[2]==0:
                hklhar.remove(hkl1)
    hklnohar = np.array(hklnohar)
    multihklnohar = []
    for i in range(hklnohar.shape[0]):
        multihklnohar = multihklnohar + multi_hkl(hklnohar[i]).tolist()
    return np.array(multihklnohar)

def SimCCD(hklnohar, Angs = (0,0,0), dd = 59.836, xcen = 1365.74, ycen = 943.05, xbet = 0.373, xgam = 0.504, pix = 0.031, sx = 2774, sy = 2594):
    varphi1, phi0, varphi2 = Angs
    mat = Euler(varphi1, phi0, varphi2)
    hklnohar_ = hklnohar.copy()

    vecs = np.inner(hklnohar_, mat)
    xcams, ycams, ids = [], [], []
    for i in range(vecs.shape[0]):
        vec = vecs[i]
        norm = sqrt((vec ** 2).sum())
        vecs[i] = vec/norm

        uflab = vecs[i].T

        cosbeta = cos(PI / 2. - xbet * D2R)
        sinbeta = sin(PI / 2. - xbet * D2R)

        IOlab = dd * np.array([0.0, cosbeta, sinbeta])
        unlab = IOlab / sqrt(np.dot(IOlab, IOlab))
        norme_uflab = sqrt(np.sum(uflab ** 2))
        uflab = uflab / norme_uflab

        scal = np.dot(uflab, unlab)
        if scal == 0:
            continue
        if scal < 0:
            uflab = -uflab
            hklnohar_[i] = -hklnohar_[i]
            scal = -scal
        normeIMlab = dd / scal

        IMlab = uflab * normeIMlab
        OMlab = IMlab - IOlab
    
        xca0 = OMlab[0]
        yca0 = OMlab[1] / sinbeta

        cosgam = cos(-xgam * D2R)
        singam = sin(-xgam * D2R)

        xcam1 = cosgam * xca0 + singam * yca0
        ycam1 = -singam * xca0 + cosgam * yca0

        xcam = xcen + xcam1 / pix
        ycam = ycen + ycam1 /pix

        if xcam < sx and xcam >0 and ycam < sy and ycam >0:
            xcams.append(xcam)
            ycams.append(ycam)
            ids.append(i)
    plt.plot(xcams, ycams, '.')
    for i in range(len(ids)):
        #print hklnohar[i]
        k = ids[i]
        plt.text(xcams[i], ycams[i], '%i %i %i'%(hklnohar_[k,0],hklnohar_[k,1],hklnohar_[k,2]))
    plt.axis([0, sx, 0, sy])
    plt.show()
    return hklnohar
t=SimCCD(hkl_diamond())

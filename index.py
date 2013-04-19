planck = 4.135667516e-18
lightspeed = 299792458
C = planck * lightspeed / 1e-10
import numpy as np

from common import Euler

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
        multihklnohar = multihklnohar + multi_hkl(hklnohar[i])
    return hklnohar
"""
def SimCCD(hklnobar, Euler = (1,2,3), dd, dd = 59.836, xcen = 1365.74, ycen = 943.05, xbet = 0.373, xgam = 0.504, pix = 0.031, sx = 2774, sy = 2594):
    varphi1, phi0, varphi2 = euler
    mat = Euler(varphi1, phi0, varphi2)
"""
    

import numpy as np
from numpy import sin, cos, arccos, sqrt

PI = np.pi
D2R = PI/180.
R2D = 180./PI
alpha = 40.# alpha is treated as the incident angle.
u_i = np.array([0., 1., 0.])# direction of incident angle.
pix = 0.031

def stringint(k, n):
    strint = str(k)
    res = '0' * (n - len(strint)) + strint
    return res

def angle2xy(twicetheta, chi, dd = 59.836, xcen = 1365.74, ycen = 943.05, xbet = 0.373, xgam = 0.504):
    ctw = cos(twicetheta * D2R)
    stw = sin(twicetheta * D2R)
    cchi = cos(chi * D2R)
    schi = sin(chi * D2R)

    xuflab = -stw * schi
    yuflab = ctw
    zuflab = stw * cchi

    uflab = np.array([xuflab, yuflab, zuflab]).T

    cosbeta = cos(PI / 2. - xbet * D2R)
    sinbeta = sin(PI / 2. - xbet * D2R)

    IOlab = dd * np.array([0.0, cosbeta, sinbeta])
    unlab = IOlab / sqrt(np.dot(IOlab, IOlab))
    norme_uflab = sqrt(np.sum(uflab ** 2))
    uflab = uflab / norme_uflab

    scal = np.dot(uflab, unlab)
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

    return xcam, ycam

def xy2angle(xcam, ycam, dd = 59.836, xcen = 1365.74, ycen = 943.05, xbet = 0.373, xgam = 0.504):
    cosbeta = cos(PI / 2. - xbet * D2R)
    sinbeta = sin(PI / 2. - xbet * D2R)

    cosgam = cos(-xgam * D2R)
    singam = sin(-xgam * D2R)

    xcam1 = (xcam - xcen) * pix
    ycam1 = (ycam - ycen) * pix

    xca0 = cosgam * xcam1 - singam * ycam1
    yca0 = singam * xcam1 + cosgam * ycam1

    xO, yO, zO = dd * np.array([0.0, cosbeta, sinbeta])

    xOM = xca0
    yOM = yca0 * sinbeta
    zOM = -yca0 * cosbeta

    xM = xO + xOM
    yM = yO + yOM
    zM = zO + zOM

    IMlab = np.array([xM, yM, zM]).T

    nIMlab = 1.*sqrt(xM ** 2 + yM ** 2 + zM ** 2)
    uflab = IMlab / nIMlab

    twicetheta = arccos(uflab[1])
    coschi = uflab[2] / sin(twicetheta)
    twicetheta = twicetheta * R2D
    chi = arccos(coschi) * R2D
    if uflab[0] > 0:
        chi = -chi
    return twicetheta, chi

def Euler(ang1, ang2, ang3):
    c1, c2, c3, s1, s2, s3 = cos(ang1), cos(ang2), cos(ang3), sin(ang1), sin(ang2), sin(ang3)

    mat = np.array([[c1 * c3 - c2 * s1 * s3, c3 * s1 + c1 * c2 * s3, s2 * s3],
                    [-c1 * s3 - c2 * c3 * s1, c1 * c2 * c3 - s1 * s3, c3 * s2],
                    [s1 * s2, -c1 * s2, c2]
                    ])
    return mat

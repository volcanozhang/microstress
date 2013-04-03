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

def angle2xy(two_Theta, Chi, dd = 59.836, Xcen = 1365.74, Ycen = 943.05, Beta = 0.373, Gamma = 0.504):
    two_theta, chi, beta, gamma = np.array([two_Theta, Chi, Beta, Gamma]) * D2R
    DD = dd / pix
    
    sin2theta, sinchi, sinbeta, singamma = sin(np.array([two_theta, chi, beta, gamma]))
    cos2theta, coschi, cosbeta, cosgamma = cos(np.array([two_theta, chi, beta, gamma]))

    uf1, uf2, uf3 = -sin2theta * sinchi, cos2theta, sin2theta * coschi
    xL, yL = np.array([uf1, uf2]) * DD / uf3
    return cosgamma * xL - singamma * yL /cosbeta + Xcen, singamma * xL + cosgamma * yL /cosbeta + Ycen

def xy2angle(X, Y, dd = 59.836, Xcen = 1365.74, Ycen = 943.05, Beta = 0.373, Gamma = 0.504):
    beta, gamma = np.array([Beta, Gamma]) * D2R
    DD = dd / pix
    
    sinbeta, singamma = sin(np.array([beta, gamma]))
    cosbeta, cosgamma = cos(np.array([beta, gamma]))

    Xr, Yr = X - Xcen, Y - Ycen
    f1, f2, f3 = cosgamma * Xr + singamma * Yr, (-singamma * Xr + cosgamma * Yr) * cosbeta, DD

    norm = sqrt(f1 ** 2 + f2 ** 2 + f3 **2)
    two_theta = arccos(f2/norm)
    sin2theta = sin(two_theta)
    sinchi, coschi = -f1/(norm*sin2theta), f3/(norm*sin2theta)
    chi = arccos(coschi)
    if sinchi < 0:
        chi = -chi
    return two_theta * R2D, chi * R2D

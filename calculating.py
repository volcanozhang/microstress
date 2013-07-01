dd0, xcen0, ycen0, bet0, gam0 = 59.832, 1365.71, 942.85, 0.377, 0.505
dd1, xcen1, ycen1, bet1, gam1 = 59.832, 1365.71, 942.85, 0.377, 0.505
import numpy as np
import scipy as sp
from numpy import sin, cos, pi
from scipy import linalg
from scipy.optimize import leastsq
from common import point_D2L, Normalize
import matplotlib.pyplot as plt
DEG2RAD = np.pi/180.0
k_i = np.array([0.0, 1.0, 0.0])
pix_size = 0.03100
D2R = pi/180.
R2D = 180./pi
base_strain = np.array(
                      [[1., 0., 0., 0. ,0. ,0. ,0. ,0. ,-1.],
                       [0., 1., 0., 0. ,0. ,0. ,0. ,0. , 0.],
                       [0., 0., 1., 0. ,0. ,0. ,0. ,0. , 0.],
                       [0., 0., 0., 1. ,0. ,0. ,0. ,0. , 0.],
                       [0., 0., 0., 0. ,1. ,0. ,0. ,0. ,-1.],
                       [0., 0., 0., 0. ,0. ,1. ,0. ,0. , 0.],
                       [0., 0., 0., 0. ,0. ,0. ,1. ,0. , 0.],
                       [0., 0., 0., 0. ,0. ,0. ,0. ,1. , 0.]])
base_strain = base_strain.reshape((8,3,3))

def residuals(Fres, q0, q1):
    F = np.ones(9)
    F[1:9] = Fres
    num = q0.shape[0]
    res = np.zeros(num)
    for i in range(num):
        q1_ = Normalize(np.tensordot(F.reshape(3,3),q0[i],(1,0)))
        res[i] = 1 - np.dot(q1[i],q1_)
    return res

def Get_F(xys0_pix, xys1_pix, dd0, xcen0, ycen0, bet0, gam0, dd1, xcen1, ycen1, bet1, gam1, pix_size = 0.031, Iter = True):
    num = xys0_pix.shape[0]
    #k0s, k1s = np.zeros((num,3)), np.zeros((num,3))
    A, b = np.zeros((num*3,8+num)), np.zeros(num*3)
    q0, q1 = np.zeros((num,3)), np.zeros((num,3))
    for i in range(num):
        x0, y0 = xys0_pix[i]
        x1, y1 = xys1_pix[i]
        p0, p1 = point_D2L((x0-xcen0, y0-ycen0, 0), bet0, gam0, dd0/pix_size), point_D2L((x1-xcen1, y1-ycen1, 0), bet1, gam1, dd1/pix_size)
        p0, p1 = np.array(p0), np.array(p1)
        #k0s[i], k1s[i] = Normalize(p0), Normalize(p1)
        q0[i], q1[i] = Normalize(Normalize(p0)-k_i), Normalize(Normalize(p1)-k_i)
        A[3*i,0:8] = q0[i,1], q0[i,2], 0, 0, 0, 0, 0, 0
        A[3*i+1,0:8] = 0, 0, q0[i,0], q0[i,1], q0[i,2], 0, 0, 0
        A[3*i+2,0:8] = 0, 0, 0, 0, 0, q0[i,0], q0[i,1], q0[i,2]
        A[3*i,8+i] = -q1[i,0]
        A[3*i+1,8+i] = -q1[i,1]
        A[3*i+2,8+i] = -q1[i,2]
        b[3*i:3*(i+1)] = -q0[i,0], 0, 0
    F = np.ones(9)
    #F[0] = 1
    F[1:9] = linalg.solve(np.dot(A.T,A), np.dot(A.T,b))[0:8]
    #F = F.reshape(3,3)
    if Iter:
        F_ = np.ones(9)
        F_[1:9], msg = leastsq(residuals, F[1:9], args=(q0,q1))
        if msg == 1 or msg == 2 or msg == 3 or msg == 4:
            F_ = F_.reshape(3,3)
            det = linalg.det(F_)
            return F_/det**(1./3)
        else:
            F = F.reshape(3,3)
            det = linalg.det(F)
            return F/det**(1./3)
    else:
        F = F.reshape(3,3)
        det = linalg.det(F)
        return F/det**(1./3)

solve = np.zeros((285,3,3))
for i in range(609, 894):
    xy0, xy1 = np.loadtxt('/home/fengguo/microstress/%iref'%i, skiprows=1)[:,2:4], np.loadtxt('/home/fengguo/microstress/%i200'%i, skiprows=1)[:,2:4]
    xys0, xys1 = np.zeros(xy0.shape), np.zeros(xy1.shape)
    xys0[:,0], xys1[:,0] = -xy0[:,1], -xy1[:,1]
    xys0[:,1], xys1[:,1] = xy0[:,0], xy1[:,0]
    solve[i-609]=Get_F(xys0, xys1, dd0, xcen0, ycen0, bet0, gam0, dd1, xcen1, ycen1, bet1, gam1, pix_size = 0.031)#,Iter=False)
np.save('solve070113',solve)
#F=Get_F(xys0, xys1, dd0, xcen0, ycen0, bet0, gam0, dd1, xcen1, ycen1, bet1, gam1, pix_size = 0.031)

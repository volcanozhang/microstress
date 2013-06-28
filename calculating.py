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
k_i_u = np.array([0.0, 1.0, 0.0])
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

def residuals(solve, k_f_u0, k_f_u1):
    num = k_f_u0.shape[0]
    #print solve
    strain = np.zeros((3,3))
    for i in range(0,8):
        strain = strain+solve[i]*base_strain[i]
    res = np.zeros(num)
    for i in range(num):
        temp = k_f_u0[i]-k_i_u
        temp = temp+np.dot(strain, temp)
        q = Normalize(temp)
        cos = -np.dot(q, k_i_u)
        k_f_ue = k_i_u + 2*cos*q
        res[i] = 1 - np.dot(k_f_ue,k_f_u1[i])
    return res
"""
def Get_LS_DGradient(xys0_pix, xys1_pix, dd0, xcen0, ycen0, bet0, gam0, dd1, xcen1, ycen1, bet1, gam1, pix_size = 0.031, eta=0.5, Iter = True):
    dd0, bet0, gam0 = dd0/pix_size, bet0*D2R, gam0*D2R
    dd1, bet1, gam1 = dd1/pix_size, bet1*D2R, gam1*D2R
    cosbet0, sinbet0, cosgam0, singam0 = cos(bet0), sin(bet0), cos(gam0), sin(gam0)
    cosbet1, sinbet1, cosgam1, singam1 = cos(bet1), sin(bet1), cos(gam1), sin(gam1)
    num = len(xys0_pix)
    xys0, xys1 = xys0_pix*pix_size, xys1_pix*pix_size
    dxy = (xys1-xys0).reshape(2*num)
    ed_0, ed_1, ed_2 = np.array([cosgam,-singam*cosbet,singam*sinbet]), np.array([singam,cosgam*cosbet,-cosgam*sinbet]), np.array([0,sinbet,cosbet])
    omegax, omegay = np.outer(ed_0,ed_2)-np.outer(ed_2,ed_0), np.outer(ed_1,ed_2)-np.outer(ed_2,ed_1)
    coef = np.zeros((8, 2*num))
    k_f_u0, k_f_u1 = np.zeros((num, 3)), np.zeros((num, 3))
    for i in range(0, num):
        p_x0, p_y0, k_f_u0[i]= Get_P_Tensor(xys0_pix[i])
        p_x1, p_y1, k_f_u1[i]= Get_P_Tensor(xys1_pix[i])
        p_x, p_y = p_x0*eta+p_x1*(1-eta), p_y0*eta+p_y1*(1-eta)
        #p_x, p_y = Get_P_Tensor(xys0_pix[i])
        xcoef, ycoef = coef[2*i], coef[2*i+1]
    #base = base_strain[i]
        for j in range(0, 8):
            base = base_strain[j]
            xcoef[j], ycoef[j] = np.tensordot(base, p_x), np.tensordot(base, p_y)
    #coef = np.matrix(coef)
    #return coef, dxy
    ls_coef, ls_dxy = np.dot(coef.T, coef), np.dot(coef.T, dxy)
    #solve = np.dot(ls_coef.I, ls_dxy)
    solve = linalg.solve(ls_coef, ls_dxy)
    #print dxy
    #strain = np.zeros((3,3))
    #return solve
    if Iter == False:
        strain = np.tensordot(solve, base_strain, [0,0])
        return strain
    else:
        solve_, msg = leastsq(residuals, solve, args=(k_f_u0,k_f_u1))
        if msg == 1 or msg == 2 or msg == 3 or msg == 4:
            strain = np.tensordot(solve_, base_strain, [0,0])
            return strain
        else:
            strain = np.tensordot(solve, base_strain, [0,0])
            return strain
"""

def Get_F(xys0_pix, xys1_pix, dd0, xcen0, ycen0, bet0, gam0, dd1, xcen1, ycen1, bet1, gam1, pix_size = 0.031):
    num = xys0_pix.shape[0]
    #k0s, k1s = np.zeros((num,3)), np.zeros((num,3))
    A, b = np.zeros((num*3,8+num)), np.zeros(num*3)
    for i in range(num):
        x0, y0 = xys0_pix[i]
        x1, y1 = xys1_pix[i]
        p0, p1 = point_D2L((x0-xcen0, y0-ycen0, 0), bet0, gam0, dd0/pix_size), point_D2L((x1-xcen1, y1-ycen1, 0), bet1, gam1, dd1/pix_size)
        p0, p1 = np.array(p0), np.array(p1)
        #k0s[i], k1s[i] = Normalize(p0), Normalize(p1)
        k0, k1 = Normalize(p0), Normalize(p1)
        A[3*i,0:8] = k0[1], k0[2], 0, 0, 0, 0, 0, 0
        A[3*i+1,0:8] = 0, 0, k0[0], k0[1], k0[2], 0, 0, 0
        A[3*i+2,0:8] = 0, 0, 0, 0, 0, k0[0], k0[1], k0[2]
        A[3*i,8+i] = -k1[0]
        A[3*i+1,8+i] = -k1[1]
        A[3*i+2,8+i] = -k1[2]
        b[3*i:3*(i+1)] = -k0[0], 0, 0
    F = np.zeros(9)
    F[0] = 1
    F[1:9] = linalg.solve(np.dot(A.T,A), np.dot(A.T,b))[0:8]
    F = F.reshape(3,3)
    det = linalg.det(F)
    return F/det**(1./3)

solve = np.zeros((285,3,3))
for i in range(609, 894):
    xy0, xy1 = np.loadtxt('/home/fengguo/microstress/%iref'%i, skiprows=1)[:,2:4], np.loadtxt('/home/fengguo/microstress/%i200'%i, skiprows=1)[:,2:4]
    xys0, xys1 = np.zeros(xy0.shape), np.zeros(xy1.shape)
    xys0[:,0], xys1[:,0] = -xy0[:,1], -xy1[:,1]
    xys0[:,1], xys1[:,1] = xy0[:,0], xy1[:,0]
    solve[i-609]=Get_F(xys0, xys1, dd0, xcen0, ycen0, bet0, gam0, dd1, xcen1, ycen1, bet1, gam1, pix_size = 0.031)
#F=Get_F(xys0, xys1, dd0, xcen0, ycen0, bet0, gam0, dd1, xcen1, ycen1, bet1, gam1, pix_size = 0.031)

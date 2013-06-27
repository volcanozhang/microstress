dd, xcen, ycen, beta, gamma = 59.832, 1365.71, 942.85, 0.377, 0.505
import numpy as np
import scipy as sp
from scipy import linalg
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
DEG2RAD = np.pi/180.0
k_i_u = np.array([0.0, 1.0, 0.0])
pix_size = 0.03100
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
def Det_Lab_Mat(beta_deg = beta, gamma_deg = gamma):
    beta, gamma = beta_deg*DEG2RAD, gamma_deg * DEG2RAD
    sin_beta, cos_beta, sin_gamma, cos_gamma = np.sin(beta), np.cos(beta), np.sin(gamma), np.cos(gamma)
    beta_mat = np.array([[1.0, 0.0, 0.0], [0.0, cos_beta, sin_beta], [0.0, -sin_beta, cos_beta]])
    gamma_mat = np.array([[cos_gamma, sin_gamma, 0.0], [-sin_gamma, cos_gamma, 0.0], [0.0, 0.0, 1.0]])
    return np.dot(beta_mat, gamma_mat)
    
def Get_Omega(det2lab_mat):
    e_x, e_y, e_z = det2lab_mat
    omega_x = np.outer(e_x, e_z) - np.outer(e_z, e_x)
    omega_y = np.outer(e_y, e_z) - np.outer(e_z, e_y)
    return omega_x, omega_y
    
def Get_Length(vector):
    return np.sqrt(vector.__pow__(2).sum())

def Normalize(vector):
    length = Get_Length(vector)
    return vector/length

def Detector_to_Lab(xy_pix, det2lab_mat, dd = dd, xcen_pix = xcen, ycen_pix = ycen, pix_size = pix_size):
    x_pix, y_pix = xy_pix
    temp = np.array([0.0, 0.0, dd])
    temp[0: 2] = np.array([x_pix-xcen_pix, y_pix-ycen_pix])*pix_size
    k_f = np.dot(det2lab_mat, temp)
    k_f_u = Normalize(k_f)
    return k_f_u

def Get_P_Tensor(xy_pix, beta_deg = beta, gamma_deg = gamma, dd = dd, xcen_pix = xcen, ycen_pix = ycen):
    det2lab_mat = Det_Lab_Mat(beta_deg = beta_deg, gamma_deg = gamma_deg)
    e_x, e_y, e_z = det2lab_mat
    omega_x, omega_y = Get_Omega(det2lab_mat = det2lab_mat)
    k_f_u = Detector_to_Lab(xy_pix = xy_pix, det2lab_mat = det2lab_mat, dd = dd, xcen_pix = xcen_pix, ycen_pix = ycen_pix, pix_size = pix_size)
    sin_theta = Get_Length(k_f_u-k_i_u)/2
    q_u = (k_f_u-k_i_u)/(2*sin_theta)

    w_x, w_y = np.dot(omega_x, k_f_u), np.dot(omega_y, k_f_u)
    r_x, r_y = np.dot(q_u, w_x), np.dot(q_u, w_y)
    
    coef = dd/pow(np.dot(k_f_u, e_z), 2)
    temp1 = np.outer(k_i_u, q_u)
    temp2 = np.outer(q_u, q_u)*sin_theta
    p_x = coef*(2*sin_theta*np.outer(w_x, q_u)-2*r_x*temp1-4*r_x*temp2)
    p_y = coef*(2*sin_theta*np.outer(w_y, q_u)-2*r_y*temp1-4*r_y*temp2)
    return p_x, p_y, k_f_u

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

def Get_LS_DGradient(xys0_pix, xys1_pix, eta=0.5, Iter = True):
    num = len(xys0_pix)
    xys0, xys1 = xys0_pix*pix_size, xys1_pix*pix_size
    dxy = (xys1-xys0).reshape(2*num)
    
    coef = np.zeros((2*num, 8))
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

    #return coef, dxy
def Get_Effective_Strain(strain):
    iso = (strain[0,0]+strain[1,1]+strain[2,2])/3
    temp = strain - iso*np.eye(3)
    return np.sqrt(temp.__pow__(2).sum())*np.sqrt(2.0/3.0)

def Judge_By_Determinant(xys0_pix, xys1_pix, det2lab_mat, dd = 69.579, xcen_pix = 931.377, ycen_pix = 1063.829):
    num = len(xys0_pix)
    q0s_u, q1s_u = np.zeros((num, 3)), np.zeros((num, 3))
    #det2lab_mat = Det_Lab_Mat(beta_deg = beta_deg, gamma_deg = gamma_deg)
    for i in range(0, num):
        k_f_u = Detector_to_Lab(xy_pix = xys0_pix[i], det2lab_mat = det2lab_mat, dd = dd, xcen_pix = xcen_pix, ycen_pix = ycen_pix, pix_size = pix_size)  
        q0s_u[i] = k_f_u-k_i_u
        two_sin_theta = Get_Length(q0s_u[i])
        q0s_u[i] = q0s_u[i]/(two_sin_theta)

        k_f_u = Detector_to_Lab(xy_pix = xys1_pix[i], det2lab_mat = det2lab_mat, dd = dd, xcen_pix = xcen_pix, ycen_pix = ycen_pix, pix_size = pix_size)  
        q1s_u[i] = k_f_u-k_i_u
        two_sin_theta = Get_Length(q1s_u[i])
        q1s_u[i] = q1s_u[i]/(two_sin_theta)
    coef = np.zeros((3*num, num+9))
    for i in range(0, num):
        seg0 = coef[3*i]
        seg1 = coef[3*i+1]
        seg2 = coef[3*i+2]
        u0 = q0s_u[i]
        u1 = q1s_u[i]
        i_9 = i+9
        seg0[0:3] = u0
        seg0[i_9] = -u1[0]

        seg1[3:6] = u0
        seg1[i_9] = -u1[1]
    
        seg2[6:9] = u0
        seg2[i_9] = -u1[2]
    coef = np.mat(coef)
    return linalg.det(np.dot(coef.T, coef))

def Simulate_Shifted_Spots(xys_pix, Gradient, det2lab_mat, dd = dd, xcen_pix = xcen, ycen_pix = ycen):
    F = Gradient + np.eye(3)
    num = len(xys_pix)
    qs = np.zeros((num, 3))
    for i in range(0, num):
        k_f_u = Detector_to_Lab(xy_pix = xys_pix[i], det2lab_mat = det2lab_mat, dd = dd, xcen_pix = xcen_pix, ycen_pix = ycen_pix, pix_size = pix_size)  
        temp = k_f_u-k_i_u
        temp = np.dot(F, temp)
        qs[i] = Normalize(temp)
    xys_ = np.zeros((num, 2))
    for i in range(0, num):
        q = qs[i]
        cos = -np.dot(q, k_i_u)
        k_f_u = k_i_u + 2*cos*q
        k_f_u_t = np.dot(det2lab_mat.T, k_f_u)
        factor = dd/k_f_u_t[2]
        xys_[i] = factor*k_f_u_t[0:2]/pix_size
    return xys_+np.array([xcen_pix, ycen_pix])
    #return xys_+np.array([xcen_pix, ycen_pix])
def Iterate_LS_Strain(xys0_pix, xys1_pix, num = 10):
    solve = Get_LS_DGradient(xys0_pix, xys1_pix)
    for j in range(0,num):
        xy_pix = Simulate_Shifted_Spots(xys0_pix, linalg.expm(solve), det2lab_mat)
        dsolve = Get_LS_DGradient(xy_pix, xys1_pix)
        solve = solve+dsolve
    return linalg.expm(solve)

#def Get_LS_DGradient(xys0_pix, xys1_pix, eta=0.5, Iter = True):
solve = np.zeros((285,3,3))
for i in range(609, 894):
    xy0, xy1 = np.loadtxt('/home/fengguo/microstress/%iref'%i, skiprows=1)[:,2:4], np.loadtxt('/home/fengguo/microstress/%i200'%i, skiprows=1)[:,2:4]
    xys0, xys1 = np.zeros(xy0.shape), np.zeros(xy1.shape)
    xys0[:,0], xys1[:,0] = -xy0[:,1], -xy1[:,1]
    xys0[:,1], xys1[:,1] = xy0[:,0], xy1[:,0]
    solve[i-609]=Get_LS_DGradient(xys0, xys1, eta=0.5, Iter = True)
"""
num = len(xys0)
eta, Iter = 0.5, True
xys0_pix, xys1_pix = xys0*pix_size, xys1*pix_size
dxy = (xys1_pix-xys0_pix).reshape(2*num)
    
coef = np.zeros((2*num, 8))
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
ls_coef, ls_dxy = np.dot(coef.T, coef), np.dot(coef.T, dxy)
solve = linalg.solve(ls_coef, ls_dxy)
print solve
strain = np.zeros((3,3))
if Iter == False:
    for i in range(0,8):       
        strain = strain+solve[i]*base_strain[i]
else:
    solve = leastsq(residuals, solve, args=(k_f_u0,k_f_u1))[0]
    for i in range(0,8):     
        strain = strain+solve[i]*base_strain[i]
"""
#t=Get_LS_DGradient(xys0_pix, xys1_pix, eta=0.5, Iter = True)

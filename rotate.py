import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, sin, cos
from numpy.linalg import eig
ang = 40

F = np.load('./Ge1_100N/det50/solve_lab.npy')
F_s, strain_s = np.zeros(F.shape), np.zeros((F.shape[0],6))

cos_, sin_ = cos(ang*pi/180), sin(ang*pi/180)
mat = np.array([[1,0,0],[0,cos_,-sin_],[0,sin_,cos_]])

c11, c12, c44 = 1.66, 0.64, 0.79
ela = np.zeros((6,6))

for i in range(3):
    ela[i,i] = c11
    ela[i+3,i+3] = c44
ela[1,2],ela[2,1],ela[0,2],ela[2,0],ela[0,1],ela[1,0] = c12, c12, c12, c12, c12, c12

for i in range(F.shape[0]):
    temp = np.tensordot(F[i],mat,(0,0))
    F_s[i] = np.tensordot(temp,mat,(0,0))
    dil = (c11+2*c12)/(c12*(F_s[i,0,0]+F_s[i,2,2])+c11*F_s[i,1,1])
    print dil
    strain_s[i,0:3] = (dil*F_s[i]-1).diagonal()
    strain_s[i,3:6] = F_s[i,1,2]+F_s[i,2,1], F_s[i,0,2]+F_s[i,2,0], F_s[i,1,0]+F_s[i,0,1]
plt.plot(strain_s[:,2],'.')

F = np.load('./Ge1_140N/det50/solve_lab.npy')
F_s, strain_s = np.zeros(F.shape), np.zeros((F.shape[0],6))

cos_, sin_ = cos(ang*pi/180), sin(ang*pi/180)
mat = np.array([[1,0,0],[0,cos_,-sin_],[0,sin_,cos_]])

c11, c12, c44 = 1.66, 0.64, 0.79
ela = np.zeros((6,6))

for i in range(3):
    ela[i,i] = c11
    ela[i+3,i+3] = c44
ela[1,2],ela[2,1],ela[0,2],ela[2,0],ela[0,1],ela[1,0] = c12, c12, c12, c12, c12, c12

for i in range(F.shape[0]):
    temp = np.tensordot(F[i],mat,(0,0))
    F_s[i] = np.tensordot(temp,mat,(0,0))
    dil = (c11+2*c12)/(c12*(F_s[i,0,0]+F_s[i,2,2])+c11*F_s[i,1,1])
    print dil
    strain_s[i,0:3] = (dil*F_s[i]-1).diagonal()
    strain_s[i,3:6] = F_s[i,1,2]+F_s[i,2,1], F_s[i,0,2]+F_s[i,2,0], F_s[i,1,0]+F_s[i,0,1]
plt.plot(strain_s[:,2],'.')
plt.show()
"""
ela = np.zeros((6,6))

for i in range(3):
    ela[i,i] = 1.66
    ela[i+3,i+3] = 0.79
ela[1,2],ela[2,1],ela[0,2],ela[2,0],ela[0,1],ela[1,0] = 0.64, 0.64, 0.64, 0.64, 0.64, 0.64



solve6 = np.zeros(6)
stress6 = np.zeros((solve.shape[0],6))

for i in range(solve.shape[0]):
    for j in range(3):
        solve6[j] = solve[i,j,j]
    solve6[3] = solve[i,1,2]
    solve6[4] = solve[i,0,2]
    solve6[5] = solve[i,0,1]
    stress6[i] = np.dot(ela,solve6)
"""
    

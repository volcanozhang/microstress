import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from common import stringint

pix_size = 0.03100
DEG2RAD = np.pi/180.0

def Normalize(vector):
    length = np.sqrt(vector.__pow__(2).sum())
    return vector/length

def Detector_to_Lab(xy_pix, det2lab_mat, dd, xcen_pix, ycen_pix, pix_size):
    x_pix, y_pix = xy_pix
    temp = np.array([0.0, 0.0, dd])
    temp[0: 2] = np.array([x_pix-xcen_pix, y_pix-ycen_pix])*pix_size
    k_f = np.dot(det2lab_mat, temp)
    k_f_u = Normalize(k_f)
    return k_f_u

def Det_Lab_Mat(beta, gamma):
    beta, gamma = beta*DEG2RAD, gamma * DEG2RAD
    sin_beta, cos_beta, sin_gamma, cos_gamma = np.sin(beta), np.cos(beta), np.sin(gamma), np.cos(gamma)
    beta_mat = np.array([[1.0, 0.0, 0.0], [0.0, cos_beta, sin_beta], [0.0, -sin_beta, cos_beta]])
    gamma_mat = np.array([[cos_gamma, sin_gamma, 0.0], [-sin_gamma, cos_gamma, 0.0], [0.0, 0.0, 1.0]])
    return np.dot(beta_mat, gamma_mat)
"""
def Trace(seq0, seq1, tolerance = 4000):
    num0, num1 = seq0.shape[0], seq1.shape[0]
    paths = {i: [i] for i in range(num0)}
    inv_path = dict()
    distancesq = lambda i, j: (seq0[i, :]-seq1[j, :]).__pow__(2).sum()
    for i in range(num0):
        jdp = 0
        mindissq = distancesq(i,0)
        for j in range(1, num1):
            dissq = distancesq(i,j)
            if dissq < mindissq:
                jdp = j
                mindissq = dissq
        #print mindissq, jdp
        if mindissq > tolerance:
            del paths[i]
            continue
        if inv_path.has_key(jdp):
            if paths.has_key(inv_path[jdp]):
                del paths[inv_path[jdp]]
            del paths[i]
        else:
            inv_path[jdp] = i
            paths[i].append(jdp)
            
    return paths
"""
def Trace(k_ref, k_cur, tolerance = 0.9977):
    num0, num1 = k_ref.shape[0], k_cur.shape[0]
    paths = {i: [i] for i in range(num0)}
    inv_path = dict()
    cosine = lambda i, j: np.dot(k_ref[i], k_cur[j])
    for i in range(num0):
        jdp = 0
        max_cosine = cosine(i,0)
        for j in range(1, num1):
            cosij = cosine(i,j)
            if cosij > max_cosine:
                jdp = j
                max_cosine = cosij
        #print max_cosine, jdp
        if max_cosine < tolerance:
            del paths[i]
            continue
        if inv_path.has_key(jdp):
            if paths.has_key(inv_path[jdp]):
                del paths[inv_path[jdp]]
            del paths[i]
        else:
            inv_path[jdp] = i
            paths[i].append(jdp)
            
    return paths

#ref = np.load('pos0738_.npy')
"""
for i in range(0, 1):
    ref = np.load('pos%s_5.npy'%stringint(i,4))
    cur = np.load('pos%s_200.npy'%stringint(i+8,4))
    mch = np.array(Trace(ref, cur).values())
    if len(mch)<5:
        print i
    f_cur, f_ref = open('%i200'%i,'w'),open('%iref'%i,'w')
    f_cur.write('%i %i %i'%(len(mch),0,len(mch)))
    f_ref.write('%i %i %i'%(len(mch),0,len(mch)))
    for j in range(len(mch)):
        f_cur.write('\n%i %i %f %f'%(j+1,0,cur[mch[j,1],0],-cur[mch[j,1],1]))
        f_ref.write('\n%i %i %f %f'%(j+1,0,ref[mch[j,0],0],-ref[mch[j,0],1]))
    f_cur.close()
    f_ref.close()
"""

dd5, xcen5, ycen5, beta5, gamma5 = 59.836, 1365.74, 943.05, 0.373, 0.504
dd200, xcen200, ycen200, beta200, gamma200 = 59.836, 1365.74, 943.05, 0.373, 0.504

for n in range(0,1):
    ref = np.load('pos%s_5.npy'%stringint(n,4))
    cur = np.load('pos%s_50.npy'%stringint(n+7,4))
    k_ref, k_cur = np.zeros((ref.shape[0], 3)), np.zeros((cur.shape[0], 3))
    det2lab_mat5 = Det_Lab_Mat(beta5, gamma5)
    for i in range(ref.shape[0]):
        k_ref[i] = Detector_to_Lab(ref[i], det2lab_mat5, dd = dd5, xcen_pix = xcen5, ycen_pix = ycen5, pix_size = pix_size)
    det2lab_mat200 = Det_Lab_Mat(beta200, gamma200)
    for i in range(cur.shape[0]):
        k_cur[i] = Detector_to_Lab(cur[i], det2lab_mat200, dd = dd200, xcen_pix = xcen200, ycen_pix = ycen200, pix_size = pix_size)
    mch = np.array(Trace(k_ref, k_cur).values())

    f_cur, f_ref = open('%i200'%n,'w'),open('%iref'%n,'w')
    print len(mch)
    f_cur.write('%i %i %i'%(len(mch),0,len(mch)))
    f_ref.write('%i %i %i'%(len(mch),0,len(mch)))
    for j in range(len(mch)):
        f_cur.write('\n%i %i %f %f'%(j+1,0,cur[mch[j,1],0],-cur[mch[j,1],1]))
        f_ref.write('\n%i %i %f %f'%(j+1,0,ref[mch[j,0],0],-ref[mch[j,0],1]))
    f_cur.close()
    f_ref.close()

plt.axis([0,2594,0,2774])
#for i in range(len(mch)):
#    plt.text(ref[mch[i,0],0],ref[mch[i,0],1],'%i'%i,color ='r')
#    plt.text(cur[mch[i,1],0],cur[mch[i,1],1],'%i'%i,color='b')
plt.plot(ref[:,0],ref[:,1],'r.')
plt.plot(cur[:,0],cur[:,1],'b.')
plt.show()

"""
xy0, xy1 = np.load('pos0141_.npy'), np.load('pos0005_.npy')
path = np.array(Trace(xy0, xy1).values())

offset = 4096
bkg_offset = 110
framedim = (2594, 2774)
nb_elem = framedim[0]*framedim[1]
formatdata = np.uint16

path0 = '/home/fengguo/Data/Si1g_5N/scan/S1gscan_0005_mar.tiff'
path1 = '/home/fengguo/Data/Si1g_5N/scan/S1gscan_0141_mar.tiff'
f = open(path0, 'rb')
f.seek(offset)
image0 = np.fromfile(f, dtype = formatdata, count = nb_elem).reshape(framedim)
f.close()
f = open(path1, 'rb')
f.seek(offset)
image1 = np.fromfile(f, dtype = formatdata, count = nb_elem).reshape(framedim)
f.close()

gs = gridspec.GridSpec(1, 2)
ax0 = plt.subplot(gs[0, 0])
ax1 = plt.subplot(gs[0, 1])

ax0.imshow(np.log(image0))
ax1.imshow(np.log(image1))
for i in range(len(path)):
    ax0.text(xy0[path[i,0],0],xy0[path[i,0],1],'%i'%i)
    ax1.text(xy1[path[i,1],0],xy1[path[i,1],1],'%i'%i)
plt.show()
"""


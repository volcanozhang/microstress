import numpy as np
offset = 4096
bkg_offset = 110
framedim = (2594, 2774)
nb_elem = framedim[0]*framedim[1]
formatdata = np.uint16
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.cm as cm
#from PeakSearch import PeakSearch

path = '/home/fengguo/Data/Si1g_200N/scan/S1gscan_0100_mar.tiff'
bkg_path = '/home/fengguo/Data/back_5s_0001.tif'
fit = np.loadtxt('/home/fengguo/Data/Si1g_200N/scan/dat_S1gscan_0100_mar_LT_0.fit', skiprows = 5)
dat = np.loadtxt('/home/fengguo/Data/Si1g_200N/scan/S1gscan_0100_mar_LT_0.dat', skiprows = 1)
cor = np.loadtxt('/home/fengguo/Data/Si1g_200N/scan/dat_S1gscan_0100_mar_LT_0.cor', skiprows = 1)
f = open(path, 'rb')
f.seek(offset)
bkg_f = open(bkg_path, 'rb')
bkg_f.seek(bkg_offset)
image = np.fromfile(f, dtype = formatdata, count = nb_elem).reshape(framedim)
bkg = np.fromfile(bkg_f, dtype = formatdata, count = nb_elem).reshape(framedim)
raw_image = image+bkg-100

hkl_xy = dict()
hkl_theta = dict()
theta = cor[:,0]/2
for i in range(fit.shape[0]):
    xy = dat[fit[i, 0], 0: 2]
    hkl_xy[(int(fit[i, 2]), int(fit[i, 3]), int(fit[i, 4]))] = tuple(xy.tolist())
    hkl_theta[(int(fit[i, 2]), int(fit[i, 3]), int(fit[i, 4]))] = theta[i]
for key in hkl_xy.keys():
    #plt.text(hkl_xy[key][0],hkl_xy[key][1], '(%i%i%i)'%key)
    plt.text(hkl_xy[key][0],hkl_xy[key][1], '%f'%(np.sin(hkl_theta[key])*np.sqrt(key[0]*key[0]+key[1]*key[1]+key[2]*key[2])))
plt.imshow(raw_image)
plt.show()

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
fit = np.loadtxt('/home/fengguo/Data/Si1g_200N/scan/dat_S1gscan_0100_mar_LT_1.fit', skiprows = 5)
cor = np.loadtxt('/home/fengguo/Data/Si1g_200N/scan/S1gscan_0100_mar_LT_1.dat', skiprows = 1)
f = open(path, 'rb')
f.seek(offset)
bkg_f = open(bkg_path, 'rb')
bkg_f.seek(bkg_offset)
image = np.fromfile(f, dtype = formatdata, count = nb_elem).reshape(framedim)
bkg = np.fromfile(bkg_f, dtype = formatdata, count = nb_elem).reshape(framedim)
raw_image = image+bkg-100

hkl_xy = dict()
for i in range(fit.shape[0]):
    xy = cor[fit[i, 0], 2: 4]
    hkl_xy[(int(fit[i, 2]), int(fit[i, 3]), int(fit[i, 4]))] = tuple(xy.tolist())
for key in hkl_xy.keys():
    plt.text(hkl_xy[key][0],hkl_xy[key][1], '(%i%i%i)'%key)
plt.imshow(raw_image)
plt.show()

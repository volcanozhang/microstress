import numpy as np
import numpy.ma as ma
import os
#from scipy import optimize, stats
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

offset = 4096
bkg_offset = 110
framedim = (2594, 2774)
nb_elem = framedim[0]*framedim[1]
formatdata = np.uint16
path = os.path.join('/', 'home', 'fengguo', 'Data', 'Si1g_5N', 'scan', 'S1gscan_0034_mar.tiff')
bkg_path = os.path.join('/', 'home', 'fengguo', 'Data', 'back_5s_0001.tif')

f = open(path, 'rb')
f.seek(offset)
bkg_f = open(bkg_path, 'rb')
bkg_f.seek(bkg_offset)
image = np.fromfile(f, dtype = formatdata, count = nb_elem).reshape(framedim)
bkg = np.fromfile(bkg_f, dtype = formatdata, count = nb_elem).reshape(framedim)

if image.min() == 0:
    image = ma.masked_equal(image, 0)
    m = image.mask
    bkg = ma.array(bkg, mask = m)

raw_image = image+bkg-100
plt.imshow(np.log(raw_image))
plt.show()

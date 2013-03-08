import numpy as np
import os
#from scipy import optimize, stats
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

offset = 4096
framedim = (2594, 2774)
nb_elem = framedim[0]*framedim[1]
formatdata = np.uint16
path = '/home/fengguo/Data/21Feb13/Si1g_200N/Getdet/'
nom = 'Getdet_0034_mar.tiff'
path = os.path.join('/home', 'fengguo', 'Data', '21Feb13', 'Si1g_200N', 'Getdet', 'Getdet_0034_mar.tiff')

f = open(path, 'rb')
f.seek(offset)
raw_image = np.fromfile(f, dtype = formatdata, count = nb_elem).reshape(framedim)
plt.imshow(np.log(raw_image))
plt.show()

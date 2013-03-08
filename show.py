import numpy as np
#from scipy import optimize, stats
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

offset = 4096
framedim = (2594, 2774)
nb_elem = framedim[0]*framedim[1]
formatdata = np.uint16
path = '/home/fengguo/Data/21Feb13/Si1g_200N/Getdet/'
nom = 'Getdet_0031_mar.tiff'

f = open(path+nom, 'rb')
f.seek(offset)
raw_image = np.fromfile(f, dtype = formatdata, count = nb_elem).reshape(framedim)
plt.imshow(np.log(raw_image))

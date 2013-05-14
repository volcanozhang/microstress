import numpy as np
import finding as FIND
from common import stringint

offset = 4096
bkg_offset = 110
framedim = (2594, 2774)
nb_elem = framedim[0]*framedim[1]
formatdata = np.uint16
bkg_path = os.path.join('/', 'home', 'fengguo', 'Data', 'back_5s_0001.tif')

ref_path = '/home/fengguo/Data/Si1g_5N/scan/' + 'S1gscan_0440_mar.tiff'
path = '/home/fengguo/Data/Si1g_200N/scan/' + 'S2gscan_0001_mar.tiff'

bkg_f = open(bkg_path, 'rb')
bkg_f.seek(bkg_offset)
bkg = np.fromfile(bkg_f, dtype = formatdata, count = nb_elem).reshape(framedim)
bkg_f.close()

for i in range(0, 881):
    path = '/home/fengguo/Data/Si1g_200N/scan/' + 'S2gscan_%s_mar.tiff'%stringint(4,i)
    f = open(path, 'rb')
    f.seek(offset)
    image = np.fromfile(f, dtype = formatdata, count = nb_elem).reshape(framedim)
    raw_image = image+bkg-100
    pos=FindingPosition(raw_image, xdim = 2594, ydim = 2774)
    np.save('pos%s'%stringint(4,i),pos)
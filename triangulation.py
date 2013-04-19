from PeakSearch import PeakSearch

filename = '/home/fengguo/Data/Si1g_5N/tdet/S1gtdet_0000_mar.tiff'
res = PeakSearch(filename, IntensityThreshold = 100)[0]
f = open('XY.pts', 'w')
f.write('%i %i %i'%(res.shape[0], 1, res.shape[0]))
for i in range(res.shape[0]):
    f.write('\n%i %i %.06f %.06f'%(i+1, 0, res[i, 0], res[i,1]))
f.close()

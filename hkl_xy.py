import numpy as np

fit = np.loadtxt('dat_peak_20.fit', skiprows = 5)
cor = np.loadtxt('dat_peak_20.cor', skiprows = 1)

hkl_xy = dict()
for i in range(fit.shape[0]):
    xy = cor[fit[i, 0], 2: 4]
    hkl_xy[(int(fit[i, 2]), int(fit[i, 3]), int(fit[i, 4]))] = xy.tolist()

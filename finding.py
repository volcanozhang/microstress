import numpy as np
import numpy.ma as ma

import fitting as FIT
#import matplotlib.pyplot as plt
displacements = np.array([[1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1], [1, -1]])


def Pos_in_Circle(Radius):
    pos = [[0,0]]
    for i in range(1, Radius+1):
        pos.append((i, 0))
    for i in range(-Radius, 0):
        pos.append((i, 0))
    for j in range(1, Radius+1):
        imax = np.int(np.sqrt(Radius*Radius-j*j))
        pos.append((0, j))
        for i in range(1, imax+1):
            pos.append((i, j))
        for i in range(-imax, 0):
            pos.append((i, j))          
        pos.append((0, -j))
        for i in range(1, imax+1):
            pos.append((i, -j))
        for i in range(-imax, 0):
            pos.append((i, -j))
    return np.array(pos)

#int_image = ma.masked_equal(raw_image, 0)
def FindingPosition(raw_image, PixelNearRadius = 10, NumOfPoints = 80, xdim = 2048, ydim = 2048, Boxsize = 31, Format = np.float):
    image = ma.masked_equal(raw_image, 0).astype(Format)#/ raw_image.max()
    mask = image.mask.copy()
    
    Circle = Pos_in_Circle(PixelNearRadius)
    IniCenters = []
    #print "begin finding"
    for ipoint in range(0, NumOfPoints):
        pos = np.array([image.argmax()/ydim, image.argmax()%ydim])
        #print image.max()
        while True:
            Coverage = pos + Circle
            while Coverage.max() >= min(xdim,ydim) or Coverage.min() <= 0:
                Coverage = pos + Pos_in_Circle(PixelNearRadius-1)
            sum_cen = image[Coverage[:,0],Coverage[:,1]].sum()
            sum_max, imax = 0, 0
            for i in range(0,8):
                tmp_Coverage = Coverage + displacements[i]
                if tmp_Coverage.max() >= min(xdim,ydim) or tmp_Coverage.min() <= 0:
                    continue
                sum_tmp = image[tmp_Coverage[:,0],tmp_Coverage[:,1]].sum()
                if sum_tmp>sum_max:
                    sum_max = sum_tmp
                    imax = i
            if sum_cen>=sum_max:
                break
            else:
                pos = pos + displacements[imax]
# Check whether the spot is on the edge of the screen. If so, delete it
        trial_pos = np.zeros((2), Format)
        if pos[0]>=xdim/2:
            trial_pos[0] = pos[0]+Boxsize/2
        else:
            trial_pos[0] = pos[0]-Boxsize/2
        if pos[1]>=ydim/2:
            trial_pos[1] = pos[1]+Boxsize/2
        else:
            trial_pos[1] = pos[1]-Boxsize/2
        if pow(trial_pos[0]-xdim/2,2)+pow(trial_pos[1]-ydim/2,2)<pow(xdim/2, 2):
            IniCenters.append(pos.tolist())
            
        for i in range(0, len(Coverage)):
            image[Coverage[i, 0], Coverage[i, 1]] = ma.masked
            
    IniCenters = np.array(IniCenters)

    image.mask = mask
    #print "begin fitting"

    FittedPos = []
#PreFittedPos=[]

    for i in range(0, len(IniCenters)):
        IniX = IniCenters[i,0]
        IniY = IniCenters[i,1]
        loc_image = image[IniX - Boxsize/2: IniX + Boxsize/2 + 1, IniY - Boxsize/2: IniY + Boxsize/2 + 1].copy()
    #PreFittedPos.append([loc_image.argmax()/Boxsize,loc_image.argmax()%Boxsize])
        result = FIT.gaussfit(loc_image)
        if result != None:
            FittedPos.append((result[2: 4]-np.array([Boxsize/2,Boxsize/2])+IniCenters[i]).tolist())

    FittedPos = np.array(FittedPos)
    
    return FittedPos


#loc_image = image[]


#plt.imshow(np.log(image))
#plt.plot(IniCenters[:,1],IniCenters[:,0],'o')
#plt.show()


#remember that the datatype must be float32 (and range from 0.0 to 1.0) or uint8.
#view only pixel whose grey level is above threshold
#pro_image = np.zeros(framedim, dtype = formatfloat)
#for iraw, ipro in zip(np.nditer(raw_image), np.nditer(pro_image, op_flags=['writeonly'])):
#    ipro[...] = formatfloat(iraw)/65536.0
#plt.imshow(pro_image)
#plt.show()

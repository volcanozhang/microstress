import numpy as np
#import numpy.ma as ma
import scipy as sp
from scipy import optimize, stats

"""
dirpath = 'E:/data/Si_load0/'
nom = 'Si0023.tif'
#outnom = 'ref.tiff'
xdim = 2048
ydim = 2048
offset = 4096
formatdata = np.uint16
formatfloat = np.float

IntensityThreshold = 500
xtol = 0.001
FitPixelDev = 2.0

framedim = (xdim, ydim)
nb_elem = xdim * ydim
path = dirpath + nom
f = open(path, 'rb')
f.seek(offset)
raw_image = np.fromfile(f, dtype = formatdata, count = nb_elem).reshape(framedim)
"""
def twodgaussian(inpars):
    """Returns a 2d gaussian function of the form:
        x' = cos(rota) * x - sin(rota) * y
        y' = sin(rota) * x + cos(rota) * y
        (rota should be in degrees)
        g = b + a exp ( - ( ((x-center_x)/width_x)**2 +
        ((y-center_y)/width_y)**2 ) / 2 )

        where x and y are the input parameters of the returned function,
        and all other parameters are specified by this function

        However, the above values are passed by list.  The list should be:
        inpars = (height,amplitude,center_x,center_y,width_x,width_y,rota)

        You can choose to ignore / neglect some of the above input parameters using the following options:
            circle=0 - default is an elliptical gaussian (different x, y widths), but can reduce
                the input by one parameter if it's a circular gaussian
            rotate=1 - default allows rotation of the gaussian ellipse.  Can remove last parameter
                by setting rotate=0
            vheight=1 - default allows a variable height-above-zero, i.e. an additive constant
                for the Gaussian function.  Can remove first parameter by setting this to 0
        """
    inpars_old = inpars
    inpars = list(inpars)
    
    height = inpars.pop(0)
    height = float(height)

    amplitude, center_x, center_y = inpars.pop(0), inpars.pop(0), inpars.pop(0)
    amplitude = float(amplitude)
    center_x = float(center_x)
    center_y = float(center_y)

    width_x, width_y = inpars.pop(0), inpars.pop(0)
    width_x = float(width_x)
    width_y = float(width_y)

    rota = inpars.pop(0)
    
    rcen_x = center_x * np.cos(rota) - center_y * np.sin(rota)
    rcen_y = center_x * np.sin(rota) + center_y * np.cos(rota)

    if len(inpars) > 0:
        raise ValueError("There are still input parameters:" + str(inpars) + \
                " and you've input: " + str(inpars_old) + " circle=%d, rotate=%d, vheight=%d" % (circle, rotate, vheight))
            
    def rotgauss(x, y):
        xp = x * np.cos(rota) - y * np.sin(rota)
        yp = x * np.sin(rota) + y * np.cos(rota)
        g = height + amplitude * np.exp(
 -(((rcen_x - xp) / width_x) ** 2 + 
            ((rcen_y - yp) / width_y) ** 2) / 2.)
        return g
    
    return rotgauss

def IniParams(data):
    """
    Prepare the initial parameters for iterations
    """
    start_baseline = np.min(data)
    start_i, start_j = np.argmax(data) / data.shape[1], np.argmax(data) % data.shape[1]
    start_amplitude_1 = np.max(data) - start_baseline
    start_sigma1, start_sigma2 = 2, 2
    start_anglerot = 0
    startingparams = [start_baseline,
                    start_amplitude_1,
                    start_i,
                    start_j,
                    start_sigma1, start_sigma2, start_anglerot]
    return startingparams

def gaussfit(data, xtol=0.0000001):
    """
    Gaussian fitter with the ability to fit a variety of different forms of 2-dimensional gaussian.
    
    Input Parameters:
        data - 2-dimensional data array
        err=None - error array with same size as data array
        params=[] - initial input parameters for Gaussian function.
            (height, amplitude, x, y, width_x, width_y, rota)

    Output:
        Default output is a set of Gaussian parameters with the same shape as the input parameters
        Can also output the covariance matrix, 'infodict' that contains a lot more detail about
            the fit (see scipy.optimize.leastsq), and a message from leastsq telling what the exit
            status of the fitting routine was

        Warning: Does NOT necessarily output a rotation angle between 0 and 360 degrees.
    """
    params = IniParams(data)
    errorfunction = lambda p: np.ravel((twodgaussian(p)(*np.indices(data.shape)) - data))

    p, cov, infodict, errmsg, ier = optimize.leastsq(errorfunction, params, full_output=1, xtol=xtol)

    if ier in range(1,5):
        return p
    else:
        return None
"""
image = ma.masked_equal(raw_image, 0).astype(formatfloat)
dat=image[478-30:478+31,1149-30:1149+31]
"""

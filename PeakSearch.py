#built-in modules
import os
import time as ttt

#third party modules
import numpy as np
from scipy import optimize
import scipy.spatial.distance as ssd

# lauetools modules

def PeakSearch(filename,
                framedim=(2594, 2774),
                offset=4096,
                formatdata="uint16",
                fliprot="no",
                center=None,
                boxsizeROI=(200, 200), #use only if center != None
                PixelNearRadius=5,
                removeedge=2,
                IntensityThreshold=400,
                thresholdConvolve=200,
                paramsHat=(4, 5, 2),
                boxsize=15,
                verbose=0,
                position_definition=0,
                local_maxima_search_method=1,
                fit_peaks_gaussian=1,
                xtol=0.0000001,
                return_histo=1,
                FitPixelDev=25, #  to_reject3 parameter
                write_execution_time=1,
                Saturation_value=65535,
                Saturation_value_flatpeak=65535,
                oldversion=False
                ):
    """

    Find local intensity maxima as starting position for fitting. Return peaklist.

    From a image MARCCD:

    center = position of the ROI center in CCD frame
    boxsizeROI = dimensions of the ROI to crop the data array

    boxsize                    :  half length of the selected array:
                                    for fitting a peak
                                    for estimating the background around a peak
                                    for shifting array in second method (shifted arrays)

    IntensityThreshold        :     select as germ for potential peak pixels with intensity higher than this threshold
                                start with high value
                                If too high, few peaks are found (only the most important)
                                If too low, too much germs lead to time consuming

    PixelNearRadius:         :    pixel distance between two regions considered as peaks
                                start rather with a large value
                                if too low, there are very much peaks duplicates and
                                this is very time consuming

    local_maxima_search_method     Select method for find the local maxima, each of them will fitted
                            : 0   extract all pixel above intensity threshold
                            : 1   find pixels are highest than their neighbours in horizontal, vertical
                                    and diagonal direction (up to a given pixel distance)
                            : 2   find local hot pixels which after numerical convolution give high intensity
                                above threshold (thresholdConvolve)

    fit_peaks_gaussian      :     0   no fit
                            :    1  gaussian fit
                            :    2  lorentzian fit

    position_definition: due to various conventional habits when reading array, add some offset to fit XMAS or fit2d peak search values
                         = 0    no offset (python numpy convention)
                         = 1   XMAS offset
                         = 2   fit2d offset

    xtol  : relative error on solution (x vector)  see args for leastsq in scipy.optimize

    return_histo        : 0   3 output elements
                        : 1   4 elemts, last one is histogram of data
                        : 2   4 elemts, last one is the nb of raw blob found after convolution and threshold

    returns peak list sorted by decreasing (integrated intensity - fitted bkg)
    peak_X,peak_Y,peak_I,peak_fwaxmaj,peak_fwaxmin,peak_inclination,Xdev,Ydev,peak_bkg


    """

    if return_histo in (0, 1):
        return_nb_raw_blobs = 0
    if return_histo in (2,):
        return_nb_raw_blobs = 1
    if write_execution_time:
        t0 = ttt.time()

    # numpy read of data
    d = readoneimage(filename,
                     framedim=framedim, offset=offset, formatdata=formatdata)

    Data = np.reshape(d, framedim)

    if fliprot == "spe":
        Data = np.rot90(Data, k=1)

    elif fliprot == "vhr":
        Data = np.rot90(Data, k=3)
        Data = np.fliplr(Data)

    elif fliprot == "VHR_Feb13":
#        Data = np.rot90(Data, k=3)
        # TODO: fliplr useful ?
        Data = np.fliplr(Data)

    elif fliprot == "vhrdiamond":
        Data = np.rot90(Data, k=3)
        Data = np.fliplr(Data)

    elif fliprot == "frelon2":
        Data = np.flipud(Data)

    # peak search in a particular region of image
    if center != None:

#        imin, imax, jmin, jmax = getindices2cropArray(center, boxsizeROI, framedim)
#
#        Data = Data[imin: imax, jmin: jmax]

        framedim = (framedim[1], framedim[0])

        imin, imax, jmin, jmax = getindices2cropArray(center, boxsizeROI, framedim)

        Data = Data[jmin: jmax, imin: imax]

    if write_execution_time:
        dtread = ttt.time() - t0
        ttread = ttt.time()
        print "Read Image. Execution time : %.3f seconds" % dtread

    if return_histo:
        # from histogram, deduces
        min_intensity = max(np.amin(d), 1)  # in case of 16 integer
        max_intensity = min(np.amax(d), Saturation_value)
        print "min_intensity", min_intensity
        print "max_intensity", max_intensity
        histo = np.histogram(d,
                             bins=np.logspace(np.log10(min_intensity),
                                           np.log10(max_intensity), num=30))

    #maxhisto=np.argmax(histo[0])
    ## cumulative sum starting from hottest pixels
    #cs=cumsum(histo[0][::-1])
    #number_pixels=5000
    #pos=searchsorted(cs,number_pixels,side='right')
    #intensitythreshold=histo[1][-(pos+1)]
    #print "number_pixels",number_pixels
    #print "intensitythreshold",intensitythreshold

    #--- PRE SELECTION OF HOT PIXELS as STARTING POINTS FOR FITTING ---------

    if local_maxima_search_method == 0:
        #first method ---------- "Basic Intensity Threshold"
        print "Raw Search of local intensity maxima (pixels intensity thresholding)"
        peaklist = LocalMaxima_fromarray(Data, IntensityThreshold=IntensityThreshold)

        print "len(peaklist)", len(peaklist)

        Ipixmax = np.ones(len(peaklist)) * IntensityThreshold

    if local_maxima_search_method == 1:

        #second method ----------- "Local Maxima in a box"
        print "Using shift arrays to detect local maxima"
        peaklist, Ipixmax = LocalMaxima_ShiftArrays(Data,
                            framedim=framedim,
                            IntensityThreshold=IntensityThreshold,
                            Saturation_value=Saturation_value_flatpeak,
                            boxsize_for_probing_minimal_value_background=boxsize, # 30
                            pixeldistance_remove_duplicates=PixelNearRadius, # 25
                            nb_of_shift=boxsize)  #25

    if local_maxima_search_method == 2:
        # third method: ------------ "Convolution by a gaussian kernel"
        print "Using mexican hat convolution to detect local maxima"

        peakValConvolve, boxsizeConvolve, central_radiusConvolve = paramsHat

        Candidates = LocalMaxima_GoodCandidates(Data,
                            framedim=framedim,
                            peakValConvolve=peakValConvolve,
                            boxsizeConvolve=boxsizeConvolve,
                            central_radiusConvolve=central_radiusConvolve,
                            thresholdConvolve=thresholdConvolve, # 600 for CdTe
                            connectivity=1,
                            IntensityThreshold=IntensityThreshold,
                            boxsize_for_probing_minimal_value_background=PixelNearRadius,
                            return_nb_raw_blobs=return_nb_raw_blobs)

        if return_nb_raw_blobs == 1:
            peaklist, Ipixmax, nbrawblobs = Candidates
        else:
            peaklist, Ipixmax = Candidates
        # -------------------------------------------------------------


    if peaklist in (None, [], np.array([])) or len(peaklist) == 0:
        print "No local maxima found, change peak search parameters !!!"
        return False
    # pixel origin correction due to ROI croping
    if center != None:
        x1, y1 = center # TODO : to ne checked !!
        peaklist = peaklist + np.array([x1, y1])

    if write_execution_time:
        dtsearch = ttt.time() - float(ttread)
        ttsearch = ttt.time()

        print "Local maxima search. Execution time : %2.3f seconds" % dtsearch

    # removing some duplicates ------------
    if len(peaklist) >= 2:
        print "%d peaks in peaklist before purge" % len(peaklist), peaklist
        Xpeaklist, Ypeaklist, tokeep = removeClosePoints(peaklist[:, 0], peaklist[:, 1],
                                                            dist_tolerance=1)

        peaklist = np.array([Xpeaklist, Ypeaklist]).T
        Ipixmax = np.take(Ipixmax, tokeep)

        print "%d peaks in peaklist after purge before fitting" % len(peaklist)
    # -----------------------------------------------

    #---- ---------------------------- no FITTING ----------------------------

    # NO FIT  and return raw list of local maxima
    if fit_peaks_gaussian == 0:

        if position_definition == 1:  # XMAS offset
            peaklist[:, :2] = peaklist[:, :2] + np.array([1, 1])
        if position_definition == 2:  # fit2D offset
            peaklist[:, 0] = peaklist[:, 0] + 0.5
            peaklist[:, 1] = framedim[0] - peaklist[:, 1] + 0.5

        if verbose:
            print "%d local maxima found" % len(peaklist)
            print "20 first peaks", peaklist[:20]

        # tabpeak mimics the array built after fitting procedures
        tabpeak = np.zeros((len(peaklist[:, 0]), 10))
        tabpeak[:, 0] = peaklist[:, 0]
        tabpeak[:, 1] = peaklist[:, 1]
        tabpeak[:, 2] = Ipixmax
        #return tabpeak, peaklist, peaklist, peaklist  # no fit return raw list of local maxima
        lastelem = peaklist
        if return_nb_raw_blobs == 1:
            lastelem = nbrawblobs

        return tabpeak, peaklist, peaklist, lastelem

    #----  ----------------------------FITTING ----------------------------

    # gaussian fit
    if fit_peaks_gaussian == 1:
        type_of_function = 'gaussian'

    # lorentzian fit
    elif fit_peaks_gaussian == 2:
        type_of_function = 'lorentzian'

    else:
        raise ValueError, "optional fit_peaks_gaussian value is not understood! Must be 0,1 or 2"


    print "\n*****************\n\n"
    print "%d local maxima found" % len(peaklist)
    print "\nGaussian Fitting of local maxima\n"

#    print "framedim", framedim
#    print "offset", offset
#    print "formatdata", formatdata
#    print "fliprot", fliprot

    if center != None:
        position_start = 'centers'
    else:
        position_start = 'max'

    params, cov, info, message, baseline = readoneimage_multiROIfit(filename,
                                        peaklist,
                                        boxsize,
                                        baseline='auto', # min in ROI box
                                        startangles=0.,
                                        start_sigma1=1.,
                                        start_sigma2=1.,
                                        position_start=position_start, # 'centers' or 'max'
                                        showfitresults=0,
                                        offsetposition=0, #offset are applied after fit
                                        fitfunc=type_of_function,
                                        xtol=xtol,
                                        framedim=framedim,
                                        offset=offset,
                                        formatdata=formatdata,
                                        fliprot=fliprot)

    par = np.array(params)

    if par == []:
        print 'no fitted peaks'
        return

    peak_bkg = par[:, 0]
    peak_I = par[:, 1]
    peak_X = par[:, 2]
    peak_Y = par[:, 3]
    peak_fwaxmaj = par[:, 4]
    peak_fwaxmin = par[:, 5]
    peak_inclination = par[:, 6] % 360

    # pixel deviations from guessed initial position before fitting
    Xdev = peak_X - peaklist[:, 0]
    Ydev = peak_Y - peaklist[:, 1]
#    print 'Xdev', Xdev
#    print "Ydev", Ydev

    if write_execution_time:
        dtfit = ttt.time() - ttsearch
        tfit = ttt.time()
        print "Fit all maxima. Execution time : %.3f seconds" % dtfit

    #--- --- PEAKS REJECTION -------------------------------
    #print "peaklist[:20]",peaklist[:]
    # number of iteration screening
    to_reject = []
    k = 0
    for inf in info:
        if inf['nfev'] > 1550:
            print "k= %d   too much iteration" % k
            to_reject.append(k)
        k += 1

    # negative intensity screening
    to_reject2 = np.where((peak_bkg - baseline) < 0)[0]

    # too far found peak screening
    to_reject3 = np.where(np.sqrt(Xdev ** 2 + Ydev ** 2) > FitPixelDev)[0]

    if verbose:
        print "to_reject ...(len)", len(to_reject)
        print np.take(peaklist, to_reject, axis=0)
        print "to_reject2 ...(len)", len(to_reject2)
        print np.take(peaklist, to_reject2, axis=0)
        print np.take(peaklist, to_reject3, axis=0)

    print "After fitting, %d peaks have been rejected (final - initial position)> FitPixelDev" % \
                                                            len(to_reject3)

    ToR = set(to_reject) | set(to_reject2) | set(to_reject3)  # to reject

    ToTake = set(np.arange(len(peaklist))) - ToR

    print "ToTake", ToTake
    if len(ToTake) < 1:
        # TODO improve output branching
        if return_histo == 1:  # return 4 elements, last is the image histogram
            return None, par, peaklist, histo
        elif return_histo == 2:  # return 4 elements, last is the nb of raw blobs found after convolution and threshold
            return None, par, peaklist, nbrawblobs
        else:  # return 3 elements
            return None, par, peaklist

#    for elem in [peak_X, peak_Y,
#                               peak_I, peak_fwaxmaj, peak_fwaxmin,
#                               peak_inclination, Xdev, Ydev, peak_bkg, Ipixmax]:
#        print len(elem)

    # all peaks list building
    tabpeak = np.array([peak_X, peak_Y,
                               peak_I, peak_fwaxmaj, peak_fwaxmin,
                               peak_inclination, Xdev, Ydev, peak_bkg, Ipixmax]).T

    tabpeak = np.take(tabpeak, list(ToTake), axis=0)

    intense_rank = np.argsort(tabpeak[:, 2])[::-1]  # sort by decreasing intensity-bkg
#    print "intense_rank", intense_rank
    tabIsorted = tabpeak[intense_rank]

    if position_definition == 1:  # XMAS offset
        tabIsorted[:, :2] = tabIsorted[:, :2] + np.array([1, 1])

    elif position_definition == 2:  # fit2D offset
        tabIsorted[:, 0] = tabIsorted[:, 0] + 0.5
        tabIsorted[:, 1] = framedim[0] - tabIsorted[:, 1] + 0.5

    if verbose:
        print "\n\nIntensity sorted\n\n"
        print tabIsorted[:10]
        print "X,Y", tabIsorted[:10, :2]

    print "\n%d fitted peak(s)\n" % len(tabIsorted)

    if write_execution_time:
        dtreject = ttt.time() - tfit
        print "Reject bad peaks. Execution time : %.3f seconds " % dtreject

    if write_execution_time:
        ttotal = ttt.time() - t0
        print "\nTotal Execution time : %.3f seconds \n" % ttotal

    if return_histo == 1:  # return 4 elements, last is the image histogram
        return tabIsorted, par, peaklist, histo
    elif return_histo == 2:  # return 4 elements, last is the nb of raw blobs found after convolution and threshold
        return tabIsorted, par, peaklist, nbrawblobs
    else:  # return 3 elements
        return tabIsorted, par, peaklist
	
def readoneimage_multiROIfit(filename,
                             centers,
                             boxsize,
                             baseline='auto',
                             startangles=0.,
                             start_sigma1=1.,
                             start_sigma2=1.,
                             position_start='max',
                             fitfunc='gaussian',
                             showfitresults=1,
                             offsetposition=0,
                             verbose=0,
                             xtol=0.00000001,
                             framedim=(2048, 2048),
                             offset=4096,
                             formatdata="uint16",
                             fliprot="no"):
    """

    returns results of gaussian 2D fitting of several ROI of one image:
    params_sol = bkg,  amp  (gaussian height-bkg),
                X , Y ,
                major axis standard deviation ,minor axis standard deviation,
                major axis tilt angle / Ox

    from:
    filename:           image in mccd format (2048*2048 pixels)
    centers:             array of centers of selected ROI (MUST be an iterable object)
    boxsize:             array of boxsizes [in x, in y] direction or 
                        integer to set a square ROI for all ROI
    baseline:            'auto' (ie minimum intensity in ROI) or array of floats
    startangles:         elliptic gaussian angle (major axis with respect to X direction) ,
                        one value or array of values
    start_sigma1,start_sigma2: gaussian standard deviation (major and minor axis) in pixel,
                        one value or array of values
    position_start:     starting gaussian center:
        'max' (position of maximum intensity in ROI),
        'centers' (centre of each ROI)

    offsetposition:
        0 for no offset
        1  XMAS compatible, since XMAS consider first pixel as index 1 (in array, index starts with 0)
        2  fit2d, since fit2d for peaksearch put pixel labelled n at the position n+0.5 (between n and n+1)


    # TODO: setting list of initial guesses can be improve with
    scipy.ndimages of a concatenate array of multiple slices?
    """
    # read data
    Data = readoneimage_manycrops(filename, centers, boxsize,
                                framedim=framedim, offset=offset,
                                formatdata=formatdata, fliprot=fliprot)

    # setting initial guessed values for each center
    nb_Images = len(Data)
    print "nb of images to fit ... in  readoneimage_multiROIfit()", nb_Images
    if baseline == 'auto':  # background height
        list_min = []
        for k, dd in enumerate(Data):
            print "k, dd.shape", k, dd.shape
            list_min.append(np.amin(dd))
        start_baseline = list_min
    else:  # input numerical value array
        start_baseline = baseline

    if type(startangles) == type(1.2):
        start_anglerot = startangles * np.ones(nb_Images)
    else:
        start_anglerot = startangles

    if type(start_sigma1) == type(1.2):
        start_sigma1 = start_sigma1 * np.ones(nb_Images)
    if type(start_sigma2) == type(1.2):
        start_sigma2 = start_sigma2 * np.ones(nb_Images)

    if isinstance(boxsize, (int, float)):
        xboxsize, yboxsize = int(boxsize), int(boxsize)
    else:
        xboxsize, yboxsize = boxsize

    xboxsize = xboxsize * np.ones(nb_Images)
    yboxsize = yboxsize * np.ones(nb_Images)

    start_j = []
    start_i = []
    start_amplitude = []

    if position_start in ('centers', 'center'):  # starting position  from input center
        start_j, start_i = yboxsize, xboxsize
        start_amplitude = []
        d = 0
        for dd in Data:
            start_amplitude.append(dd[int(start_j[d]), int(start_i[d])] - start_baseline[d])
            d += 1
    elif position_start == 'max':  #starting position  from maximum intensity in dat

        d = 0
        for dd in Data:
            start_j.append(np.argmax(dd) / dd.shape[1])
            start_i.append(np.argmax(dd) % dd.shape[1])
            start_amplitude.append(np.amax(dd) - start_baseline[d])
            d += 1

    ##        else:
    ##                start_j,start_i=posxy
    ##                start_amplitude=start_amplitude=dat[start_j,start_i]-start_baseline

    startingparams_zip = np.array([start_baseline,
                                   start_amplitude,
                                   start_j, start_i,
                                   start_sigma1, start_sigma2,
                                   start_anglerot])


    RES_params = []
    RES_cov = []
    RES_infodict = []
    RES_errmsg = []

    if verbose:
        print "startingparams_zip", startingparams_zip

    k_image = 0
    for startingparams in startingparams_zip.T:
        #if (k_image%25) == 0: print "%d/%d"%(k_image,nb_Images)
        if verbose:
            print "startingparams", startingparams
        if fitfunc == 'gaussian':
            params, cov, infodict, errmsg = \
                        gaussfit(Data[k_image],
                                            err=None,
                                            params=startingparams,
                                            autoderiv=1,
                                            return_all=1,
                                            circle=0,
                                            rotate=1,
                                            vheight=1,
                                            xtol=xtol)

        elif fitfunc == 'lorentzian':
            params, cov, infodict, errmsg = \
                        lorentzfit(Data[k_image],
                                            err=None,
                                            params=startingparams,
                                            autoderiv=1,
                                            return_all=1,
                                            circle=0,
                                            rotate=1,
                                            vheight=1,
                                            xtol=xtol)

        if showfitresults:
            #print "startingparams"
            #print startingparams
            print "\n *****fitting results ************\n"
            print params
            print "background intensity:                        %.2f" % params[0]
            print "Peak amplitude above background              %.2f" % params[1]
            print "pixel position (X)                   %.2f" % (params[3] - xboxsize[k_image] + centers[k_image][0])  # WARNING Y and X are exchanged in params !
            print "pixel position (Y)                   %.2f" % (params[2] - yboxsize[k_image] + centers[k_image][1])
            print "std 1,std 2 (pix)                    ( %.2f , %.2f )" % (params[4], params[5])
            print "e=min(std1,std2)/max(std1,std2)              %.3f" % (min(params[4], params[5]) / max(params[4], params[5]))
            print "Rotation angle (deg)                 %.2f" % (params[6] % 360)
            print "************************************\n"
        bkg_sol, amp_sol, Y_sol, X_sol, std1_sol, std2_sol, ang_sol = params

        RES_cov.append(cov)
        RES_infodict.append(infodict)
        RES_errmsg.append(errmsg)

        params_sol = np.array([bkg_sol,
                            amp_sol,
                            X_sol - xboxsize[k_image] + centers[k_image][0],
                            Y_sol - yboxsize[k_image] + centers[k_image][1],
                            std1_sol,
                            std2_sol,
                            ang_sol])  # now X,Y in safest order

        if offsetposition == 1:
            # PATCH: To match peak positions given by XMAS
            # (confusion coming from counting array indices from 0 or 1...)
            # array fitted by python module see pixels at lower position
            params_sol[3] = params_sol[3] + 1.
            params_sol[2] = params_sol[2] + 1.
            # End of PATCH

        elif offsetposition == 2:  # see Compute_data2thetachi() in find2thetachi.py
            # PATCH: To match peak positions given by peaksearch of fit2D
            # in fit2D graphics window first pixel labelled 1 is for the
            # peaksearch located at position in  between 0 and 1 (ie 0.5)
            params_sol[3] = params_sol[3] + 0.5
            params_sol[2] = (framedim[0] - params_sol[2]) + 0.5  # TODO: tocheck dim[0] or dim[1]
            # End of PATCH

        RES_params.append(params_sol)

        k_image += 1

    return RES_params, RES_cov, RES_infodict, RES_errmsg, start_baseline
def readoneimage(filename, framedim=(2048, 2048), dirname=None,
                 offset=4096, formatdata="uint16"):
    """
    returns a 1d array of integers from a binary image file (full data)
    """

    nb_elem = framedim[0] * framedim[1]

    if dirname == None:
        dirname = os.curdir

    f = open(os.path.join(dirname, filename), 'rb')
    f.seek(offset)

    #d=scipy.io.fread(f,2048*2048,np.oldnumeric.Int16)
    #d = scipy.io.fread(f,nb_elem,np.oldnumeric.UInt16)
    d = np.fromfile(f, dtype=formatdata, count=nb_elem)
    f.close()
    return d
def LocalMaxima_ndimage(Data,
                        peakVal=4,
                        boxsize=5,
                        central_radius=2,
                        threshold=1000,
                        connectivity=1,
                        returnfloatmeanpos=0):

    """
    returns (float) i,j positions in array of each blob
    (peak, spot, assembly of hot pixels or whatever)

    input:

    peakVal, boxsize, central_radius    :
        parameters for numerical convolution with a mexican-hat-like kernel

    threshold :
        intensity threshold of filtered Data (by convolution with the kernel)
        above which blob signal will be considered
        if = 0 : take all blobs at the expense of processing time

    connectivity : 
        1 for filled square 3*3 connectivity
        0 for 3*3 star like connectivity
    
    output:
    array (n,2): array of 2 indices
    """
    aa = LocalMaxima_Convolve(Data,
                        peakVal=peakVal,
                        boxsize=boxsize,
                        central_radius=central_radius)

    print "Histogram after convolution with Mexican Hat"
    print np.histogram(aa)

    thraa = np.where(aa > threshold, 1, 0)

    if connectivity == 0:
        star = np.eye(3)
        ll, nf = ndimage.label(thraa, structure=star)
    elif connectivity == 1:
        ll, nf = ndimage.label(thraa, structure=np.ones((3, 3)))
    elif connectivity == 2:
        star = np.eye(3)
        ll, nf = ndimage.label(thraa, structure=np.array([[1, 1, 1], [0, 1, 0], [1, 1, 1]]))

    meanpos = np.array(ndimage.measurements.center_of_mass(thraa,
                                                                ll,
                                                                np.arange(1, nf + 1)),
                                                                dtype=np.float)
    if returnfloatmeanpos:
        return meanpos
    else:
        return np.array(meanpos, dtype=np.int)


def LocalMaxima_Convolve(Data,
                        peakVal=4,
                        boxsize=5,
                        central_radius=2
                        ):

    """
    Convolve Data array witn mexican-hat kernel 

    inputs:
    Data                            : 2D array containing pixel intensities
    peakVal > central_radius        : defines pixel distance from box center where weights are positive
                                    (in the middle) and negative farther to converge back to zero
    boxsize                            : size of the box

    ouput:
    array  (same shape as Data)
    """

    #from scipy import ndimage

    #outa = np.zeros(Data.shape)
    #ndimage.filters.gaussian_laplace(d,(5,5),output=outa)

    #whole_structure= createstructure(10, 10)-2*createstructure(10, 7)+4*createstructure(10, 5)
    #bb= ndimage.convolve(d,whole_structure)

    #bb=ndimage.morphology.white_tophat(Data,(boxsize,boxsize))
    #mexicanhat = array(LoGArr((10,10),r0=6,peakVal=4),dtype= int16)

    mexicanhat = LoGArr((boxsize, boxsize), r0=central_radius, peakVal=peakVal)
    mexicanhat = mexicanhat - sum(mexicanhat) / mexicanhat.size
    bb = ndimage.convolve(np.array(Data, dtype=np.float32), mexicanhat)

    return bb

def LocalMaxima_GoodCandidates(Data,
                        framedim=(2048, 2048),
                        peakValConvolve=4,
                        boxsizeConvolve=5,
                        central_radiusConvolve=2,
                        thresholdConvolve=1000,
                        connectivity=1,
                        IntensityThreshold=500,
                        boxsize_for_probing_minimal_value_background=30,
                        return_nb_raw_blobs=0):  # full side length

    """
    return local maxima position and amplitude as good candidates for further fit
    Use convolution by a mexican hat like kernel
    Determine blobs position
    Then select blobs according to their absolute intensity (with background...)

    inputs:
    Data                            : 2D array containing pixel intensities

    peakValConvolve, boxsizeConvolve, central_radiusConvolve : convolution parameter
    thresholdConvolve                    : minimum threshold (expressed in convoluted array intensity)
                                        under which convoluted blob is rejected
                                        can be zero (all blobs are accepted but time consuming)
    connectivity                    : shape of connectivity pattern to consider pixels belonging to the
                                    same blob
                                    1 = filled square  (1 pixel connected to 8 neighbours)
                                    0 = star (4 neighbours in vertical and horizontal direction)

    IntensityThreshold                : minimum array intensity threshold to accept blob
    boxsize_for_probing_minimal_value_background        : boxsize to evaluate the background and the blob amplitude 
    """

    dataimage_ROI = Data

    peak = LocalMaxima_ndimage(dataimage_ROI,
                        peakVal=peakValConvolve,
                        boxsize=boxsizeConvolve,
                        central_radius=central_radiusConvolve,
                        threshold=thresholdConvolve,
                        connectivity=connectivity)

    peak = (peak[:, 0], peak[:, 1])

    intensity_localmaxima = dataimage_ROI[peak]

    peaki = peak[0]
    peakj = peak[1]

    # building an array of hot pixels (2 coordinates)
    Yarray = peakj
    Xarray = peaki
    peaklist = np.array([Xarray, Yarray]).T

#    print "peaklist", peaklist
    print "%d local maxima have been found from convolution method" % len(peaklist)

    # probing background and maximal intensity in boxsize
    #
    ptp_boxsize = boxsize_for_probing_minimal_value_background

    # first method ---------------------------  
    tabptp = []  # tab of min and max around each peak
    tabposmax = []  ## tab of position of hottest pixel close to that found after convolution
    for k in range(len(peaklist)):
        print "k in LocalMaxima_GoodCandidates", k
        print "dataimage_ROI.shape", dataimage_ROI.shape
        print "peaklist[k]", peaklist[k]
        minimaxi, maxpos = minmax(dataimage_ROI,
                                 peaklist[k],
                                ptp_boxsize,
                                framedim=framedim,
                                withmaxpos=1)
        tabptp.append(minimaxi)  # TODO check order of framedim
        tabposmax.append(maxpos)  ## new

    ar_ptp = np.array(tabptp)
    ar_posmax = np.array(tabposmax)  ## new

#    print "ar_ptp", ar_ptp
    # -------------------------------------------

    # second method ---------------
    #tabptp = minmax_fast(dataimage_ROI,tuple(transpose(peaklist)),boxsize=(boxsize_for_probing_minimal_value_background,boxsize_for_probing_minimal_value_background))
    #ar_ptp = array(tabptp).T
    # ---------------------------------


    #ar_amp = np.subtract(ar_ptp[:,1],ar_ptp[:,0])
    ar_amp = np.subtract(intensity_localmaxima, ar_ptp[:, 0])
    amp_rank = np.argsort(ar_amp)[::-1]

    peaklist_sorted = peaklist[amp_rank]
    ptp_sorted = ar_ptp[amp_rank]
    amp_sorted = ar_amp[amp_rank]
    posmax_sorted = ar_posmax[amp_rank]  ## new
    # thresholding on peak-to-peak amplitude
    threshold_amp = IntensityThreshold

    cond = np.where(amp_sorted > threshold_amp)
    th_peaklist = peaklist_sorted[cond]
    th_ar_ptp = ptp_sorted[cond]
    th_ar_amp = amp_sorted[cond]
    th_ar_pos = posmax_sorted[cond]  ## new

    ##### peak positions that will be returned are the hottest pixels
    th_peaklist = th_ar_pos

    print "%d local maxima found after thresholding above %d amplitude above local background" % (len(th_ar_amp), threshold_amp)

    npeaks = np.shape(th_peaklist)[0]
    Ipixmax = np.zeros(npeaks, dtype=int)
    #print np.shape(Data)
    for i in range(npeaks):
        #Ipixmax[i]=Data[th_peaklist[i,0],th_peaklist[i,1]]
        Ipixmax[i] = th_ar_ptp[i][1]

    if return_nb_raw_blobs == 1:
        return np.fliplr(th_peaklist), Ipixmax, len(peaklist)
    else:
        return np.fliplr(th_peaklist), Ipixmax


def LocalMaxima_ShiftArrays(Data,
                        framedim=(2048, 2048),
                        IntensityThreshold=500,
                        Saturation_value=65535,
                        boxsize_for_probing_minimal_value_background=30, # full side length
                        nb_of_shift=25,
                        pixeldistance_remove_duplicates=25,
                        verbose=0):

    import networkx as NX
    #time_0 = ttt.time()
    #pilimage,dataimage=readoneimage_full(filename)

    xminfit2d, xmaxfit2d, yminfit2d, ymaxfit2d = 1, framedim[1], 1, framedim[0]

        # warning i corresponds to y
        # j corresponds to x
        # change nom xminf2d => xminfit2d pour coherence avec le reste

    #imin,imax,jmin,jmax=2048-ymaxfit2d,2048-yminfit2d,xminfit2d,xmaxfit2d
    imin, imax, jmin, jmax = framedim[0] - ymaxfit2d, \
                            framedim[0] - yminfit2d, \
                            xminfit2d - 1, \
                            xmaxfit2d - 1

    #dataimage_ROI=dataimage[imin:imax,jmin:jmax]# array index   i,j
    ## fit2d index:  X=j Y=2048-i

    dataimage_ROI = Data

    print "searching local maxima for non saturated consecutive pixels"

    peak = localmaxima(dataimage_ROI, nb_of_shift, diags=1)

    print "Done...!"
    print peak
    #print "execution time : %f  secondes"%(ttt.time()-time_0)

    intensity_localmaxima = dataimage_ROI[peak]
    #print intensity_localmaxima

    # SATURATION handling ------------------------------
    # if the top of the local maximum has at least two pixels with the same intensity, this maximum is not detected
    # this generally the case for saturated peaks
    # saturation value : saturation above which we will take into account the pixel
    # this value may be lower than the 2^n bits value to handle unfortunately very flat weak peak with 2 neighbouring pixels
    # Saturation_value = 65535 for mccd 
    Size_of_pixelconnection = 20
    print "Saturation value for flat top peak handling", Saturation_value
    sat_pix = np.where(dataimage_ROI >= Saturation_value)

    if verbose:
        print "positions of saturated pixels \n", sat_pix
    print "nb of saturated pixels", len(sat_pix[0])
    sat_pix_mean = None

    # there is at least one peak above or equal to the Saturation_value threshold
    # loop over saturated pixels
    if len(sat_pix[0]) > 0:

        if verbose:
            print "positions of saturated pixels \n", sat_pix

        if 1:  # use of graph algorithms
            sat_pix = np.column_stack(sat_pix)

            disttable_sat = ssd.pdist(sat_pix, 'euclidean')
            sqdistmatrix_sat = ssd.squareform(disttable_sat)
            # building adjencymat

            a, b = np.indices(sqdistmatrix_sat.shape)
            indymat = np.triu(b) + np.tril(a)
            cond2 = np.logical_and(sqdistmatrix_sat < Size_of_pixelconnection,
                                sqdistmatrix_sat > 0)
            adjencymat = np.where(cond2, indymat, 0)

            #print "before networkx"
            print "networkx version", NX.__version__
            GGraw = NX.to_networkx_graph(adjencymat, create_using=NX.Graph())
            list_of_cliques = NX.find_cliques(GGraw)
            #print "after networkx"

            # now find average pixel of each clique
            sat_pix_mean = []
            for clique in list_of_cliques:
                ii, jj = np.mean(sat_pix[clique], axis=0)
                sat_pix_mean.append([int(ii), int(jj)])

            sat_pix_mean = np.array(sat_pix_mean)
            print "Mean position of saturated pixels blobs = \n", sat_pix_mean

        if 0:  # of scipy.ndimage

            df = ndimage.gaussian_filter(dataimage_ROI, 10)

            #histo = np.histogram(df)
            #print "histogram",histo
            #print "maxinten",np.amax(df)
            threshold_for_measurements = np.amax(df) / 10. #histo[1][1]# 1000  pour CdTe # 50 pour Ge

            tG = np.where(df > threshold_for_measurements, 1, 0)
            ll, nf = ndimage.label(tG)  #, structure = np.ones((3,3)))
            meanpos = \
            np.array(ndimage.measurements.center_of_mass(tG,
                                                          ll,
                                                          np.arange(1, nf + 1)),
                                                            dtype=float)
            #meanpos = np.fliplr(meanpos)  # this done later

            #print "meanpos",meanpos

            sat_pix_mean = meanpos
    else:
        print "No pixel saturation"
    # SATURATION handling -(End) --------------------------------------------------------

    # x,y from localmaxima is a matter of convention

    peaki = peak[0] + imin
    peakj = peak[1] + jmin

    # building an array of hot pixels (2 coordinates)
    Yarray = peakj
    Xarray = peaki
    peaklist = np.array([Xarray, Yarray]).T

    #print peaklistfit2D
    #print peaklist[100:150]
    print "%d local maxima have been found" % len(peaklist)

    # probing background and maximal intensity in boxsize
    #
    ptp_boxsize = boxsize_for_probing_minimal_value_background

    tabptp = []
    for k in range(len(peaklist)):
        tabptp.append(minmax(dataimage_ROI, peaklist[k], ptp_boxsize, framedim=framedim))

    ar_ptp = np.array(tabptp)
    #ar_amp = np.subtract(ar_ptp[:,1],ar_ptp[:,0])
    ar_amp = np.subtract(intensity_localmaxima, ar_ptp[:, 0])
    amp_rank = np.argsort(ar_amp)[::-1]

    peaklist_sorted = peaklist[amp_rank]
    #ptp_sorted = ar_ptp[amp_rank]
    amp_sorted = ar_amp[amp_rank]
    # thresholding on peak-to-peak amplitude
    threshold_amp = IntensityThreshold

    cond = np.where(amp_sorted > threshold_amp)
    th_peaklist = peaklist_sorted[cond]
    th_ar_amp = amp_sorted[cond]

    print "%d local maxima found after thresholding above %d amplitude above local background" % (len(th_ar_amp), threshold_amp)

    # remove duplicates (close points), the most intense pixel is kept
    # minimum distance between hot pixel
    # it corresponds both to distance between peaks and peak size ...
    pixeldistance = pixeldistance_remove_duplicates

    disttable = ssd.pdist(th_peaklist, 'euclidean')
    maxdistance = np.amax(disttable)
    sqdistmatrix = ssd.squareform(disttable)
    distmatrix = sqdistmatrix + np.eye(sqdistmatrix.shape[0]) * maxdistance

    si, fi = np.where(distmatrix < pixeldistance)

    index_todelete = np.where(fi > si, fi, si)

    purged_pklist = np.delete(th_peaklist, index_todelete, axis=0)  # np.delete
    purged_amp = np.delete(th_ar_amp, index_todelete)

    print "%d local maxima found after removing duplicates (minimum intermaxima distance = %d)" % (len(purged_amp), pixeldistance)

    #print "execution time : %f  secondes"%( ttt.time() - time_0)

    # merging different kind of peaks
    if sat_pix_mean != None:
        print "Merging saturated and normal peaks"
        print "number of saturated peaks : ", np.shape(sat_pix_mean)[0]
        purged_pklist = np.vstack((sat_pix_mean, purged_pklist))

    if 0:  # check if there are still close hot pixels
        disttable_c = ssd.pdist(purged_pklist, 'euclidean')
        maxdistance_c = np.amax(disttable_c)
        sqdistmatrix_c = ssd.squareform(disttable_c)
        distmatrix_c = sqdistmatrix_c + np.eye(sqdistmatrix_c.shape[0]) * maxdistance_c

        print np.where(distmatrix_c < pixeldistance)  # must be (array([], dtype=int64), array([], dtype=int64))

    #print "purged_pklist", purged_pklist
    print "shape(purged_pklist)", np.shape(purged_pklist)
    npeaks = np.shape(purged_pklist)[0]
    Ipixmax = np.zeros(npeaks, dtype=int)
    #print np.shape(Data)

    for i in range(npeaks):
        Ipixmax[i] = Data[purged_pklist[i, 0], purged_pklist[i, 1]]
        #print "Ipixmax = ", Ipixmax

    return np.fliplr(purged_pklist), Ipixmax

def LocalMaxima_fromarray(Data, IntensityThreshold=400):

    """
    return center of mass of each blobs composes by pixels above 
    IntensityThreshold
    
    !warning!: in this procedure alone axes X,Y are swaped with respect to 
    XMAS, fit2D convention and lauetools display !!
    """
    thraa = np.where(Data > IntensityThreshold, 1, 0)

#    star = array([[0,1,0],[1,1,1],[0,1,0]])
    #ll, nf = ndimage.label(thraa, structure=np.ones((3,3)))
    #ll, nf = ndimage.label(thraa, structure=star)
    ll, nf = ndimage.label(thraa)

    print "nb of peaks in LocalMaxima_fromarray", nf
    #meanpos = np.zeros((nf,2))
    #for k in range(nf):
        #meanpos[k] = np.mean(np.where((ll == k),axis=1)

    #ndimage.find_objects(ll)
    meanpos = \
    np.array(ndimage.measurements.center_of_mass(thraa,
                                                       ll,
                                                       np.arange(1, nf + 1)),
                                                        dtype=float)
    if len(np.shape(meanpos)) > 1:
        meanpos = np.fliplr(meanpos)
    else:
        meanpos = np.roll(meanpos, 1)

    return meanpos

def gaussfit(data, err=None, params=[], autoderiv=1, return_all=0, circle=0, rotate=1, vheight=1, xtol=0.0000001):
    """
    Gaussian fitter with the ability to fit a variety of different forms of 2-dimensional gaussian.
    
    Input Parameters:
        data - 2-dimensional data array
        err=None - error array with same size as data array
        params=[] - initial input parameters for Gaussian function.
            (height, amplitude, x, y, width_x, width_y, rota)
            if not input, these will be determined from the moments of the system, 
            assuming no rotation
        autoderiv=1 - use the autoderiv provided in the lmder.f function (the alternative
            is to us an analytic derivative with lmdif.f: this method is less robust)
        return_all=0 - Default is to return only the Gaussian parameters.  See below for
            detail on output
        circle=0 - default is an elliptical gaussian (different x, y widths), but can reduce
            the input by one parameter if it's a circular gaussian
        rotate=1 - default allows rotation of the gaussian ellipse.  Can remove last parameter
            by setting rotate=0
        vheight=1 - default allows a variable height-above-zero, i.e. an additive constant
            for the Gaussian function.  Can remove first parameter by setting this to 0

    Output:
        Default output is a set of Gaussian parameters with the same shape as the input parameters
        Can also output the covariance matrix, 'infodict' that contains a lot more detail about
            the fit (see scipy.optimize.leastsq), and a message from leastsq telling what the exit
            status of the fitting routine was

        Warning: Does NOT necessarily output a rotation angle between 0 and 360 degrees.
    """
    if params == []:
        params = (momentsr(data, circle, rotate, vheight))
    if err == None:
        errorfunction = lambda p: np.ravel((twodgaussian(p, circle, rotate, vheight)(*np.indices(data.shape)) - data))
    else:
        errorfunction = lambda p: np.ravel((twodgaussian(p, circle, rotate, vheight)(*np.indices(data.shape)) - data) / err)
    if autoderiv == 0:
        # the analytic derivative, while not terribly difficult, is less efficient and useful.  I only bothered
        # putting it here because I was instructed to do so for a class project - please ask if you would like 
        # this feature implemented
        raise ValueError("I'm sorry, I haven't implemented this feature yet.")
    else:
        p, cov, infodict, errmsg, success = optimize.leastsq(errorfunction, params, full_output=1, xtol=xtol)
    if  return_all == 0:
        return p
    elif return_all == 1:
        return p, cov, infodict, errmsg
	
def twodgaussian(inpars, circle, rotate, vheight):
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
    if vheight == 1:
        height = inpars.pop(0)
        height = float(height)
    else:
        height = float(0)
    amplitude, center_x, center_y = inpars.pop(0), inpars.pop(0), inpars.pop(0)
    amplitude = float(amplitude)
    center_x = float(center_x)
    center_y = float(center_y)
    if circle == 1:
        width = inpars.pop(0)
        width_x = float(width)
        width_y = float(width)
    else:
        width_x, width_y = inpars.pop(0), inpars.pop(0)
        width_x = float(width_x)
        width_y = float(width_y)
    if rotate == 1:
        rota = inpars.pop(0)
        rota = np.pi / 180. * float(rota)
        rcen_x = center_x * np.cos(rota) - center_y * np.sin(rota)
        rcen_y = center_x * np.sin(rota) + center_y * np.cos(rota)
    else:
        rcen_x = center_x
        rcen_y = center_y
    if len(inpars) > 0:
        raise ValueError("There are still input parameters:" + str(inpars) + \
                " and you've input: " + str(inpars_old) + " circle=%d, rotate=%d, vheight=%d" % (circle, rotate, vheight))

    def rotgauss(x, y):
        if rotate == 1:
            xp = x * np.cos(rota) - y * np.sin(rota)
            yp = x * np.sin(rota) + y * np.cos(rota)
        else:
            xp = x
            yp = y
        g = height + amplitude * np.exp(
 -(((rcen_x - xp) / width_x) ** 2 +
            ((rcen_y - yp) / width_y) ** 2) / 2.)
        return g

    return rotgauss
    
def lorentzfit(data, err=None, params=[], autoderiv=1, return_all=0, circle=0, rotate=1, vheight=1, xtol=0.0000001):
    """
    Lorentzian fitter with the ability to fit a variety of different forms of 2-dimensional gaussian.

    Input Parameters:
        data - 2-dimensional data array
        err=None - error array with same size as data array
        params=[] - initial input parameters for Gaussian function.
            (height, amplitude, x, y, width_x, width_y, rota)
            if not input, these will be determined from the moments of the system,
            assuming no rotation
        autoderiv=1 - use the autoderiv provided in the lmder.f function (the alternative
            is to us an analytic derivative with lmdif.f: this method is less robust)
        return_all=0 - Default is to return only the Gaussian parameters.  See below for
            detail on output
        circle=0 - default is an elliptical gaussian (different x, y widths), but can reduce
            the input by one parameter if it's a circular gaussian
        rotate=1 - default allows rotation of the gaussian ellipse.  Can remove last parameter
            by setting rotate=0
        vheight=1 - default allows a variable height-above-zero, i.e. an additive constant
            for the Gaussian function.  Can remove first parameter by setting this to 0

    Output:
        Default output is a set of Gaussian parameters with the same shape as the input parameters
        Can also output the covariance matrix, 'infodict' that contains a lot more detail about
            the fit (see scipy.optimize.leastsq), and a message from leastsq telling what the exit
            status of the fitting routine was

        Warning: Does NOT necessarily output a rotation angle between 0 and 360 degrees.
    """
    if params == []:
        params = (momentsr(data, circle, rotate, vheight))
    if err == None:
        errorfunction = lambda p: ravel((twodlorentzian(p, circle, rotate, vheight)(*indices(data.shape)) - data))
    else:
        errorfunction = lambda p: ravel((twodlorentzian(p, circle, rotate, vheight)(*indices(data.shape)) - data) / err)
    if autoderiv == 0:
        # the analytic derivative, while not terribly difficult, is less efficient and useful.  I only bothered
        # putting it here because I was instructed to do so for a class project - please ask if you would like
        # this feature implemented
        raise ValueError("I'm sorry, I haven't implemented this feature yet.")
    else:
        p, cov, infodict, errmsg, success = optimize.leastsq(errorfunction, params, full_output=1, xtol=xtol)
    if  return_all == 0:
        return p
    elif return_all == 1:
        return p, cov, infodict, errmsg
	
def twodlorentzian(inpars, circle, rotate, vheight):
    """Returns a 2d gaussian function of the form:
        x' = cos(rota) * x - sin(rota) * y
        y' = sin(rota) * x + cos(rota) * y
        (rota should be in degrees)
        g = b + a /( 1 + 4*(x-center_x)**2/width_x +  4*(y-center_y)**2/width_y  )

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
    if vheight == 1:
        height = inpars.pop(0)
        height = float(height)
    else:
        height = float(0)
    amplitude, center_x, center_y = inpars.pop(0), inpars.pop(0), inpars.pop(0)
    amplitude = float(amplitude)
    center_x = float(center_x)
    center_y = float(center_y)
    if circle == 1:
        width = inpars.pop(0)
        width_x = float(width)
        width_y = float(width)
    else:
        width_x, width_y = inpars.pop(0), inpars.pop(0)
        width_x = float(width_x)
        width_y = float(width_y)
    if rotate == 1:
        rota = inpars.pop(0)
        rota = pi / 180. * float(rota)
        rcen_x = center_x * cos(rota) - center_y * sin(rota)
        rcen_y = center_x * sin(rota) + center_y * cos(rota)
    else:
        rcen_x = center_x
        rcen_y = center_y
    if len(inpars) > 0:
        raise ValueError("There are still input parameters:" + str(inpars) + \
                " and you've input: " + str(inpars_old) + " circle=%d, rotate=%d, vheight=%d" % (circle, rotate, vheight))

def localmaxima(DataArray, n, diags=1, verbose=0):
    """
    from DataArray 2D  returns (array([i1,i2,...,ip]),array([j1,j2,...,jp]))
    of indices where pixels value is higher in two direction up to n pixels

    this tuple can be easily used after in the following manner:
    DataArray[tupleresult] is an array of the intensity of the hottest pixels in array

    in similar way with only four cardinal directions neighbouring (found in the web):
    import numpy as N
    def local_minima(array2d):
        return ((array2d <= np.roll(array2d,  1, 0)) &
                (array2d <= np.roll(array2d, -1, 0)) &
                (array2d <= np.roll(array2d,  1, 1)) &
                (array2d <= np.roll(array2d, -1, 1)))

    WARNING: flat top peak are not detected !!
    """
    dim = len(np.shape(DataArray))

    if diags:
        c, alll, allr, alld, allu, diag11, diag12, diag21, diag22 = \
            shiftarrays_accum(DataArray, n, dimensions=dim, diags=diags)
        flag = np.greater(c, alll[0])
        for elem in alll[1:] + allr + alld + allu + \
                                    diag11 + diag12 + diag21 + diag22:
            flag = flag * np.greater(c, elem)
    else:
        c, alll, allr, alld, allu = shiftarrays_accum(DataArray,
                                                      n,
                                                      dimensions=dim,
                                                      diags=diags)
        flag = np.greater(c, alll[0])
        for elem in alll[1:] + allr + alld + allu:
            flag = flag * np.greater(c, elem)

    peaklist = np.nonzero(flag)  # in c frame index

    if verbose:
            print "value local max", c[peaklist]
            print "value from original array ", DataArray[tuple(np.array(peaklist) + n)]
            print "positions of local maxima in original frame index", tuple(np.array(peaklist) + n)

    # first slow index array , then second fast index array
    return tuple(np.array(peaklist) + n)
    
def shiftarrays_accum(Data_array, n, dimensions=1, diags=0):
    """
    idem than shiftarrays() but with all intermediate shifted arrays
    1D
    returns 3 arrays corresponding to shifted arrays
    by n in two directions and original one
    2D
    returns 5 arrays corresponding to shifted arrays
    by n in two directions and original one

    these arrays are ready for comparison with eg np.greater

    Data_array must have shape (slowdim,fastdim) so that
    slowdim-2*n>=1 and fastdim-2*n>=1
    (ie central array with zero shift has some elements)

    TODO: replace append by a pre allocated array
    """
    if n <= 0:
        raise ValueError, "shift value must be positive"

    if dimensions == 2:
        if diags:
            shift_zero = Data_array[n:-n, n:-n]

            allleft = []
            allright = []
            allup = []
            alldown = []
            alldiagleftdown = []  # diag "y=x"
            alldiagrightup = []  # diag "y=x"
            alldiagrightdown = []  #  diah "y=-x"
            alldiagleftup = []  #  diah "y=-x"

            for k in np.arange(1, n + 1)[::-1]:

                allleft.append(Data_array[n - k:-(n + k), n:-n])
                alldown.append(Data_array[n:-n, n - k:-(n + k)])
                alldiagrightdown.append(Data_array[n - k:-(n + k), n - k:-(n + k)])

                if (n - k) != 0:
                    allright.append(Data_array[k + n:-(n - k), n:-n])
                    allup.append(Data_array[n:-n, k + n:-(n - k)])
                    alldiagleftdown.append(Data_array[k + n:-(n - k), n - k:-(n + k)])
                    alldiagleftup.append(Data_array[k + n:-(n - k), k + n:-(n - k)])
                    alldiagrightup.append(Data_array[n - k:-(n + k), k + n:-(n - k)])

                else:  # correct python array slicing at the end :   a[n:0]  would mean a[n:]

                    allright.append(Data_array[k + n:, n:-n])
                    allup.append(Data_array[n:-n, k + n:])
                    alldiagleftdown.append(Data_array[k + n:, n - k:-(n + k)])
                    alldiagleftup.append(Data_array[k + n:, k + n:])
                    alldiagrightup.append(Data_array[n - k:-(n + k), k + n:])

            return shift_zero, \
                allleft, allright, alldown, allup, \
                alldiagleftdown, alldiagrightup, alldiagrightdown, alldiagleftup

        else:
            shift_zero = Data_array[n:-n, n:-n]

            allleft = []
            allright = []
            allup = []
            alldown = []

            allleft.append(Data_array[:-2 * n, n:-n])
            alldown.append(Data_array[n:-n, :-2 * n])

            for k in np.arange(1, n)[::-1]:
                allleft.append(Data_array[n - k:-(n + k), n:-n])
                allright.append(Data_array[k + n:-(n - k), n:-n])
                alldown.append(Data_array[n:-n, n - k:-(n + k)])
                allup.append(Data_array[n:-n, k + n:-(n - k)])

            allright.append(Data_array[2 * n:, n:-n])
            allup.append(Data_array[n:-n, 2 * n:])

            return shift_zero, allleft, allright, alldown, allup

    elif dimensions == 1:
        shift_zero = Data_array[n:-n]
        allleft = []
        allright = []
        allleft.append(Data_array[:-2 * n])
        for k in np.arange(1, n)[::-1]:
            allright.append(Data_array[k + n:-(n - k)])
            allleft.append(Data_array[n - k:-(n + k)])
        allright.append(Data_array[2 * n:])

        return shift_zero, allleft, allright
	
def minmax(D_array, center, boxsize,
           framedim=(2048, 2048), withmaxpos=0):
    """
    extract min and max from a 2d array

    # TODO: replace by scipy.ndimage.extrema
    # framedim = from dictionary of CCDs
    D_array shape is flip(framedim)
    """
    if not isinstance(boxsize, int):
        boxsize = (boxsize, boxsize)
#    print "framedim in minmax", framedim
#    print "D_array.shape", D_array.shape
    halfbox = int(boxsize / 2)

#    xc, yc = center
#    imin, imax, jmin, jmax = max(0, yc - halfbox), \
#                            min(yc + halfbox, framedim[0]), \
#                            max(0, xc - halfbox), \
#                            min(framedim[1], xc + halfbox)
#
#    print "imin, imax, jmin, jmax", imin, imax, jmin, jmax
    framedim = framedim[1], framedim[0]
    imin, imax, jmin, jmax = getindices2cropArray(center, boxsize, framedim, flipxy=1)

    print "imin, imax, jmin, jmax", imin, imax, jmin, jmax

    fulldata = D_array
    array_short = fulldata[imin:imax, jmin:jmax]

    if withmaxpos:
        return [np.amin(array_short), np.amax(array_short)], \
                ndimage.maximum_position(array_short) + np.array([imin, jmin])
    else:
        return [np.amin(array_short), np.amax(array_short)]
	
def getindices2cropArray(center, boxsizeROI, arrayshape, flipxy=0):
    xpic, ypic = center
    if flipxy:
        ypic, xpic = center

    xpic, ypic = int(xpic), int(ypic)



    if isinstance(boxsizeROI, int):
        boxsizex, boxsizey = boxsizeROI, boxsizeROI
    else:
        boxsizex, boxsizey = boxsizeROI

    x1 = np.maximum(0, xpic - boxsizex)
    x2 = np.minimum(arrayshape[0], xpic + boxsizex)
    y1 = np.maximum(0, ypic - boxsizey)
    y2 = np.minimum(arrayshape[1], ypic + boxsizey)

    imin, imax, jmin, jmax = y1, y2, x1, x2

    return imin, imax, jmin, jmax
    
def removeClosePoints(X, Y, dist_tolerance=0.5):
    """
    remove very close spots within dist_tolerance (cartesian distance)
    """
    coord = np.array([X, Y]).T
    dist_tab = calcdistancetab(coord, coord)

    close_pos = np.where(dist_tab < dist_tolerance)

    i, j = close_pos

#    print "close_pos", close_pos
#    print "i", i
#    print "j", j

    dict_sets = Set_dict_frompairs(np.array([i, j]).T, verbose=0)[0]

#    print "dict_sets", dict_sets

    toremove = []
    for val in dict_sets.values():
        if len(val) > 1:
            toremove += val[1:]

    tokeep = list(set(range(len(X))) - set(toremove))

    return X[tokeep], Y[tokeep], tokeep

def calcdistancetab(listpoints1, listpoints2):
    data1 = np.array(listpoints1)
    data2 = np.array(listpoints2)
    xdata1 = data1[:, 0]
    ydata1 = data1[:, 1]
    xdata2 = data2[:, 0]
    ydata2 = data2[:, 1]
    deltax = xdata1 - np.reshape(xdata2, (len(xdata2), 1))
    deltay = ydata1 - np.reshape(ydata2, (len(ydata2), 1))
    didist = np.sqrt(deltax ** 2 + deltay ** 2)

    return didist
    
def Set_dict_frompairs(pairs_index, verbose=0):
    """
    from association pairs of integers return dictionnary of associated integer

    example:
    array([[ 0,  1],[ 0,  2],[ 1,  2],[ 3,  7],[ 3, 11],[ 4,  9],[ 4, 14],[ 7, 11],[ 9, 14],
    [15, 31],[15, 47],[16, 33],[16, 50],[19, 39],[19, 59],[20, 41],[20, 62],[31, 47],[19,14],
    [39, 59],[33, 0]])
    ->
    {0: [0, 1, 2, 33, 16, 50],
    3: [3, 7, 11],
    4: [4, 9, 14, 59, 19, 39],
    15: [15, 47, 31],
    20: [20, 41, 62]}

    """
    if len(pairs_index) == 1:
        #print "pairs_index",pairs_index
        res_final = {}
        res_final[pairs_index[0][0]] = [pairs_index[0][0], pairs_index[0][1]]
        return res_final, res_final

    res_dict = {}

    for elem in set(np.ravel(np.array(pairs_index))):
        pairs = return_pair(elem, pairs_index)
        #print "\nelem ",elem," pairs ",pairs
        #print "res_dict",res_dict
        if len(pairs) > 1:
            classmembers = [elem] + pairs.tolist()
            set_min = set()
            for elem in classmembers:
                if elem in res_dict:
                    set_min = set_min.union(res_dict[elem])
                else:
                    set_min.add(elem)

                res_dict[elem] = set_min

    if verbose:
        print "res_dict", res_dict

    res_final = {}

    for key, val in res_dict.items():
        listval = list(val)
        smallest_index = min(listval)

        if smallest_index in res_final:
            res_final[smallest_index].append(key)
        else:
            res_final[smallest_index] = [key]

    return res_final, res_dict
    
def return_pair(n, pairs):
    """
    return array of integer that are in correspondence in pairs with integer n

    """
    Pairs = np.array(pairs)

    i, j = np.where(Pairs == n)
    j = (j + 1) % 2
    return Pairs[(i, j)]
    
def readoneimage_manycrops(filename, centers, boxsize,
                           framedim=(2048, 2048), offset=4096,
                            formatdata="uint16", fliprot="no"):
    """
    reads 1 image and extract many regions
    centered on center_pixel with xyboxsizedimensions in pixel unit
    
    centers   : list or array of [x,y]

    returns:
    arrayofdata: list of 2D array of intensity
    """

    if type(boxsize) in (type(5), type(5.)):
        xboxsize, yboxsize = int(boxsize), int(boxsize)
    else:
        xboxsize, yboxsize = boxsize

    fulldata = np.reshape(readoneimage(filename,
                framedim=framedim, offset=offset,
                formatdata=formatdata), framedim)  # (framedim[1],framedim[0]))

    if fliprot == "spe":
        fulldata = np.rot90(fulldata, k=1)

    elif fliprot == "vhr":
        fulldata = np.rot90(fulldata, k=3)
        fulldata = np.fliplr(fulldata)

    elif fliprot == "vhrdiamond":
        fulldata = np.rot90(fulldata, k=3)
        fulldata = np.fliplr(fulldata)
        framedim = framedim[1], framedim[0]

    elif fliprot == "VHR_Feb13":
        # TODO: fliplr useful?
        fulldata = np.fliplr(fulldata)

    elif fliprot == "frelon2":
        fulldata = np.flipud(fulldata)

    if type(boxsize) == type(5):
        boxsizex, boxsizey = boxsize, boxsize
    elif type(boxsize) == type((10, 20)):
        boxsizex, boxsizey = boxsize

    xpic, ypic = np.array(centers).T

    x1 = np.array(np.maximum(0, xpic - boxsizex), dtype=np.int)
    x2 = np.array(np.minimum(framedim[0], xpic + boxsizex), dtype=np.int)
    y1 = np.array(np.maximum(0, ypic - boxsizey), dtype=np.int)
    y2 = np.array(np.minimum(framedim[1], ypic + boxsizey), dtype=np.int)

    Data = []

    print "framedim in readoneimage_manycrops", framedim
    framedim = framedim[1], framedim[0]
    print "framedim in readoneimage_manycrops", framedim
    for center in centers:
        i1, i2, j1, j2 = getindices2cropArray(center, (boxsizex, boxsizey), framedim)
        print "i1, i2, j1, j2-----", i1, i2, j1, j2
        cropdata = fulldata[i1:i2, j1:j2]
        #print "cropdata.shape", cropdata.shape
        Data.append(cropdata)

        print np.amax(cropdata)


#
#        i1, i2, j1, j2 = getindices2cropArray(center, (boxsizex, boxsizey), framedim)
#
#        print "i1, i2, j1, j2", i1, i2, j1, j2
#        cropdata = fulldata[j1:j2, i1:i2]
#        #print "cropdata.shape", cropdata.shape
#        Data.append(cropdata)


    #print "(x1, y1, x2, y2)", (x1, y1, x2, y2)

#    for k, box in enumerate(zip(y1, y2, x1, x2)):
#        print "k, box", k, box
##        cropdata = fulldata[box[1]:box[3], box[0]:box[2]]
##        cropdata = fulldata[box[0]:box[1], box[2]:box[3]]
#        #print "cropdata.shape", cropdata.shape
##        Data.append(cropdata)
    return Data
#t = PeakSearch('/home/fengguo/Data/21Feb13/Si1g_5N/nomvt/S1gnomvt_0000_mar.tiff')

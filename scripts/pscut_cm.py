from scipy.ndimage import center_of_mass
from read_drs import event_generator, return_dtype
from numpy import fromstring, histogram, std, sign, uint16, array
from numpy import median, argmax


def get_spec(fname, roi_start, roi_width=180, nchannels=2, force_start=False,
             **kwargs):
    """return a sipm spectrum using the cm method as a cut.
    roi_start is the star of the region of interest, + roi_width channels
    ref_cm is the reference center of mass. if None then mean(cm) of all events
    will be calculated.
    dev_cm is the allowed deviation from ref_cm. if None then std(cm)
    of all events will be calculated.
    nchannels is the number of DRS channels with data. either 1 or 2"""

    st, wd = roi_start, roi_width
    my_dtype = return_dtype(nchannels)
    if not force_start:
        st, ref_cm, dev_cm = find_start(fname, roi_start, roi_width, nchannels,
                                        **kwargs)
    else:
        cmsarr = cms_(fname, roi_start, roi_width, nchannels)
        cmhist = histogram(cmsarr, bins=512)
        ref_cm = cmhist[1][argmax(cmhist[0])]
        dev_cm = dev_cm_(cmsarr)

    with open(fname, 'r') as f:
        gen = (fromstring(event, my_dtype)[0][5]
               for event in event_generator(f, nchannels))
        specdata = [sum(event[st:st + wd]) for event in gen
                    if abs(center_of_mass(- event[st:st + wd])[0] - ref_cm)
                    < dev_cm]
    return histogram(specdata, bins=2048)


def get_darks(fname, protoevent, roi_start=9, ref_cm=89.5, dev_cm=None,
              roi_width=180, nchannels=2):
    """dev_cm should be the same as in the normal cut, so it gets calculated
    in the same way"""
    print fname
    st, wd = roi_start, roi_width
    my_dtype = return_dtype(2)
    protoevent = array(protoevent, dtype=uint16)
    if dev_cm is None:
        spom, spam, dev_cm = find_start(fname, 340, roi_width, nchannels)

    with open(fname, 'r') as f:
        gen = (fromstring(event, my_dtype)[0][5]
               for event in event_generator(f, 2))
        specdata = [sum((event + protoevent)[st:st + wd]) for event in gen
                    if abs(center_of_mass(-
                           (event + protoevent)[st:st + wd])[0] - ref_cm)
                    < dev_cm]
    return histogram(specdata, bins=2048)


def find_start(fname, roi_start, roi_width, nchannels, res_=25):
    cmsarr = cms_(fname, roi_start, roi_width, nchannels)
    cmhist = histogram(cmsarr, bins=512)
    cms = cmhist[1][argmax(cmhist[0])]
    #cms = median(cmsarr)
    res = (cms - (roi_width / 2. - 0.5)) * roi_width
    if res < res_:
        roi_start = roi_start + round(res, 0)
    else:
        return roi_start, cms, dev_cm_(cmsarr)
    it = 0
    while abs(res) > 1:
        cmsarr = cms_(fname, roi_start, roi_width, nchannels)
        cmhist = histogram(cmsarr, bins=512)
        cms = cmhist[1][argmax(cmhist[0])]
        #cms = median(cmsarr)
        res = (cms - (roi_width / 2. - 0.5)) * roi_width
        if res < 25:
            roi_start = roi_start + round(res, 0)
        else:
                return roi_start, cms, dev_cm_(cmsarr)
        roi_start = roi_start + round(res, 0)
        it += 1
        if it == 5:
            roi_start -= 1*sign(res)
            break
    print ' '.join(('found start in %i iterations:' % it, str(roi_start)))
    return roi_start, cms, dev_cm_(cmsarr)


def cms_(fname, roi_start, roi_width, nchannels=2):
    st, wd = roi_start, roi_width
    my_dtype = return_dtype(nchannels)
    with open(fname, 'r') as f:
        gen = (fromstring(event, my_dtype)[0][5]
               for event in event_generator(f, nchannels))
        cms = [center_of_mass(-event[st:st + wd])[0] for event in gen]
    return cms


def dev_cm_(cms):
    return std(cms)


def cms_dark(fname, roi_start, roi_width, protoevent):
    st, wd = roi_start, roi_width
    my_dtype = return_dtype(2)
    protoevent = array(protoevent, dtype=uint16)
    with open(fname, 'r') as f:
        gen = (fromstring(event, my_dtype)[0][5]
               for event in event_generator(f, 2))
        cms = [center_of_mass(-(event + protoevent)[st:st + wd])[0]
               for event in gen]
    return cms


def moveevent(event, x, peakstart=340):
    head = list(event[peakstart - x:])
    tail = list(event[:peakstart - x])
    return array(head + tail)

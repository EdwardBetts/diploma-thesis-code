from scipy.ndimage import center_of_mass
from read_drs import event_generator, return_dtype
from numpy import fromstring, histogram, mean, std


def get_spec(fname, roi_start, ref_cm=None, dev_cm=None, roi_width=180,
             nchannels=2):
    """return a sipm spectrum using the cm method as a cut.
    roi_start is the star of the region of interest, + roi_width channels
    ref_cm is the reference center of mass. if None then mean(cm) of all events
    will be calculated.
    dev_cm is the allowed deviation from ref_cm. if None then std(cm) / 2
    of all events will be calculated.
    nchannels is the number of DRS channels with data. either 1 or 2"""

    st, wd = roi_start, roi_width
    my_dtype = return_dtype(2)
    if ref_cm is None:
        cms = cms_(fname, roi_start, roi_width)
        ref_cm = mean(cms)
    if dev_cm is None:
        dev_cm = dev_cm_(cms)

    with open(fname, 'r') as f:
        gen = (fromstring(event, my_dtype)[0][5]
               for event in event_generator(f, 2))
        specdata = [sum(event[st:st + wd]) for event in gen
                    if abs(center_of_mass(event[st:st+wd])[0]-ref_cm) < dev_cm]
    return histogram(specdata, bins=2048)


def get_darks(fname, roi_start, protoevent, ref_cm=None, dev_cm=None,
              roi_width=180, nchannels=2):

    st, wd = roi_start, roi_width
    my_dtype = return_dtype(2)
    if ref_cm is None:
        cms = cms_dark(fname, roi_start, roi_width, protoevent)
        ref_cm = mean(cms)
    if dev_cm is None:
        dev_cm = dev_cm_(cms)

    with open(fname, 'r') as f:
        gen = (fromstring(event, my_dtype)[0][5]
               for event in event_generator(f, 2))
        specdata = [sum((event + protoevent)[st:st + wd]) for event in gen
                    if abs(center_of_mass(event[st:st+wd])[0]-ref_cm) < dev_cm]
    return histogram(specdata, bins=2048)


def cms_(fname, roi_start, roi_width):
    st, wd = roi_start, roi_width
    my_dtype = return_dtype(2)
    with open(fname, 'r') as f:
        gen = (fromstring(event, my_dtype)[0][5]
               for event in event_generator(f, 2))
        cms = [center_of_mass(event[st:st + wd])[0] for event in gen]
    return cms


def dev_cm_(cms):
    return std(cms) / 2


def cms_dark(fname, roi_start, roi_width, protoevent):
    st, wd = roi_start, roi_width
    my_dtype = return_dtype(2)
    with open(fname, 'r') as f:
        gen = (fromstring(event, my_dtype)[0][5]
               for event in event_generator(f, 2))
        cms = [center_of_mass((event + protoevent)[st:st + wd])[0]
               for event in gen]
    return cms

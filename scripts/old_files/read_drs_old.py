import numpy as np
from scipy.signal import iirfilter, filtfilt
from fsum import fsum, fmin
from eventtools import to_units, extract_events
from contextlib import contextmanager


#reads binary files from DRS4 Chip software version 4


def process_file(filename, filter_params=0.1,
                 p_w=(355, 530), cutoff=32600, nchannels=1,
                 f_type='lowpass', **args):
    #reads each event, filters it and return filtered minima
    #p_w: peak window
    #cutoff: geatest pretrigger min value

    b, a = iirfilter(1, filter_params, btype=f_type, **args)
    my_dtype = return_dtype(nchannels)
    with open(filename, 'rb') as f:
        fil_min = [min(fevent[p_w[0]:p_w[1]])
                   for fevent in (filtfilt(b, a, np.fromstring(event,
                                  my_dtype)[0][5])
                                  for event in event_generator(f, nchannels))
                   if min(fevent[:p_w[0]]) > cutoff
                   and max(fevent[:10]) < 350]
    return fil_min


def read_events(filename, n_events, chnl=5, nchannels=1):
    #drs_chnls: 1->5, 2->7, 3->9, 4->11; time->3
    with trace_gen(filename, nchannels, chnl) as gen:
        events = [gen.next() for i in range(n_events)]
    return events


def event_generator(f, nchannels=1):
    rec = True
    while rec:
        rec = f.read(4120 + nchannels * 2052)
        if rec:
            yield rec


def return_dtype(n):
    #binary structure can be seen in
    #http://www.psi.ch/drs/DocumentationEN/manual_rev40.pdf
    #define dtype
    if n == 1:
        dtype = np.dtype([('header', np.str_, 4), ('serial number', np.int32),
                         ('date', np.int16, 8), ('time', np.float32, 1024),
                         ('channelheader', np.str_, 4),
                         ('voltage', np.uint16, 1024)])
    elif n == 2:
        dtype = np.dtype([('header', np.str_, 4), ('serial number', np.int32),
                         ('date', np.int16, 8), ('time', np.float32, 1024),
                         ('channelheader', np.str_, 4),
                         ('voltage', np.uint16, 1024),
                         ('channelheader_2', np.str_, 4),
                         ('voltage_2', np.uint16, 1024)])
    return dtype


@contextmanager
def trace_gen(filename, nchannels=1, chnl=5):
    """generator which yields event from binary drs4 file.
    usage: with trace_gen() as gen: events = [event for event in gen]"""
    f = open(filename, 'rb')
    dtype = return_dtype(nchannels)
    gen = (np.fromstring(event, dtype)[0][chnl]
           for event in event_generator(f, nchannels))
    try:
        yield gen
    finally:
        f.close()


def scatter(filename, thrs, nchannels=2, chnl=5):
    b, a = iirfilter(1, (0.002, 0.05))
    with trace_gen(filename, nchannels, chnl) as gen:
        data = [(sum(event[thrs:thrs + 180]),
                 min(filtfilt(b, a, event)[thrs:thrs + 180]) - 200)
                for event in gen if max(filtfilt(b, a, event)[:20]) < 400]

    int_data = zip(*data)[0]
    min_data = zip(*data)[1]
    return int_data, min_data


def lowpass_scatter(filename, thrs, nchannels=2):
    b, a = iirfilter(1, 0.05, btype='lowpass')
    with trace_gen(filename, nchannels) as gen:
        data = [(sum(event[thrs:thrs + 180]),
                 min(filtfilt(b, a, event)[thrs:thrs + 180]) - 200)
                for event in gen]

    int_data = zip(*data)[0]
    min_data = zip(*data)[1]
    return int_data, min_data


def int_spec(filename, win, nchannels=2, chnl=5, fortran=True, range=None):
    with trace_gen(filename, nchannels, chnl) as gen:
        if fortran:
            int_data = [fsum(event[win[0]:win[1]]) for event in gen]
        else:
            int_data = [event[win[0]:win[1]].sum() for event in gen]

    return np.histogram(int_data, bins=2048, range=range)


def min_spec(filename, win, nchannels=2, chnl=5, range=None):
    with trace_gen(filename, nchannels, chnl) as gen:
        min_data = [fmin(event[win[0]:win[1]]) for event in gen]

    return np.histogram(min_data, bins=2048, range=range)


def base_test(filename, win, nchannels=2, chnl=5, **kwargs):
    with trace_gen(filename, nchannels, chnl) as gen:
        min_data = [fsum(to_units(event[win[0]:win[1]], **kwargs))
                    for event in gen]

    return np.histogram(min_data, bins=2048)


def fact_dark_spec(filename, thres, nchannels=2, chnl=5, base=0, **kwargs):
    with trace_gen(filename, nchannels, chnl) as gen:
        int_data = [fsum(event) for trace in gen
                    for event in extract_events(trace, base, thres)]

    return np.histogram(int_data, bins=2048)

import numpy as np
from fsum import fsum
from contextlib import contextmanager
from eventtools import to_units, extract_events
import datetime as dt
# read drs4v5 files


def return_header(f, nchannels=1):
    # tbw: time bin width
    if nchannels == 1:
        dtype = np.dtype([('time header', np.str_, 4),
                          ('serial number', np.int16, 2),
                          ('c1 header', np.str_, 4),
                          ('c1 tbw', np.float32, 1024)])
    elif nchannels == 2:
        dtype = np.dtype([('time header', np.str_, 4),
                          ('serial number', np.int16, 2),
                          ('c1 header', np.str_, 4),
                          ('c1 tbw', np.float32, 1024),
                          ('c2 header', np.str_, 4),
                          ('c2 tbw', np.float32, 1024)])

    return np.fromstring(f.read(8 + nchannels * 4100), dtype)


def return_dtype(n):
    # binary structure can be seen in
    # http://www.psi.ch/drs/DocumentationEN/manual_rev50.pdf
    # define dtype
    if n == 1:
        dtype = np.dtype([('header', np.str_, 4), ('event number', np.int32),
                          ('date', np.int16, 7), ('range', np.int16),
                          ('serial number', np.int8, 6),
                          ('trigger cell', np.int16),
                          ('c1 header', np.str_, 4),
                          ('c1 voltage', np.uint16, 1024)])
    elif n == 2:
        dtype = np.dtype([('header', np.str_, 4), ('event number', np.int32),
                          ('date', np.int16, 7), ('range', np.int16),
                          ('serial number', np.int8, 6),
                          ('trigger cell', np.int16),
                          ('c1 header', np.str_, 4),
                          ('c1 voltage', np.uint16, 1024),
                          ('c2 header', np.str_, 4),
                          ('c2 voltage', np.uint16, 1024)])
    return dtype


def event_generator(f, nchannels=1, n=1):
    rec = True
    while rec:
        rec = f.read(n * (32 + nchannels * 2052))
        if rec:
            yield rec


def calc_time(event, bin_width, trigg_cell):
    return [fsum([bin_width[(j + trigg_cell) % 1024] for j in range(i)])
            for i, channel in enumerate(event)]


channels = ['c1', 'c2']


@contextmanager
def trace_gen(filename, nchannels=1, channel='c1', base_offset=0):
    if channel not in channels:
        raise KeyError

    f = open(filename, 'rb')
    header = return_header(f, nchannels)
    dtype = return_dtype(nchannels)
    chan_string = ''.join((channel, ' voltage'))
    gen = (to_units(np.fromstring(event, dtype)[chan_string][0], base_offset)
           for event in event_generator(f, nchannels))
    try:
        yield gen
    finally:
        f.close()


@contextmanager
def two_channel_trace_gen(filename, channels=('c1', 'c2'), base_offset=0):
    f = open(filename, 'rb')
    header = return_header(f, 2)
    dtype = return_dtype(2)
    chan_strings = [''.join((ch, ' voltage')) for ch in channels]
    gen = ([to_units(np.fromstring(event, dtype)[ch][0], base_offset)
            for ch in chan_strings]
           for event in event_generator(f, 2))
    try:
        yield gen
    finally:
        f.close()


@contextmanager
def time_trace_gen(filename, nchannels=1, channel='c1', base_offset=0):
    if channel not in channels:
        raise KeyError

    f = open(filename, 'rb')
    header = return_header(f, nchannels)
    dtype = return_dtype(nchannels)
    chan_string = ''.join((channel, ' voltage'))
    gen = ((dt.datetime(*raw['date']), to_units(raw[chan_string],
           base_offset))
           for event in event_generator(f, nchannels)
           for raw in np.fromstring(event, dtype))
    try:
        yield gen
    finally:
        f.close()


def read_events(filename, nevents, nchannels=1, channel='c1', base_offset=0):
    with trace_gen(filename, nchannels, channel, base_offset) as gen:
        return [gen.next() for i in range(nevents)]


def base_test(filename, win, nchannels=1, chnl='c1', base_offset=0):
    if chnl not in channels:
        raise KeyError

    with trace_gen(filename, nchannels, chnl, base_offset) as gen:
        min_data = [fsum(event[win[0]:win[1]]) for event in gen]

    return np.histogram(min_data, bins=2048)


def fact_dark_spec(filename, thres, nchannels=1, chnl='c1', base=0, **kwargs):
    with trace_gen(filename, nchannels, chnl, base) as gen:
        int_data = [fsum(event) for trace in gen
                    for event in extract_events(trace, thres, **kwargs)]

    return np.histogram(int_data, bins=2048)

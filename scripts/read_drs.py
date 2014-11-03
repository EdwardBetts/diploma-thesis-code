import numpy as np
from fsum import fsum
from contextlib import contextmanager
from eventtools import to_units

#read drs4v5 files


def return_header(f, nchannels=1):
    #tbw: time bin width
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
    #binary structure can be seen in
    #http://www.psi.ch/drs/DocumentationEN/manual_rev50.pdf
    #define dtype
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


def event_generator(f, nchannels=1):
    rec = True
    while rec:
        rec = f.read(32 + nchannels * 2052)
        if rec:
            yield rec


def calc_time(event, bin_width, trigg_cell):
    return [fsum([bin_width[(j + trigg_cell) % 1024] for j in range(i)])
            for i, channel in enumerate(event)]


channels = ['c1', 'c2']


@contextmanager
def trace_gen(filename, nchannels=1, channel='c1'):
    if not channel in channels:
        raise KeyError

    f = open(filename, 'rb')
    header = return_header(f, nchannels)
    dtype = return_dtype(nchannels)
    chan_string = ''.join((channel, ' voltage'))
    gen = (to_units(np.fromstring(event, dtype)[chan_string][0])
           for event in event_generator(f, nchannels))
    try:
        yield gen
    finally:
        f.close()


def read_events(filename, nevents, nchannels=1, channel='c1'):
    with trace_gen(filename, nchannels, channel) as gen:
        return [gen.next() for i in range(nevents)]

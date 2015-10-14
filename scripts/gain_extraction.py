from fsum import arm
from read_drs import trace_gen
import numpy as np


def argrelmin(event, order=30):
    a = arm(event, order)
    return a[0][:a[1]]


def alt_spec_for_gain(filename, kernel, nchannels=1, chnl='c1'):
    with trace_gen(filename, nchannels, chnl) as gen:
        min_data = []
        for event in gen:
            arr = np.convolve(event, kernel, mode='valid')
            min_data.extend(arr[argrelmin(arr, order=20)])
    hist = np.histogram(min_data, bins=1024)
    s = np.log(hist[0])
    s[s < 0] = 0
    return (s, hist[1])

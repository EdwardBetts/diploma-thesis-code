import numpy as np
from fsum import *


def sliding_average(x, window=10):
    return np.convolve(x, np.ones((window,)) / window, mode='valid')


def to_units(x, base_offset=0):
    return x / 65536. - 0.5 - base_offset


def extract_events(trace, trigger_threshold, base_offset=0, old=False):
    #clip first and last 5 bins due to errors
    #readout of old drs4 board doesnt convert to physical units automatically
    if old:
        raw_event = to_units(trace[5:-5], base_offset=base_offset)
    else:
        raw_event = trace[5:-5] - base_offset

    #perform sliding average to not be less susceptiable to noise
    event = sliding_average(raw_event)

    #trigger on 1 pe event
    events = []
    #start at 30, to search for arrival time correctly
    start = 30
    #integration window
    i_win = 35
    while start < 950:
        trigger = trig_ind_pos(trigger_threshold, event, start)
        if 950 > trigger > start:
            mx, mxi = max_and_ind(event, trigger + 5)
            halfmax = mx / 2
            arrival_time = trig_ind_pos(halfmax, event, mxi - 30)
            distance = mxi - arrival_time
            if not arrival_time:
                distance = -1
            if 0 < distance < 20:
                #+4 is the offset produced by the sliding average
                events.append(raw_event[arrival_time+4:arrival_time+4 + i_win])
                start = arrival_time + i_win
            else:
                start = mxi + 15
        else:
            start = 950

    return events

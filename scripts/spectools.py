from smooth import smooth
from numpy import log, argmin, sqrt, r_, where, mean, argmax, histogram, arange
from scipy.signal import find_peaks_cwt
from read_drs import base_test, trace_gen
from eventtools import sliding_average
from fit import find_peak_data


def find_params(filename, pos_pol=True, **kwargs):
    base = find_base(filename, **kwargs)
    threshold, n_events = find_threshold(
        filename, base_offset=base, pos_pol=pos_pol, **kwargs)
    return base, threshold, n_events


def find_base(filename, **kwargs):
    spec = base_test(filename, (300, 400), **kwargs)
    # 4 from the 9 bins difference due to filtering
    return spec[1][argmax(sliding_average(spec[0])) + 4] / 100


def find_threshold(filename, base_offset=0, pos_pol=True, **kws):
    with trace_gen(filename, base_offset=base_offset, **kws) as gen:
        if pos_pol:
            min_data = [event[100:-100].max() for event in gen]
        else:
            min_data = [event[100:-100].min() for event in gen]
    hist = histogram(min_data, bins=2048)
    # second = 2 if pos_pol else -2
    # index = argrelmax(hist[0][hist[0] != 0], order=5)[0][second]
    second = 1 if pos_pol else -1
    try:
        index = find_peaks_cwt(hist[0][hist[0] != 0], arange(10, 40))[second]
    except IndexError:
        index = find_peaks_cwt(hist[0][hist[0] != 0], arange(10, 40))[0]
    return hist[1][hist[0] != 0][index] * 3 / 4, len(min_data)


def get_gain(hist, reverse=False):
    # finds gain in arbitrary units from the difference of the first two peaks
    if reverse:
        hist = [line[::-1] for line in hist]
    try:
        params = find_peak_data(hist[0], 2)
    except IndexError:
        return -1
    return hist[1][int(params[4])] - hist[1][int(params[1])]


def get_n_mean(data, cutoff):
    pedestal = float(sum(data[cutoff:]))
    all_ = sum(data)
    return - log(pedestal / all_)


def get_cutoff(data, guess=None, window=50):
    if not guess:
        guess = get_guess(data)
    smoothdata = smooth(data)
    amin = argmin(smoothdata[guess - window: guess + window])
    return amin - 5.5 + guess - window


def get_co_exp(data, guess, window=30):
    co = guess
    for i in range(10):
        co_last = co
        co = get_cutoff(data, co, window)
        if abs(co_last - co) < 1:
            return co
    print 'no min found'
    return guess


def get_nmean_errors(data, cutoff, swin=7):
    pedestal = float(sum(data[cutoff:]))
    all_ = sum(data)
    n_mean = -log(pedestal / all_)
    # statistical error
    sterr = sqrt(pedestal) / pedestal
    # systematic error:
    syserr = max([abs(n_mean - get_n_mean(data, cutoff - swin + 2 * swin * i))
                  for i in range(2)])
    return n_mean, sterr, syserr


def get_guess(spec, window_len=11, n=2, m=20, bins=2048):
    sdata = smooth(spec[::-1], window_len)

    # find values which are higher than neighbors
    test = r_[True, sdata[1:] > sdata[:-1]] & r_[sdata[:-1] > sdata[1:], True]

    # find values which are higher than the next m neighbors
    newInds = where(test)[0]
    values = [i - 5 for i in newInds[newInds > m] if sdata[i] ==
              max(sdata[i - m:i + m + 1]) and sdata[i] > 1.0]
    mu_values = [bins - i for i in values[n - 2:n]]
    return mean(mu_values)


def getnm(hist):
    guess = get_guess(hist)
    cutoff = get_cutoff(hist, guess)
    return get_nmean_errors(hist, cutoff)

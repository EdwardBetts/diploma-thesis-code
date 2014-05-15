from smooth import smooth
from numpy import log, argmin, sqrt, r_, where, mean


def get_n_mean(data, cutoff):
    pedestal = float(sum(data[cutoff:]))
    all_ = sum(data)
    return - log(pedestal / all_)


def get_cutoff(data, guess, window=50):
    smoothdata = smooth(data)
    amin = argmin(smoothdata[guess - window: guess + window])
    return amin - 5.5 + guess - window


def get_nmean_errors(data, cutoff, swin=7):
    pedestal = float(sum(data[cutoff:]))
    all_ = sum(data)
    n_mean = -log(pedestal / all_)
    #statistical error
    sterr = sqrt(pedestal) / pedestal
    #systematic error:
    syserr = max([abs(n_mean - get_n_mean(data, cutoff - swin + 2 * swin * i))
                 for i in range(2)])
    return n_mean, sterr, syserr


def get_guess(spec, window_len=11, n=2, m=20, bins=2048):
    sdata = smooth(spec[::-1], window_len)

    #find values which are higher than neighbors
    test = r_[True, sdata[1:] > sdata[:-1]] & r_[sdata[:-1] > sdata[1:], True]

    #find values which are higher than the next m neighbors
    newInds = where(test)[0]
    values = [i - 5 for i in newInds[newInds > m] if sdata[i] ==
              max(sdata[i-m:i+m + 1]) and sdata[i] > 1.0]
    mu_values = [bins - i for i in values[n - 2:n]]
    return mean(mu_values)

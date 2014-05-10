from numpy import histogram, percentile
from pde_int import get_guess, get_cutoff, get_n_mean
from Fit import odr_gauss


def return_ps_param(intd, mind, lineparams=None, **kwargs):
    if lineparams is None:
        a, b = get_line(intd, mind, **kwargs)
    else:
        a, b = lineparams
    return [(mind[i] / evnt - (b / evnt)) / a for i, evnt in enumerate(intd)]


def get_line(intd, mind, n=None, **kwargs):
    int_hist = histogram(intd, bins=2048)
    if n is None:
        int_mn = get_n_mean(int_hist[0], get_cutoff(int_hist[0],
                            get_guess(int_hist[0])))
        n = decide_n(int_mn)
    temp = odr_gauss(n, int_hist[0][::-1])
    x1, x2 = coord_from_fit(int_hist, temp, **kwargs)
    min_hist = histogram(mind, bins=2048)
    temp = odr_gauss(n, min_hist[0][::-1])
    y1, y2 = coord_from_fit(min_hist, temp, **kwargs)

    a = (y1 - y2) / (x1 - x2)
    b = y1 - a*x1
    print 'b = ' + str(b)
    print 'a = ' + str(a)
    return a, b


def decide_n(n):
    return int(n) + 6


def coord_from_fit(hist, fitdat, peaknum=4):
    pedx = hist[1][2048 - fitdat[0][4] - 1]
    peaknumx = hist[1][2048 - fitdat[0][(peaknum - 1) * 3 + 1] - 1]
    return pedx, peaknumx


def get_hist(intd, mind, lineparams=None, percs=(25, 57), **kwargs):
    ps_param = return_ps_param(intd, mind, lineparams, **kwargs)
    cuts = [percentile(ps_param, i) for i in percs]
    spec = [event for i, event in enumerate(mind)
            if ps_param[i] > cuts[0] and ps_param[i] < cuts[1]]
    return histogram(spec, bins=2048)

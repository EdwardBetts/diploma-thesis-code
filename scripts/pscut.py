from numpy import histogram, percentile
from pde_int import get_guess, get_cutoff, get_n_mean
from Fit import odr_gauss


def return_ps_param(intd, mind, b=None, **args):
    if b is None:
        a, b = get_line(intd, mind, **args)
    return [mind[i] / ievent - (b / ievent) for i, ievent in enumerate(intd)]


def get_line(intd, mind, n=None, **args):
    int_hist = histogram(intd, bins=2048)
    if n is None:
        int_mn = get_n_mean(int_hist[0], get_cutoff(int_hist[0],
                            get_guess(int_hist[0])))
        n = decide_n(int_mn)
    temp = odr_gauss(n, int_hist[0][::-1])
    x1, x2 = coord_from_fit(int_hist, temp, **args)
    min_hist = histogram(mind, bins=2048)
    min_mn = get_n_mean(min_hist[0], get_cutoff(min_hist[0],
                        get_guess(min_hist[0])))
    temp = odr_gauss(decide_n(min_mn), min_hist[0][::-1])
    y1, y2 = coord_from_fit(min_hist, temp, **args)

    a = (y1 - y2) / (x1 - x2)
    b = y1 - a*x1
    print 'b = ' + str(b)
    return a, b


def decide_n(n):
    return int(n) + 6


def coord_from_fit(hist, fitdat, peaknum=4):
    pedx = hist[1][2048 - fitdat[0][4] - 1]
    peaknumx = hist[1][2048 - fitdat[0][(peaknum - 1) * 3 + 1] - 1]
    return pedx, peaknumx


def get_hist(intd, mind, b=None, **args):
    ps_param = return_ps_param(intd, mind, b, **args)
    cuts = [percentile(ps_param, i) for i in (25, 57)]
    spec = [event for i, event in enumerate(mind)
            if ps_param[i] > cuts[0] and ps_param[i] < cuts[1]]
    return histogram(spec, bins=2048)

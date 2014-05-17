from read_drs import event_generator, return_dtype
from numpy import fromstring, sum, histogram, exp, sqrt
from pde_int import get_guess, get_cutoff, get_n_mean, get_nmean_errors
from Fit import getFitData
from scipy.signal import iirfilter, filtfilt


def darks(filename, thrs_1, thrs_2=None, nchannels=2):

    print 'getting dark count spectum'
    if thrs_2 is None:
        thrs_2 = thrs_1
    my_dtype = return_dtype(nchannels)

    with open(filename, 'r') as f:
        gen = (fromstring(event, my_dtype)[0][5]
               for event in event_generator(f, nchannels))
        traces = []

        for event in gen:
            for i, j in enumerate(event[40:100]):
                if event[i+40] < thrs_1 and event[i+60] < thrs_2:
                    traces.append(event[i+20:i+200])
                    break

    return [sm for sm in (sum(trace) for trace in traces) if sm > 0]


def peaks(filename, int_limits, nchannels=2):

    print 'getting peak spectrum'
    my_dtype = return_dtype(nchannels)
    with open(filename, 'r') as f:
        gen = (fromstring(event, my_dtype)[0][5]
               for event in event_generator(f, nchannels))
        peaks = [sum(event[int_limits[0]:int_limits[1]]) for event in gen]

    return peaks


def get_spectra(filename, int_limits, thrs=None, nchannels=2):
    if thrs is None:
        thrs = find_threshold(filename, int_limits)

    _darks = darks(filename, thrs_1=thrs, nchannels=nchannels)
    _peaks = peaks(filename, int_limits=int_limits, nchannels=nchannels)

    mx = max(max(_darks), max(_peaks))
    mn = min(min(_darks), min(_peaks))

    dark_hist = histogram(_darks, bins=2048, range=(mn, mx))
    peak_hist = histogram(_peaks, bins=2048, range=(mn, mx))
    sub_hist = peak_hist[0] - dark_hist[0]

    return dark_hist, peak_hist, sub_hist


def get_xtalk(filename, int_limits, thrs_1=None):
    if thrs_1 is not None:
        dark, peak, sub = get_spectra(filename, int_limits, thrs_1)
    else:
        peak = peaks(filename, int_limits)
        temp_hist = histogram(peak, bins=2048)
        cutoff = 2027 - getFitData(temp_hist[0][::-1], 3)[4]
        thrs = temp_hist[1][cutoff] / 180
        dark = darks(filename, thrs)
        mx = max(max(dark), max(peak))
        mn = min(min(dark), min(peak))

        dark_hist = histogram(dark, bins=2048, range=(mn, mx))
        peak_hist = histogram(peak, bins=2048, range=(mn, mx))
        sub = peak_hist[0] - dark_hist[0]

    guess = get_guess(sub)
    cutoff = get_cutoff(sub, get_cutoff(sub, guess))
    nm = get_n_mean(sub, cutoff)

    pg0 = 1 - exp(- nm)
    pg1_p = 1 - exp(- nm) - nm * exp(- nm)
    cutoff_2 = get_cutoff(sub, get_cutoff(sub, get_guess(sub, n=3, m=33)))
    print cutoff_2
    pg1 = sum(sub[:cutoff_2]) / float(sum(sub))
    xtalk = (pg1 - pg1_p) / pg0
    print xtalk
    return xtalk, sub


def xtalk_single_spec(hist):
    nm = get_n_mean(hist, get_cutoff(hist, get_guess(hist)))
    pg1_p = 1 - exp(- nm) - nm * exp(- nm)
    pg0 = 1 - exp(- nm)
    cutoff = get_cutoff(hist, get_cutoff(hist, get_guess(hist, n=3, m=33)))
    pg1 = sum(hist[:cutoff]) / float(sum(hist))
    xtalk = (pg1 - pg1_p) / pg0
    print xtalk
    return xtalk


def xtalksspec_err(hist):
    nm, nmerr, spm = get_nmean_errors(hist, get_cutoff(hist, get_guess(hist)))
    pg1_p = 1 - exp(- nm) - nm * exp(- nm)
    pg1_perr = exp(- nm) * nm * nmerr
    pg0 = 1 - exp(- nm)
    pg0err = nmerr
    cutoff = get_cutoff(hist, get_cutoff(hist, get_guess(hist, n=3, m=33)))
    pg1 = sum(hist[:cutoff]) / float(sum(hist))
    pg1err = sqrt(sum(hist[:cutoff])) / float(sum(hist))
    xtalk = (pg1 - pg1_p) / pg0
    err = sqrt((1 / pg0 * pg1err)**2 + (1 / pg0 * pg1_perr)**2 +
               (xtalk / pg0 * pg0err)**2)
    print xtalk
    return xtalk, err


def test_get_spectra(filename, int_limits):
    peak = peaks(filename, int_limits)
    temp_hist = histogram(peak, bins=2048)
    cutoff = 2047 - getFitData(temp_hist[0][::-1], 3)[4]
    thrs = temp_hist[1][cutoff] / 180
    dark = darks(filename, thrs)
    mx = max(max(dark), max(peak))
    mn = min(min(dark), min(peak))

    dark_hist = histogram(dark, bins=2048, range=(mn, mx))
    peak_hist = histogram(peak, bins=2048, range=(mn, mx))
    sub = peak_hist[0] - dark_hist[0]
    return dark_hist, peak_hist, sub


def find_threshold(filename, int_limits, nchannels=2):
    b, a = iirfilter(1, 0.05, btype='lowpass')
    my_dtype = return_dtype(nchannels)

    with open(filename, 'r') as f:
        gen = (fromstring(event, my_dtype)[0][5]
               for event in event_generator(f, nchannels))
        peaks = [min(event[int_limits[0]:int_limits[1]])for event
                 in (filtfilt(b, a, eventt) for eventt in gen)]
    hist = histogram(peaks, bins=2048)
    idx = get_cutoff(hist[0], get_guess(hist[0]))
    threshold = hist[1][idx]
    print threshold
    return threshold

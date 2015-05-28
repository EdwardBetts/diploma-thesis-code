from read_drs import fact_dark_spec
from spectools import find_params, get_gain
from crosstalk import xtalk_dark_spec, xtalk_spec_lmfit
from scipy.optimize import curve_fit


def extract_xtalk_dcr_gain(filename, timebase=1./(2*10**9), pos_pol=True,
                           threshold=None):
    base, threshold_, n_events = find_params(filename, pos_pol=pos_pol)
    if threshold is None:
        threshold = threshold_
    dark_spec = fact_dark_spec(filename, threshold, base=base, pos_pol=pos_pol)
    if not pos_pol:
        dark_spec = [data[::-1] for data in dark_spec]
    try:
        xtalk, result = xtalk_spec_lmfit(dark_spec[0])
    except Exception:
        print 'nan'
        xtalk = float('nan')
    # dark count rate ist number of events / ontime
    # ontime = 1024 channels * timebase * number of events
    ontime = 1024. * timebase * n_events
    dcr = dark_spec[0].sum() / ontime

    # gain
    gain = get_gain(dark_spec)
    return xtalk, dcr, gain


# actually calc breakdown voltage
def calc_gain(df):
    def line(x, *p):
        a, b = p
        return a*x + b

    df = df.dropna()
    fit_data = curve_fit(line, df.index.values, df['gain'].abs().values,
                         (.05, 1))
    return -fit_data[0][1] / fit_data[0][0]

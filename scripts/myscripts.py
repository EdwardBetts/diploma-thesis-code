from read_drs import fact_dark_spec
from spectools import find_base
from crosstalk import xtalk_dark_spec


def extract_xtalk_dcr(filename, thres, nevents=80000, timebase=1./(2*10**9),
                      p_func='erlang'):
    dark_spec = fact_dark_spec(filename, thres, base=find_base(filename))
    xtalk, fit_data = xtalk_dark_spec(dark_spec[0], 4, m=70, prop_func=p_func)

    #dark count rate ist number of events / ontime
    #ontime = 1024 channels * timebase * number of events
    ontime = 1024. * timebase * nevents
    dcr = dark_spec[0].sum() / ontime
    return xtalk, dcr

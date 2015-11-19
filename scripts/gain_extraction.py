from fsum import argrelmin as arm, correlate
from read_drs import trace_gen
import numpy as np
from scipy.stats import scoreatpercentile
from scipy.signal import argrelmax
import lmfit


def argrelmin(event, order=30):
    a = arm(event, order)
    return a[0][:a[1]]


def alt_spec_for_gain(filename, kernel, nchannels=1, chnl='c1',
                      rng=(-1.5, .3)):
    with trace_gen(filename, nchannels, chnl) as gen:
        min_data = []
        for event in gen:
            arr = correlate(event, kernel)
            min_data.extend(arr[argrelmin(arr, order=20)])
    hist = np.histogram(min_data, bins=1024, range=rng)
    return correct_spec(hist)


def approx_gain(hist, padding=5000):
    # use log scale
    s = np.log(hist[0])
    s[s < 0] = 0
    hist = (s, hist[1])

    to_fft = np.lib.pad(hist[0], (padding, padding), 'minimum')
    timestep = hist[1][1] - hist[1][0]
    freqs = np.fft.fftfreq(len(to_fft), timestep)
    # fft'd
    fftd = np.abs(np.fft.rfft(to_fft))
    index = experimental_find_fft_max(fftd)
    return 1. / freqs[index]


def experimental_find_fft_max(data):
    indices = argrelmax(data)[0]
    bg_score = scoreatpercentile(data, 10)
    ratios = [((data[i] / bg_score), i) for i in indices]
    return max(ratios)[1]


def fit_find_gain(spec):
    # maxindices = find_likely_max(spec, np.arange(10, 100, 200))
    maxindices = find_alt_max(spec)
    gauss = lmfit.models.GaussianModel(prefix='g0') + \
        lmfit.models.GaussianModel(prefix='g1')
    p = build_params(spec, maxindices)
    result = gauss.fit(spec[0], params=p, x=spec[1], weights=np.sqrt(spec[0]),
                       fit_kws={'maxfev': 200})
    return abs(result.params['g1center'] - result.params['g0center'])


def build_params(spec, arr):
    """for the later fit"""
    p = lmfit.Parameters()
    sigma = .01
    for i, index in enumerate(arr):
        p.add('g{}center'.format(i), value=spec[1][index])
        p.add('g{}amplitude'.format(i),
              value=spec[0][index] * sigma * np.sqrt(2*np.pi))
        p.add('g{}sigma'.format(i), value=sigma)
    return p


def find_likely_max(spec, sizes):
    mxs = np.zeros(1000)
    for size in sizes:
        new_mxs = argrelmax(np.convolve(spec[0], np.ones(size), mode='same'),
                            order=size)[0]
        if 2 <= len(new_mxs) < len(mxs):
            mxs = new_mxs
    if not mxs.any():
        return [0, 0]
    return mxs[-2:]


def find_alt_max(spec):
    sizes = np.linspace(10, len(spec[0]) / 2, 30)
    indices = [argrelmax(np.convolve(spec[0], np.ones(size), mode='same'),
                         order=1)[0] for size in sizes]
    a = np.diff(map(len, indices))
    inds = [i+1 for i, val in enumerate(a) if abs(val) <= np.std(a)]
    result = indices[inds[0]][-2:]
    i = 0
    while len(result) < 2:
        i += 1
        result = indices[inds[0] - i][-2:]
    return result


def correct_spec(spec):
    bins = spec[1]
    return spec[0], np.array([np.mean((bins[i+1], bins[i]))
                             for i in range(len(bins) - 1)])


class DynamicHist(object):
    """docstring for DynamicHist"""

    def __init__(self, kernel, initdata=None, depth=None):
        super(DynamicHist, self).__init__()
        self.kernel = kernel
        self.counter = np.zeros(4, dtype=np.int)
        self.depth = depth
        if self.depth is None or self.depth.shape != (4,):
            self.depth = np.zeros(4) + 0.01
        self.means = initdata
        if self.means is None or self.means.shape != (4,):
            raise ValueError

    def insert_trace(self, event):
        arr = correlate(event, self.kernel)
        vals = arr[argrelmin(arr, order=20)]
        for val in vals:
            i = self.find_index(val)
            n_val = self.means[i] * (1 - self.depth[i]) + val * self.depth[i]
            if i == 0:
                diff = n_val - self.means[i]
                self.means += diff
            else:
                self.means[i] = n_val
            self.counter[i] += 1

    def find_index(self, val):
        diff = np.inf
        index = 0
        for i in range(len(self.means)):
            dist = abs(self.means[i] - val)
            if dist < diff:
                diff = dist
                index = i
        return index

    def return_gain(self):
        return abs(self.means[0] - self.means[1])

    @classmethod
    def from_spec(cls, spec):
        raise NotImplementedError

    @classmethod
    def from_hist(cls, histdata, kernel):
        dyn_hist = cls(kernel, histdata)
        return dyn_hist

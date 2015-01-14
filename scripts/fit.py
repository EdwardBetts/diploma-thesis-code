from numpy import arange, sum, where
from scipy import odr
from scipy.signal import find_peaks_cwt
from fitting_funcs import multigauss, spec_func, mod_erlang, tail_gauss
from scipy.optimize import curve_fit
from lmfit import Model
from math import sqrt, pi


def find_peak_data(data, n, spec=False):
    indices = find_peaks_cwt(data, arange(30, 80))
    values = [data[i] for i in indices]
    params = [val for i in range(n) for val in [values[i], indices[i], 30]]
    if spec:
        params.extend((0.2, 1))
    return params


def peaks_for_spec(data, n):
    indices = find_peaks_cwt(data, arange(30, 80))
    params = [val for i in range(n) for val in [indices[i], 30]]
    params.extend((0.2, 1, data[indices[0]]))
    return params


#lmfit, trying to find number of peaks
def params_for_spec_fit(data):
    indices = find_peaks_cwt(data, arange(30, 80))
    n = len(where(data[indices] > 5)[0])
    params = [(indices[i], 30) for i in range(n)]
    params.append((0.2, 1, data[indices[0]]))
    return n, params


def odr_fit(n, data, func, params=None, fixed_params=None, nopeaks=4, **kw):
    """#n is the number of gaussians to fit
    params is an array of parameters for each gaussian (A, mu, sigma)
    for fixed params:
    'A value of 0 fixes the parameter, a value > 0 makes the parameter free.'
    """

    if params is None:
        params = find_peak_data(data, n)

    #set up odr module
    model = odr.Model(func)
    odr_data = odr.Data(range(len(data)), data)
    fit = odr.ODR(odr_data, model, params, ifixb=fixed_params, **kw)

    #run the fit
    output = fit.run()

    return output.beta, output.sd_beta


def odr_gauss(n, data, params=None, fixed_params=None, nopeaks=4, **kw):

    def func(parameters, x):
        return multigauss(x, n, *parameters)

    return odr_fit(n, data, func, params, fixed_params, nopeaks, **kw)


def odr_spec(n, data, params, prop_func='erlang', fixed_params=None, **kw):

    def func(parameters, x):
        return spec_func(x, n, prop_func, *parameters)

    return odr_fit(n, data, func, params, fixed_params, **kw)


def cf_spec(n, data, params, prop_func='erlang'):

    def func(x, *parameters):
        return spec_func(x, n, prop_func, *parameters)

    return curve_fit(func, arange(len(data)), data, params)


def build_func(n, pfunc, gfunc):
    def func(x, mu, sigma, p, nu, t):
        return pfunc(n+1, p, nu) * gfunc(x, mu, sigma, t)
    return func


def fit_dark_spec(n, data, params):
    assert len(params) == n+1
    models = [Model(build_func(i, mod_erlang, tail_gauss),
                    prefix='f{}'.format(i))for i in range(n)]

    def amp(x, a):
        return a
    ampmodel = Model(amp)

    model = sum(models) * ampmodel

    p, nu, A = params[-1]
    for i, model_ in enumerate(models):
        mu, sigma = params[i]
        model.set_param_hint('f{}mu'.format(i), value=mu)
        model.set_param_hint('f{}sigma'.format(i), value=sigma)
        model.set_param_hint('f{}t'.format(i), value=sigma)

        if i:
            model.set_param_hint('f{}p'.format(i), value=p, expr='f0p')
            model.set_param_hint('f{}nu'.format(i), value=nu, expr='f0nu')
        else:
            model.set_param_hint('f{}p'.format(i), value=p)
            model.set_param_hint('f{}nu'.format(i), value=nu)
    model.set_param_hint('a', value=A*sigma*sqrt(2*pi))

    result = model.fit(data, x=arange(2048))

    return model, result

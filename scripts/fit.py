import smooth
import numpy
from numpy import mean, array, r_
from scipy import odr
from fitting_funcs import multigauss, spec_func


def getFitData(data, n, window_len=11, m=30):
    sdata = smooth.smooth(data, window_len)

    #find values which are higher than neighbors
    test = r_[True, sdata[1:] > sdata[:-1]] & r_[sdata[:-1] > sdata[1:], True]

    #find values which are higher than the next m neighbors
    newInds = numpy.where(test)[0]
    values = [i - 5 for i in newInds[newInds > m]
              if sdata[i] == max(sdata[i-m:i+m + 1]) and sdata[i] > 1.0]
    mu_values = values[:n]

    a_values = [sdata[i] for i in mu_values]

    #shape output in correct way: A,mu,sigma
    return array([[a_values[i], mu_values[i], 30] for i in range(n)]).flatten()


def get_first_in_range(data, approx, w_len, window_len=11):
    #find minima
    test = r_[True, data[1:] < data[:-1]] & r_[data[:-1] < data[1:], True]

    inds = numpy.where(test)[0]
    #why next line?
    #inds = [i for i in inds if i > w_len]
    values = [i for i in inds if data[i] == min(data[i-w_len:i+w_len + 1])]

    if not values == []:
        return values[0] - 5.5 + approx - w_len
    else:
        raise Exception


def get_fit_data(data, n, window_len=11, m=30, nopeaks=4, spec=False):
    #new get fit data.
    #will obtain first three maxima by the old method and the extrapolate
    #to get the other peaks
    if n < 2:
        raise ValueError('n must be greater than 1')
    if n < nopeaks:
        nopeaks = n

    smoothed_data = smooth.smooth(data, window_len)

    params = list(getFitData(data, nopeaks, m=m))
    peak_space = mean([params[(i+1)*3 + 1] - params[i*3 + 1]
                      for i in range(1, nopeaks - 1)])

    last_peak = [params[-2]]
    for i in range(nopeaks, n):
        new_peak = last_peak[0] + peak_space
        params.append(smoothed_data[new_peak])
        params.append(new_peak)
        params.append(40)
        last_peak[0] = new_peak

    if spec:
        params.append(0.2)
        params.append(1)

    return array(params)


def odr_fit(n, data, func, params=None, fixed_params=None, nopeaks=4, **kw):
    """#n is the number of gaussians to fit
    params is an array of parameters for each gaussian (A, mu, sigma)
    for fixed params:
    'A value of 0 fixes the parameter, a value > 0 makes the parameter free.'
    """

    if params is None:
        params = get_fit_data(data, n, nopeaks=nopeaks)

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


def odr_spec(n, data, params, prop_func='erlang', fixed_params=None):

    def func(parameters, x):
        return spec_func(x, n, prop_func, *parameters)

    return odr_fit(n, data, func, params, fixed_params)

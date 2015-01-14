import numpy as np
from scipy.special import gamma
from math import sqrt, pi, exp, factorial


def amp_gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


def frac_gauss(x, mu, sigma, t=None):
    return np.exp(-0.5*pow((x-mu)/sigma, 2)) / (sigma * sqrt(2*pi))


def tail_gauss(x, mu, sigma, t):
    x = x.astype('float')
    arr = np.piecewise(x, [x >= (mu+t)], [tail, frac_gauss], *(mu, sigma, t))
    return arr


def tail(x, mu, sigma, t):
    return np.exp(0.5*(-t*(2*x-2*mu-t))/sigma**2) / (sigma * sqrt(2*pi))


def gen_gauss(x, *p):
    mu, alpha, beta = p
    return beta/(2*alpha*gamma(1./beta)) * np.exp(-pow(abs(x-mu)/alpha, beta))


def multigauss(x, n, *p):
    A = [p[i*3] for i in range(n)]
    mu = [p[i*3+1] for i in range(n)]
    sigma = [p[i*3+2] for i in range(n)]

    gaussArray = [amp_gauss(x, *[A[i], mu[i], sigma[i]])for i in range(n)]
    return np.array(gaussArray).sum(axis=0)


def mod_erlang(n, p, nu):
    """modified erlang distribution used by FACT to fit sipm dark spectra
    as seen in: http://arxiv.org/abs/1403.5747"""

    q = p * exp(-p)
    num = pow(n*q, n-1)
    denum = pow(factorial(n-1), nu)
    return num / denum


def poisson(n, p, nu=0):
    return pow(p, n-1) / factorial(n)


def line(x, *p):
    a, b = p
    return a*x + b


funcs = {'erlang': mod_erlang, 'poisson': poisson}


def spec_func(x, n, func='erlang', *params):
    try:
        function = funcs[func]
    except KeyError:
        print func
        raise KeyError

    mu = [params[i*2] for i in range(n)]
    sigma = [params[i*2+1] for i in range(n)]
    p, nu, A = params[-3:]

    return A*sigma[0]*sqrt(2*pi) * np.array(
        [function(i+1, p, nu) * frac_gauss(x, mu[i], sigma[i])
            for i in range(n)]).sum(axis=0)

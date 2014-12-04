import numpy as np
#from cfit_funcs import factorial
from math import sqrt, pi, exp, factorial


def gauss(x, *p):
    A, mu, sigma = p

    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


def multigauss(x, n, *p):
    A = [p[i*3] for i in range(n)]
    mu = [p[i*3+1] for i in range(n)]
    sigma = [p[i*3+2] for i in range(n)]

    gaussArray = [gauss(x, *[A[i], mu[i], sigma[i]])for i in range(n)]
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


funcs = {'erlang': mod_erlang, 'poisson': poisson}


def spec_func(x, n, func='erlang', *params):
    if not func in funcs:
        print func
        raise KeyError
    else:
        function = funcs[func]

    A = [params[i*3] for i in range(n)]
    mu = [params[i*3+1] for i in range(n)]
    sigma = [params[i*3+2] for i in range(n)]
    p, nu = params[-2:]

    a = [sig * sqrt(2*pi) for sig in sigma]

    return A[0] * a[0] * np.array([function(i+1, p, nu) *
                                  np.exp(-0.5*pow((x - mu[i]) / sigma[i], 2))
                                  / a[i] for i in range(n)]).sum(axis=0)

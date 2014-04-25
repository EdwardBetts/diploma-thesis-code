import numpy as np
from scipy.optimize import curve_fit

def gauss(x,*p):
    A,mu,sigma = p

    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def multigauss(x, n, *p):
    A = [p[i*3] for i in range(n)]
    mu = [p[i*3+1] for i in range(n)]
    sigma = [p[i*3+2] for i in range(n)]

    gaussArray = [gauss(x, *[A[i],mu[i],sigma[i]])for i in range(n)]
    return np.array(gaussArray).sum(axis=0)

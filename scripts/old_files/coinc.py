import numpy as np
import math

def poisson(N,nMean):
    return nMean**N * np.exp(- nMean)/math.factorial(N)

def coinc(rate,_delta_t, n = 2):
    delta_t = _delta_t * 10**(-9)
    coinc_rate = n * rate * delta_t**(n-1)
    print 'coincidence rate = %d' %coinc_rate

    return coinc_rate

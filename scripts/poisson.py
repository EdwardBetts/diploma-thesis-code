import numpy as np
import math
import multigauss

def poisson(N,nMean):
    return nMean**N * np.exp(- nMean)/math.factorial(N)

def getNMean(N_0, NAll):
    return - np.log(N_0 / NAll)

def getPoisson(n, fitData, xdata):
    x = [i for i in range(n)]
    #get histogram values from gaussians
    temp = [multigauss.gauss(xdata,*[fitData[0][i*3],fitData[0][i*3+1], fitData[0][i*3+2]]).sum() for i in range(n)]
    N_0 = temp[0]
    NAll = np.array(temp).sum()
    nMean = getNMean(N_0, NAll)

    y = [poisson(i, nMean) * NAll for i in range(n)]
    
    temp = np.array(temp)
    p1 = temp[2:].sum() / temp.sum()
    p0 = temp[1:].sum() / temp.sum()
    p1Poisson = 1 - (1 + nMean) * np.exp( - nMean)
  
    #absolute crosstalk prob
    CT = (p1 - p1Poisson) / p0
    return np.array(y),CT

def getCrossTalk(n,fitData,xdata):
    temp = [multigauss.gauss(xdata,*[fitData[0][i*3],fitData[0][i*3+1],fitData[0][i*3+2]]).sum() for i in range(n)]

    #get n mean poisson from 0 pe bin
    n_0 = temp[0]
    n_all = np.array(temp).sum()
    n_mean_poisson = getNMean(n_0,n_all)
    print n_mean_poisson

    #get n mean from histogram
    n_mean_spec = 0
    for i in range(n):
        n_mean_spec += i * temp[i] / n_all
    print n_mean_spec

    ct_prob = n_mean_spec / n_mean_poisson
    ct = (ct_prob - 1) * 100
 
    return ct  

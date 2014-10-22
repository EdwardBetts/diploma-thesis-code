import sys

def getGain(mu, k = 31.23,mode = 'highres'):
    e = 1.6021766 * 10**(-19)
    if mode == 'highres':
        M_QDC = 25 * 10**(-15)
    elif mode == 'lowres':
        M_QDC = 8 * 25 * 10**(-15)
    else:
        print 'wrong mode' 

    return (mu * M_QDC / e / k)

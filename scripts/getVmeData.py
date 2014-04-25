import sys
import numpy as np
import glob

def getData(filename, channelNo = 18):
    temp = np.loadtxt(filename)
    xData = [temp[i][0] for i in range(len(temp)) if i > 0]
    yData = [temp[i][channelNo] for i in range(len(temp)) if i > 0]

    return xData, yData

def get_data(filename, channelNo=18):
    data = np.loadtxt(filename, unpack=True)
    return data[1], data[channelNo+1]

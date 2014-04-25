#! /usr/bin/env python
import subprocess
import numpy as np
import time
from Fit import getFitData, get_first_in_range
from datetime import datetime
from smooth import smooth

def get_prop(n,cutoff):
    prop_data = []
    errors = []
    process_string = '/home/astro/jammer/Programme_vme/v965_1'
    time_unit = 5
    times = [time_unit*(1+i) for i in range(n)]
    #times = [15 for i in range(30)]

    for i in times:
        for j in range(3):
            subprocess.check_call([process_string,str(i)])
            temp = np.loadtxt('/home/astro/jammer/Programme_vme/data_v965.dat')
            data = np.array([temp[l][18] for l in range(len(temp)) if l > 0])
            sum_pedestal = data[:cutoff].sum()
            sum_all = data.sum()

            prop_data.append(sum_pedestal / sum_all)
            errors.append(np.sqrt(sum_pedestal) / sum_all)

    return times, prop_data, errors

#old
def data_to_nmean(prop_data, errors):
    nmean = -np.log(prop_data)
    err_nmean = -1. / np.array(prop_data) * errors

    return nmean, err_nmean

def nmean_n_times(time, filename, n, channels=(20,), cutoffs=False, pmt=None, dir_=None):
    process_string = '/home/astro/jammer/Programme_vme/v965_1'
    if dir_ is None:
        dir_ = '/home/astro/jammer/diplom/led/pde_test/'
    from datetime import datetime
    f = open(dir_ + filename,'w') #open file to overwrite anything
    f.close()

    i = 0
    if n == 0:
        n = np.infty
    while i<n:
        subprocess.check_call([process_string,str(time)])
        temp = np.loadtxt('/home/astro/jammer/Programme_vme/data_v965.dat', unpack=True)

        #list of tuples (n_mean, err_n_mean)
        prop_err = []
        #pmt-data
        pmt_data = None

        #get all sipm-channels
        for index, channel in enumerate(channels):
            data = temp[channel]
            if cutoffs is False:
                cutoff = get_cutoff(data)
                if cutoff is False:
                    cutoff = get_cutoff(data, guess=int(raw_input()))
            else:
                cutoff = get_cutoff(data, guess=cutoffs[index])


            prop_, error_ = get_prop_error(data, cutoff)
            prop_err.append((prop_, error_))

        #check pmt-channel,
        #pmt variable format should be: (channel-#, cutoff)
        if pmt is not None:
            channel, cutoff = pmt
            data = temp[channel]
            prop_, error_ = get_prop_error(data, cutoff)
            prop_err.append((prop_, error_))


        f = open(dir_ + filename,'a')
        f.write('#' + str(datetime.now().replace(microsecond = 0)) + '\n')

        #print date and time without without microseconds, format: h:m:s; '#' for np.loadtxt
        for dat in prop_err:
            f.write(' '.join(map(str,data_to_nmean(*dat))) + ' ')
        f.write('\n')
        f.close()
        i += 1

        if n != np.infty:
            print '%i of %i' %(i,n)
        else:
            print i
        for chan, j in enumerate(channels):
            print 'N_mean = %.3f' %(-np.log(prop_err[chan][0]))
        if pmt is not None:
            print 'PMT n_mean: %.3f' %-np.log(prop_err[-1][0])

#old
def get_prop_error(data,cutoff):
        sum_pedestal = data[:cutoff].sum()
        sum_all = data.sum()
        prop = sum_pedestal / sum_all
        error = np.sqrt(sum_pedestal) / sum_all

        return prop, error

def get_threshold(data, w_len, m=30):
    fit_data = getFitData(data, 3, m=m, window_len=w_len)

    #11 is the window_len in smooth, this shifts the whole histogram to +5.5 indices
    return (fit_data[4] + fit_data[1]) / 2 - 10

def get_cutoff(data, w_len=50, m=30, guess=None):
    if guess is None:
        approx = get_threshold(data, m)
    else:
        approx = guess

    smoothData = smooth(data)
    return get_first_in_range(smoothData[approx - w_len:approx + w_len],
                              approx, w_len)

def to_file(filename, threshold, data):
    data = np.array(data[1])
    prop = data[:threshold].sum() / data.sum()
    nmean = - np.log(prop)
    append_file(filename, nmean)
    print nmean


def append_file(filename, nmean):
    f = open('/home/astro/jammer/diplom/messdaten/pde_tests/pde_first/' + filename,'a')
    np.savetxt(f,[nmean])
    f.close()

def extract_file(filename):
    data = np.loadtxt(filename, unpack=True)

    with open(filename, 'r') as f:
        temp = f.readlines()

    times = []
    for line in temp:
        if line[0] == '#':
            times.append(line[1:20])

    dates = [datetime.strptime(i, '%Y-%m-%d %H:%M:%S') for i in times]

    return data, dates


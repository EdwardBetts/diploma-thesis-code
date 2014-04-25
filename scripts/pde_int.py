from glob import glob
from read_drs import event_generator, return_dtype
from numpy import histogram, sum, log, unique, argmin, mean, inf, std
from numpy import loadtxt, where, sqrt, arctan, sin, fromstring, r_
from smooth import smooth


def extract_pde(dark_counts=False, cutoff=32675, w_start=340,
                w_len=180, small_pmt_win=False, print_pe=False):

    filelist = [fname for fname in glob('*_*') if not '.' in fname]
    wavelengths = unique(i[:-2] for i in filelist)
    QE_file = '/home/jammer/diplom/calib/pmt2c.dat'
    f = loadtxt(QE_file, unpack=True)

    pde_array = []
    wavelength_array = []
    sterr_array = []
    syserr_array = []

    my_dtype = return_dtype(2)

    if not dark_counts:
        dark_counts = get_dark_counts(filelist, w_start, w_len)

    for wavelength in wavelengths:
        A = []
        B = []
        quantum_eff = f[2][where(f[0] == float(wavelength))[0][0]]

        for i in range(2):
            #data = read_drs_binary(wavelength + '_' +str(i+1), 2)
            fname = wavelength + '_' + str(i+1)
            #create sipm spectrum
            with open(fname, 'r') as fl:
                sipm_data = (fromstring(event, my_dtype)[0][5]
                             for event in event_generator(fl, 2))

                sipm_max = [sum(event[w_start:w_start + w_len])
                            for event in sipm_data
                            if sum(event[w_start-w_len:w_start]) > cutoff]
                sipm_spec = histogram(sipm_max, bins=2048)[0]

            #create pmt spectrum
            with open(fname, 'r') as fl:
                pmt_data = (fromstring(event, my_dtype)[0][7]
                            for event in event_generator(fl, 2))

                pmt_sum = [sum(event[500:800]) for event in pmt_data]
                pmt_spec = histogram(pmt_sum, bins=2048,
                                     range=(5*10**6, 10**7))[0]

            sipm_guess = get_guess(sipm_spec)
            sipm_pe = get_nmean_errors(sipm_spec,
                                       get_cutoff(sipm_spec, guess=sipm_guess))
            pmt_pe = get_nmean_errors(pmt_spec, get_cutoff(pmt_spec,
                                      1950, window=30))
            if pmt_pe == inf or small_pmt_win:
                pmt_pe = get_nmean_errors(pmt_spec, get_cutoff(pmt_spec,
                                          1950, window=10))
            if print_pe:
                print 'sipm pe' + str(sipm_pe)
                print 'pmt pe' + str(pmt_pe)

            A.append(sipm_pe)
            B.append(pmt_pe)

        real_wl = correct_wavelength(float(wavelength))
        pde, sterr, syserr = calc_pde_with_errors(A[0], A[1], B[0], B[1],
                                                  quantum_eff,
                                                  dkcts=dark_counts)
        print ' '.join((str((pde, sterr, syserr)), wavelength))
        pde_array.append(pde)
        wavelength_array.append(real_wl)
        sterr_array.append(sterr)
        syserr_array.append(syserr)

    return pde_array, wavelength_array, sterr_array, syserr_array


def get_n_mean(data, cutoff):
    pedestal = float(sum(data[cutoff:]))
    all_ = sum(data)
    return - log(pedestal / all_)


def get_cutoff(data, guess, window=50):
    smoothdata = smooth(data)
    amin = argmin(smoothdata[guess - window: guess + window])
    return amin - 5.5 + guess - window


def calculate_pde(A_1, A_2, B_1, B_2, pde_pmt, dark_count=0.16145021896652537):
    factor = (A_1 - dark_count) * (A_2 - dark_count) / B_1 / B_2
    return sqrt(factor) * pde_pmt


def correct_wavelength(wavelength):
    angle = arctan(.75 / 3)
    #given by the manufacturor of the ir-filters
    k = 0.11
    factor = 1 - k * sin(angle) * sin(angle)

    return wavelength * factor


def calc_pde_with_errors(A1, A2, B1, B2, qe_pmt, dkcts=0.16145021896652537):
    try:
        factor = (A1[0] - dkcts[0]) * (A2[0] - dkcts[0]) / B1[0] / B2[0]
    except TypeError:
        factor = (A1[0] - dkcts) * (A2[0] - dkcts) / B1[0] / B2[0]
        dkcts = (dkcts, dkcts / 2)

    pde = sqrt(factor) * qe_pmt
    #statistical error:
    sterr = sqrt((0.5 * factor**(0.5) / A1[0] * A1[1])**2 +
                 (0.5 * factor**(0.5) / A2[0] * A2[1])**2 +
                 (- 0.5 * factor**(0.5) / B1[0] * B1[1])**2 +
                 (- 0.5 * factor**(0.5) / B2[0] * B2[1])**2 +
                 (0.5 * factor / B1[0] / B2[0] *
                  (2 * dkcts[0] - A1[0] - A2[0]) * dkcts[1])**2) * qe_pmt
    #systematic error:
    syserr = sqrt((0.5 * factor**(0.5) / A1[0] * A1[2])**2 +
                 (0.5 * factor**(0.5) / A2[0] * A2[2])**2 +
                 (- 0.5 * factor**(0.5) / B1[0] * B1[2])**2 +
                 (- 0.5 * factor**(0.5) / B2[0] * B2[2])**2) * qe_pmt
    return pde, sterr, syserr


def get_nmean_errors(data, cutoff, swin=7):
    pedestal = float(sum(data[cutoff:]))
    all_ = sum(data)
    n_mean = -log(pedestal / all_)
    #statistical error
    sterr = sqrt(pedestal) / pedestal
    #systematic error:
    syserr = max([abs(n_mean - get_n_mean(data, cutoff - swin + 2 * swin * i))
                 for i in range(2)])
    return n_mean, sterr, syserr


def get_guess(spec, window_len=11, n=2, m=20):
    sdata = smooth(spec[::-1], window_len)

    #find values which are higher than neighbors
    test = r_[True, sdata[1:] > sdata[:-1]] & r_[sdata[:-1] > sdata[1:], True]

    #find values which are higher than the next m neighbors
    newInds = where(test)[0]
    values = [i - 5 for i in newInds[newInds > m] if sdata[i] ==
              max(sdata[i-m:i+m + 1]) and sdata[i] > 1.0]
    mu_values = [2048 - i for i in values[n - 2:n]]
    return mean(mu_values)


def get_dark_counts(filelist, w_start, w_len):
    my_dtype = return_dtype(2)
    darks = []
    for fn in filelist:
        with open(fn, 'r') as fl:
            sipm_data = (fromstring(event, my_dtype)[0][5]
                         for event in event_generator(fl, 2))
            sipm_max = [sum(event[w_start - w_len:w_start])
                        for event in sipm_data]
            sipm_spec = histogram(sipm_max, bins=2048)[0]
            guess = get_guess(sipm_spec)
            cutoff = get_cutoff(sipm_spec, get_cutoff(sipm_spec, guess))
            mn = get_n_mean(sipm_spec, cutoff)
            if mn < 1:
                darks.append(mn)
        """dkcts = median(darks)
        dkcts_err = max(abs(median(darks) - percentile(darks, i))
                        for i in (80, 20))"""
    dkcts = mean(darks)
    dkcts_err = std(darks)

    print dkcts, dkcts_err
    return dkcts, dkcts_err

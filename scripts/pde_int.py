from glob import glob
from read_drs import event_generator, return_dtype
from numpy import histogram, sum, unique, mean, inf, std
from numpy import loadtxt, where, fromstring
from pdetools import calc_pde_with_errors, correct_wavelength
from spectools import get_guess, get_cutoff, get_nmean_errors, get_n_mean
from spectools import get_co_exp


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
                            for event in sipm_data]
                sipm_spec = histogram(sipm_max, bins=2048)[0]

            #create pmt spectrum
            with open(fname, 'r') as fl:
                pmt_data = (fromstring(event, my_dtype)[0][7]
                            for event in event_generator(fl, 2))

                pmt_sum = [sum(event[600:900]) for event in pmt_data]
                pmt_spec = histogram(pmt_sum, bins=2048)[0]

            sipm_guess = get_guess(sipm_spec)
            sipm_pe = get_nmean_errors(sipm_spec,
                                       get_cutoff(sipm_spec, guess=sipm_guess))
            pmt_pe = get_nmean_errors(pmt_spec, get_co_exp(pmt_spec, 1950))
            #if pmt_pe == inf or small_pmt_win:
            #    pmt_pe = get_nmean_errors(pmt_spec, get_cutoff(pmt_spec,
            #                              1950, window=10))
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

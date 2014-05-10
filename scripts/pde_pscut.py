from glob import glob
from numpy import unique, loadtxt, fromstring, histogram, inf, where
from read_drs import scatter, dark_scatter, return_dtype, event_generator
from pscut import get_line, get_hist
from spectools import get_guess, get_cutoff, get_nmean_errors
from pdetools import calc_pde_with_errors, correct_wavelength


def extract_pde(dkcts, int_start, print_pe=False, linepars=None,
                qe_loc='/home/jammer/diplom/calib/pmt2c.dat'):
    filelist = [f for f in glob('*_*') if not '.' in f]
    wavelengths = unique(i[:-2] for i in filelist)
    QE_file = qe_loc
    f = loadtxt(QE_file, unpack=True)

    pde_array = []
    wavelength_array = []
    sterr_array = []
    syserr_array = []

    for wavelength in wavelengths:
        A = []
        B = []
        quantum_eff = f[2][where(f[0] == float(wavelength))[0][0]]

        for i in range(2):
            fname = wavelength + '_' + str(i+1)
            intd, mind = scatter(fname, int_start)
            if linepars is None:
                a, b = get_line(intd, mind)
            else:
                a, b = linepars
            sipm_hist = get_hist(intd, mind, (a, b))[0]

            with open(fname, 'r') as fl:
                my_dtype = return_dtype(2)
                pmt_data = (fromstring(event, my_dtype)[0][7]
                            for event in event_generator(fl, 2))

                pmt_sum = [sum(event[500:800]) for event in pmt_data]
                pmt_spec = histogram(pmt_sum, bins=2048,
                                     range=(5*10**6, 10**7))[0]

            sipm_guess = get_guess(sipm_hist)
            sipm_pe = get_nmean_errors(sipm_hist,
                                       get_cutoff(sipm_hist, guess=sipm_guess))
            pmt_pe = get_nmean_errors(pmt_spec, get_cutoff(pmt_spec,
                                      1950, window=30))
            if pmt_pe == inf:
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
                                                  dkcts=dkcts)
        print ' '.join((str((pde, sterr, syserr)), wavelength))
        pde_array.append(pde)
        wavelength_array.append(real_wl)
        sterr_array.append(sterr)
        syserr_array.append(syserr)

    return pde_array, wavelength_array, sterr_array, syserr_array


def get_dark_counts(dark_st, int_st, filelist, protoevent, linepars=None):
    darks = []
    for f in filelist:
        if linepars is None:
            intd, mind = scatter(f, int_st)
            a, b = get_line(intd, mind)
        else:
            a, b = linepars
        intd, mind = dark_scatter(f, dark_st, protoevent)
        hist = get_hist(intd, mind, (a, b))[0]

        guess = get_guess(hist)
        sipm_pe = get_nmean_errors(hist, get_cutoff(hist, guess=guess))
        darks.append(sipm_pe)
        print sipm_pe
    return darks

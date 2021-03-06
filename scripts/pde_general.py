from glob import glob
from numpy import unique, loadtxt, fromstring, histogram, inf, where
from read_drs import return_dtype, trace_gen
from spectools import get_guess, get_cutoff, get_nmean_errors, get_co_exp
from pdetools import calc_pde_with_errors, correct_wavelength


def extract_pde(func, dkcts, filename, wavelength, print_pe=False,
                qe_loc='/home/astro/jammer/diplom/pde_measurement/QPMT3.dat',
                pmtwin=(610, 670), **kwargs):

    wls, currents, qes = loadtxt(qe_loc, unpack=True)
    qe = qes[where(wls == float(wavelength))[0][0]]

    # factors for calculating pde
    A = []
    B = []

    for i in range(2):
        fname = ''.join((filename, '_', str(i+1), '.dat'))
        sipm_hist = func(fname, **kwargs)[0]

        with trace_gen(fname, 2, 'c2') as gen:
            pmt_sum = [sum(event[pmtwin[0]:pmtwin[1]]) for event in gen]
            pmt_spec = histogram(pmt_sum, bins=2048)[0]

        sipm_guess = get_guess(sipm_hist)
        sipm_pe = get_nmean_errors(sipm_hist,
                                   get_cutoff(sipm_hist, guess=sipm_guess))
        pmt_pe = get_nmean_errors(pmt_spec, get_co_exp(pmt_spec, 1900))
        if pmt_pe == inf:
            pmt_pe = get_nmean_errors(pmt_spec, get_cutoff(pmt_spec,
                                      1900, window=10))
        if print_pe:
            print 'sipm pe' + str(sipm_pe)
            print 'pmt pe' + str(pmt_pe)

        A.append(sipm_pe)
        B.append(pmt_pe)

    wl = correct_wavelength(float(wavelength))
    pde, sterr, syserr = calc_pde_with_errors(A[0], A[1], B[0], B[1],
                                              qe, dkcts=dkcts)
    print ' '.join((str((pde, sterr, syserr)), str(wavelength)))

    return pde, wl, sterr, syserr


# old
def process_folder(func, dkcts, print_pe=False,
                   qe_loc='/home/jammer/diplom/pde_measurement/QPMT3.dat',
                   **kwargs):
    """pde calculation with arbitrary function to get a sipm hist as argument
    'func'. this function must take a filename as the first argument,
    other arguments to this function will be give with **kwargs"""

    filelist = [f for f in glob('*_*') if '.' not in f]
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
            sipm_hist = func(fname, **kwargs)[0]

            with open(fname, 'r') as fl:
                my_dtype = return_dtype(2)
                pmt_data = (fromstring(event, my_dtype)[0][7]
                            for event in event_generator(fl, 2))

                pmt_sum = [sum(event[600:900]) for event in pmt_data]
                pmt_spec = histogram(pmt_sum, bins=2048)[0]

            sipm_guess = get_guess(sipm_hist)
            sipm_pe = get_nmean_errors(sipm_hist,
                                       get_cutoff(sipm_hist, guess=sipm_guess))
            pmt_pe = get_nmean_errors(pmt_spec, get_co_exp(pmt_spec, 1950))
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

"""
from read_drs import scatter, dark_scatter
from pscut import get_line, get_hist
def dark_counts_old_pscut(dark_st, int_st, filelist, protoevent, linepars=None,
                          **kwargs):
    darks = []
    for f in filelist:
        if linepars is None:
            intd, mind = scatter(f, int_st)
            a, b = get_line(intd, mind)
        else:
            a, b = linepars
        intd, mind = dark_scatter(f, dark_st, protoevent)
        hist = get_hist(intd, mind, (a, b), **kwargs)[0]

        guess = get_guess(hist)
        sipm_pe = get_nmean_errors(hist, get_cutoff(hist, guess=guess))
        darks.append(sipm_pe)
        print sipm_pe
    return darks


def hist_old_pscut(fname, int_start, linepars=None):
    intd, mind = scatter(fname, int_start)
    if linepars is None:
        a, b = get_line(intd, mind)
    else:
        a, b = linepars
    sipm_hist = get_hist(intd, mind, (a, b))[0]
    return sipm_hist"""

from glob import glob
from read_drs import read_drs_binary, event_generator, return_dtype
from scipy.signal import iirfilter, filtfilt
from numpy import histogram, sum, log, unique, argmin
from numpy import loadtxt, where, sqrt, mean, arctan, sin, fromstring
from smooth import smooth

def extract_pde(sipm_guess=1925, sipm_cutoff=32675):

    filelist = [fname for fname in glob('*') if not '.' in fname]
    wavelengths = unique(i[:-2] for i in filelist)
    QE_file = '/home/astro/jammer/diplom/pde_measurement/QPMT2.dat'
    b, a = iirfilter(1, 0.07, btype='lowpass')
    f = loadtxt(QE_file, unpack=True)

    pde_array = []
    wavelength_array = []
    err_array = []

    my_dtype = return_dtype(2)

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

                sipm_max = [min(fevent[340:510])
                            for fevent in (filtfilt(b, a, event)
                            for event in sipm_data) if min(fevent[:340]) > sipm_cutoff]
                sipm_spec = histogram(sipm_max, bins=2048)[0]

            #create pmt spectrum
            with open(fname, 'r') as fl:
                pmt_data = (fromstring(event, my_dtype)[0][7]
                            for event in event_generator(fl, 2))

                pmt_sum = [sum(event[550:800]) for event in pmt_data]
                pmt_spec = histogram(pmt_sum, bins=2048)[0]

            sipm_pe = get_nmean_errors(sipm_spec, get_cutoff(sipm_spec, guess=sipm_guess))
            pmt_pe = get_nmean_errors(pmt_spec, get_cutoff(pmt_spec, guess=1950))

            A.append(sipm_pe)
            B.append(pmt_pe)

        real_wl = correct_wavelength(float(wavelength))
        pde, err = calc_pde_with_errors(A[0], A[1], B[0], B[1], quantum_eff)
        print ' '.join((str(pde), wavelength)) + ' calculated with a QE of ' + str(quantum_eff)
        pde_array.append(pde)
        wavelength_array.append(real_wl)
        err_array.append(err)

    return pde_array, wavelength_array, err_array


def get_n_mean(data, cutoff):
    pedestal = float(sum(data[cutoff:]))
    all_ = sum(data)
    return - log(pedestal / all_)

def get_cutoff(data, guess, window=50):
    smoothdata = smooth(data)
    amin = argmin(smoothdata[guess - window: guess + window])
    return amin - 5.5 + guess - 50

def calculate_pde(A_1, A_2, B_1, B_2, pde_pmt, dark_count=0.15676227395961931):
    factor = (A_1 - dark_count) * (A_2 - dark_count)/ B_1 / B_2
    return sqrt(factor) * pde_pmt

def correct_wavelength(wavelength):
    angle = arctan(.75 / 3)
    #given by the manufacturor of the ir-filters
    k = 0.11
    factor = 1 - k * sin(angle) * sin(angle)

    return wavelength * factor

def calc_pde_with_errors(A1, A2, B1, B2, qe_pmt, dkcts=0.15676227395961931):
    factor = (A1[0] - dkcts) * (A2[0] - dkcts)/ B1[0] / B2[0]
    pde = sqrt(factor) * qe_pmt
    err = (sqrt(0.5 * A2[0] / B1[0] / B2[0] * A1[1]) +
           sqrt(0.5 * A1[0] / B1[0] / B2[0] * A2[1]) +
           sqrt(-0.5 * A1[0] * A2[0] / B1[0]**2 / B2[0] * B1[1]) +
           sqrt(-0.5 * A1[0] * A2[0] / B1[0] / B2[0]**2 * B2[1]))* factor**(-1)
    return pde, err

def get_nmean_errors(data, cutoff):
    pedestal = float(sum(data[cutoff:]))
    all_ = sum(data)
    err = sqrt(pedestal) / pedestal
    return -log(pedestal / all_), err

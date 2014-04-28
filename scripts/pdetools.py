from numpy import sqrt, arctan, sin


def calculate_pde(A_1, A_2, B_1, B_2, pde_pmt, dkcts):
    factor = (A_1 - dkcts) * (A_2 - dkcts) / B_1 / B_2
    return sqrt(factor) * pde_pmt


def correct_wavelength(wavelength):
    angle = arctan(.75 / 3)
    #given by the manufacturor of the ir-filters
    k = 0.11
    factor = 1 - k * sin(angle) * sin(angle)

    return wavelength * factor


def calc_pde_with_errors(A1, A2, B1, B2, qe_pmt, dkcts):
    a = 0.006
    try:
        factor = (A1[0]-dkcts[0]) * (A2[0]-dkcts[0]) / (B1[0]-a) / (B2[0]-a)
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

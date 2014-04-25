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
    factor = (A1[0] - dkcts) * (A2[0] - dkcts) / B1[0] / B2[0]
    pde = sqrt(factor) * qe_pmt
    err = (sqrt(0.5 * A2[0] / B1[0] / B2[0] * A1[1]) +
           sqrt(0.5 * A1[0] / B1[0] / B2[0] * A2[1]) +
           sqrt(-0.5 * A1[0] * A2[0] / B1[0]**2 / B2[0] * B1[1]) +
           sqrt(-0.5 * A1[0] * A2[0] / B1[0] / B2[0]**2 * B2[1]))*factor**(-1)
    return pde, err

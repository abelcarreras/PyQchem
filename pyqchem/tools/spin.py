import numpy as np


def spin_matrices(s_bra, s_ket):
    """
    Get spin matrices Sx, Sy, Sz between two spin states (s,m) and (s',m') such that
    sx = < m' | Sx | m >, sy = < m' | Sy | m > and sz = < m' | Sz | m >

    :param s_bra: total spin (s) of (s',m') state
    :param s_ket: total spin (s) of (s,m) state
    :return: Sx, Sy, Sz
    """

    def get_multiplicity(s):
        return int(2 * s + 1)

    def are_equal(a, b,  thresh=1e-4):
        return abs(a - b) < thresh

    # spin-multiplicities
    mul_bra = get_multiplicity(s_bra)
    mul_ket = get_multiplicity(s_ket)

    # initialize Sx, Sy, Sz
    Sx = np.zeros((mul_bra, mul_ket), dtype=complex)
    Sy = np.zeros((mul_bra, mul_ket), dtype=complex)
    Sz = np.zeros((mul_bra, mul_ket), dtype=complex)

    # build spin matrices
    for i, sz_bra in enumerate(np.arange(-s_bra, s_bra + 1)):
        for j, sz_ket in enumerate(np.arange(-s_ket, s_ket + 1)):

            if are_equal(sz_bra, sz_ket):
                Sz[i, j] = sz_ket

            elif are_equal(sz_bra,  sz_ket + 1):
                s = s_bra
                Sx[i, j] = 0.5 * np.sqrt(s * (s + 1) - sz_bra * sz_ket)
                Sy[i, j] = -0.5j * np.sqrt(s * (s + 1) - sz_bra * sz_ket)

            elif are_equal(sz_bra, sz_ket - 1):
                s = s_bra
                Sx[i, j] = 0.5 * np.sqrt(s * (s + 1) - sz_bra * sz_ket)
                Sy[i, j] = 0.5j * np.sqrt(s * (s + 1) - sz_bra * sz_ket)

    return Sx, Sy, Sz


def s2_to_s(s2):
    """
    get total spin (s) from s2

    :param s2: s2
    :return: total spin (s)
    """
    return 0.5*(-1 + np.sqrt(1 + 4 * s2))


def s_to_s2(s):
    """
    get s2 from total spin (s)

    :param s: total spin (s)
    :return: s2
    """
    return s * (s + 1)


def get_s_sz_from_configuration(n_alpha, n_beta, s2):
    """
    Get s and sz from electronic configuration

    :param n_alpha: number of alpha electrons
    :param n_beta:  number of beta electrons
    :param s2: s2 of state
    :return: s, sz
    """
    n_tot = n_alpha + n_beta
    n_dif = n_alpha - n_beta

    s = s2_to_s(s2)

    sz = 0.5 * n_dif
    if np.mod(n_tot, 2) == 0:
        s = int(s)
        sz = int(sz)
    else:
        s = int(s + 0.5) - 0.5
        sz = int(sz + 0.5) - 0.5

    return s, sz


if __name__ == '__main__':

    # example
    n_alpha = 3
    n_beta = 2
    s2 = 3.7500
    # s2 = s_to_s2(1.5)

    s_bra, Sz_bra = get_s_sz_from_configuration(n_alpha, n_beta, s2)
    s_ket, Sz_ket = get_s_sz_from_configuration(n_alpha, n_beta, s2)

    print("KET: S ,Sz = {}, {}".format(s_ket, Sz_ket))
    print("BRA: S ,Sz = {}, {}".format(s_bra, Sz_bra))

    Sx, Sy, Sz = spin_matrices(s_bra, s_ket)

    print('Spin matrices')
    print(Sx.real)
    print(Sy.imag)
    print(Sz.real)

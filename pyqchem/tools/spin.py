import numpy as np


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


def spin_matrices(s):
    """
    Get spin matrices Sx, Sy, Sz between two spin states (s,m) and (s,m') such that
    sx = < m' | Sx | m >, sy = < m' | Sy | m > and sz = < m' | Sz | m >

    :param s: total spin (s)
    :return: Sx, Sy, Sz
    """

    def are_equal(a, b,  thresh=1e-4):
        return abs(a - b) < thresh

    def sz_values(s):
        return np.arange(-s, s + 1)

    # spin-multiplicities
    multiplicity = len(sz_values(s))

    # initialize Sx, Sy, Sz
    Sx = np.zeros((multiplicity, multiplicity), dtype=complex)
    Sy = np.zeros((multiplicity, multiplicity), dtype=complex)
    Sz = np.zeros((multiplicity, multiplicity), dtype=complex)

    # build spin matrices
    for i, sz_bra in enumerate(sz_values(s)):
        for j, sz_ket in enumerate(sz_values(s)):

            if are_equal(sz_bra, sz_ket):
                Sz[i, j] = sz_ket

            if are_equal(sz_bra,  sz_ket + 1):
                Sx[i, j] = 0.5 * np.sqrt(s_to_s2(s) - sz_bra * sz_ket)
                Sy[i, j] = -0.5j * np.sqrt(s_to_s2(s) - sz_bra * sz_ket)

            if are_equal(sz_bra, sz_ket - 1):
                Sx[i, j] = 0.5 * np.sqrt(s_to_s2(s) - sz_bra * sz_ket)
                Sy[i, j] = 0.5j * np.sqrt(s_to_s2(s) - sz_bra * sz_ket)

    return Sx, Sy, Sz


# Alternative implementation (based on https://easyspin.org/easyspin/documentation/spinoperators.html)
def spin_matrices_2(s):
    """
    Get spin matrices Sx, Sy, Sz between two spin states (s,m) and (s,m') such that
    sx = < m' | Sx | m >, sy = < m' | Sy | m > and sz = < m' | Sz | m >

    :param s: total spin (s)
    :return: Sx, Sy, Sz
    """

    def sz_values(s):
        return np.arange(-s, s + 1)

    # spin-multiplicities
    multiplicity = len(sz_values(s))

    delta = np.identity(multiplicity)
    delta_plus = np.diag(np.ones(multiplicity), k=1)[:-1, :-1]
    delta_minus = np.diag(np.ones(multiplicity), k=-1)[:-1, :-1]

    factor = np.sqrt(s_to_s2(s) - np.outer(sz_values(s), sz_values(s)))

    Sz = delta * sz_values(s)
    Sx = 0.5 * (delta_plus + delta_minus) * factor
    Sy = 0.5j * (delta_plus - delta_minus) * factor

    Sp = delta_plus * factor  # S_+
    Sm = delta_minus * factor  # S_-
    S2 = delta * s_to_s2(s)

    return Sx, Sy, Sz


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

    s, sz = get_s_sz_from_configuration(n_alpha, n_beta, s2)
    print("S ,Sz = {}, {}".format(s, sz))

    Sx, Sy, Sz = spin_matrices(s)
    # Sx, Sy, Sz = spin_matrices_2(s)

    print('Spin matrices')
    print(Sx.real)
    print(Sy.imag)
    print(Sz.real)

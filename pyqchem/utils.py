import re
import numpy as np
from scipy.optimize import leastsq


def standardize_vector(vector):
    import numpy as np
    if vector[0] != 0:
        if vector[0] < 0:
            vector = np.array(vector) * -1
            vector = vector.tolist()
    elif vector[1] != 0:
        if vector[1] < 0:
            vector = np.array(vector) * -1
            vector = vector.tolist()
    else:
        if vector[2] < 0:
            vector = np.array(vector) * -1
            vector = vector.tolist()

    for i in range(3):
        vector[i] = vector[i] + 0

    return vector


def search_bars(output, from_position=0, bar_type='---'):
    output = output[from_position:]
    positions = []
    previous = 0
    for m in re.finditer(bar_type, output):
        if m.start() > previous + 1:
            positions.append(m.start() + from_position)
        previous = m.end()

    return positions


def is_transition(configuration, reference, n_electron=1, max_jump=10):
    """
    Determine if a configuration correspond to a transition of n_electron
    :param configuration: dictionary containing the configuration to be analyzed
    :param reference: reference configuration (in general lowest energy Slater determinant)
    :param n_electron:
    :param max_jump: Restrict to transitions with jumps less or equal to max_jump orbitals
    :return: True if conditions are met, otherwise False
    """

    alpha_diff = [int(i) - int(j) for i, j in zip(configuration['alpha'], reference['alpha'])]
    beta_diff = [int(i) - int(j) for i, j in zip(configuration['beta'], reference['beta'])]

    ini_alpha = np.where(np.array(alpha_diff) < 0)[0]
    fin_alpha = np.where(np.array(alpha_diff) > 0)[0]

    ini_beta = np.where(np.array(beta_diff) < 0)[0]
    fin_beta = np.where(np.array(beta_diff) > 0)[0]

    try:
        jump_alpha = np.max(fin_alpha) - np.min(ini_alpha)
    except ValueError:
        jump_alpha = 0

    try:
        jump_beta = np.max(fin_beta) - np.min(ini_beta)
    except ValueError:
        jump_beta = 0

    n_alpha = len(fin_alpha)
    n_beta = len(fin_beta)

    elec_condition = n_alpha + n_beta == n_electron
    jump_condition = jump_alpha <= max_jump and jump_beta <= max_jump

    return elec_condition and jump_condition


def get_ratio_of_condition(state, n_electron=1, max_jump=10):
    # reference = {'hole': '', 'alpha': '111000', 'beta': '111000', 'part': '', 'amplitude': 0.9777044}

    alpha_int = [int(i) for i in state['configurations'][0]['alpha']]
    ref_alpha = ''
    for i in range(len(alpha_int)):
        if i < np.sum(alpha_int):
            ref_alpha += '1'
        else:
            ref_alpha += '0'

    reference = {'hole': '', 'alpha': ref_alpha, 'beta': ref_alpha, 'part': ''}

    p = 0
    for configuration in state['configurations']:
        if is_transition(configuration, reference, n_electron, max_jump):
            p += configuration['amplitude']**2

    return p


def set_zero_coefficients(basis, mo_coeff, range_atoms):
    """
    set 0.0 to coefficients of functions centered in atoms 'range_atoms'

    :param basis: full basis dictionary
    :param mo_coeff: full molecular orbitals coefficients to be modiffied
    :param range_atoms: list containing the atom numbers whose coefficients will be set to zero
    :return:
    """

    functions_to_atom = []
    nat = len(basis['atoms'])
    for i in range(0, nat):
        nf = 0
        for shell in basis['atoms'][i]['shells']:
            nf += shell['functions']
        functions_to_atom.append(nf)
    functions_to_atom = functions_to_atom

    # print(funtions_to_atom)
    mo_coeff_a = np.array(mo_coeff['alpha'])
    for i in range_atoms:
        ini = np.sum(functions_to_atom[:i], dtype=int)
        fin = np.sum(functions_to_atom[:i+1], dtype=int)
        # print('ini', ini, 'fin', fin)
        mo_coeff_a[:, ini: fin] *= 0.0
    mo_coeff_zero = {'alpha': mo_coeff_a.tolist()}

    if 'beta' in mo_coeff:
        mo_coeff_b = np.array(mo_coeff['beta'])
        for i in range_atoms:
            ini = np.sum(functions_to_atom[:i], dtype=int)
            fin = np.sum(functions_to_atom[:i + 1], dtype=int)
            # print('ini', ini, 'fin', fin)
            mo_coeff_b[:, ini: fin] *= 0.0
        mo_coeff_zero['beta'] = mo_coeff_b.tolist()

    return mo_coeff_zero


def get_plane(coords, direction=None):
    """Generate initial guess"""

    coords = np.array(coords)
    p0 = np.cross(coords[0] - coords[2], coords[0] - coords[-1]) / np.linalg.norm(np.cross(coords[0] - coords[2],
                                                                                           coords[0] - coords[-1]))

    # Fitting function to a plane
    def fitfunc(p, coords):
        average = np.average(coords, axis=0)
        return np.array([np.dot(p, average - c) for c in coords])

    # Error function (including force norm(normal) = 1)
    errfunc = lambda p, x: fitfunc(p, x)**2 + (np.linalg.norm(p) - 1.0)**2

    p1, flag = leastsq(errfunc, p0, args=(coords,))


    # Check final result
    point = np.average(coords, axis=0).tolist()
    normal = np.array(p1)/np.linalg.norm(p1)

    if direction is not None:
        vector = coords[direction[1]] - coords[direction[0]]
        # proj = vector - np.dot(vector, normal)*normal
        projected = np.cross(normal,  np.cross(vector, normal))
        projected /= np.linalg.norm(projected)

        normal = standardize_vector(normal)
        projected = standardize_vector(projected)

        return point, normal, projected

    normal = standardize_vector(normal)
    return point, normal


def reorder_coefficients(occupations, coefficients):
    """
    Reorder the coefficients according to occupations. Occupated orbitals will be grouped at the beginning
    non occupied will be attached at the end

    :param occupations: list on integers (0 or 1) or list of Boolean
    :param coefficients:  dictionary containing the molecular orbitals coefficients {'alpha': coeff, 'beta:' coeff}.
                          coeff should be a list of lists (Norb x NBas)
    :return:
    """

    alpha_coefficients = []
    non_occupied = []
    for occ, coeff in zip(occupations['alpha'], coefficients['alpha']):
        if occ:
            alpha_coefficients.append(coeff)
        else:
            non_occupied.append(coeff)

    # non occupated attached at the end
    alpha_coefficients += non_occupied

    if not 'beta' in coefficients:
        coefficients['beta'] = coefficients['alpha']

    beta_coefficients = []
    non_occupied = []
    for occ, coeff in zip(occupations['beta'], coefficients['beta']):
        if occ:
            beta_coefficients.append(coeff)
        else:
            non_occupied.append(coeff)

    # non occupated attached at the end
    beta_coefficients += non_occupied

    return {'alpha': alpha_coefficients, 'beta': beta_coefficients}


def get_occupied_electrons(configuration, structure):
    # not working for particle, hole enabled configurations!
    alpha_e = np.sum([int(c) for c in configuration['alpha']])
    beta_e = np.sum([int(c) for c in configuration['beta']])
    return (structure.number_of_electrons + structure.charge - (alpha_e + beta_e))//2


if __name__ == '__main__':
    state = {'configurations': [{'hole': '', 'alpha': '110100', 'beta': '111000', 'part': '', 'amplitude': 0.5}]}
    b = get_ratio_of_condition(state)
    print('test', b)
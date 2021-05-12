from scipy.optimize import leastsq
from copy import deepcopy
import numpy as np
import warnings


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


def is_rasci_transition(configuration, reference, n_electron=1, max_jump=10):
    """
    Determine if a configuration corresponds to a transition of n_electron

    :param configuration: dictionary containing the configuration to be analyzed
    :param reference: reference configuration (in general lowest energy Slater determinant)
    :param n_electron:
    :param max_jump: Restrict to transitions with jumps less or equal to max_jump orbitals

    :return: True if conditions are met, otherwise False
    """

    warnings.warn('This function will be deprecated, use "is_transition" instead', DeprecationWarning)


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


def get_ratio_of_condition_rasci(state, n_electron=1, max_jump=10):
    # reference = {'hole': '', 'alpha': '111000', 'beta': '111000', 'part': '', 'amplitude': 0.9777044}

    warnings.warn('This function will be deprecated, use "get_ratio_of_condition" instead', DeprecationWarning)
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
        if is_rasci_transition(configuration, reference, n_electron, max_jump):
            p += configuration['amplitude']**2

    return p

def is_transition(configuration, reference, n_electron=1, max_jump=10):
    """
    Determine if a configuration corresponds to a transition of n_electron

    :param configuration: dictionary containing the configuration to be analyzed
    :param reference: reference configuration (in general lowest energy Slater determinant)
    :param n_electron: number of electrons in the transition
    :param max_jump: Restrict to transitions with jumps less or equal to max_jump orbitals

    :return: True if conditions are met, otherwise False
    """

    alpha_diff = [int(i) - int(j) for i, j in zip(configuration['occupations']['alpha'], reference['alpha'])]
    beta_diff = [int(i) - int(j) for i, j in zip(configuration['occupations']['beta'], reference['beta'])]


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


def get_ratio_of_condition(state, n_electron=1, max_jump=10, ground_state_configuration=None):
    # reference = {'hole': '', 'alpha': '111000', 'beta': '111000', 'part': '', 'amplitude': 0.9777044}

    n_alpha = np.sum(state['configurations'][0]['occupations']['alpha'])
    n_beta = np.sum(state['configurations'][0]['occupations']['beta'])
    n_orb = len(state['configurations'][0]['occupations']['alpha'])

    if ground_state_configuration is None:
        reference = {'alpha': [1] * n_alpha + [0] * (n_orb - n_alpha),
                     'beta': [1] * n_beta + [0] * (n_orb - n_beta)}
    else:
        reference = ground_state_configuration

    p = 0
    for configuration in state['configurations']:
        if is_transition(configuration, reference, n_electron, max_jump):
            p += configuration['amplitude']**2

    return p


def _set_zero_to_coefficients(basis, mo_coeff, range_atoms):
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


def crop_electronic_structure(electronic_structure, atom_list):
    n_atoms = electronic_structure['structure'].get_number_of_atoms()

    complementary_list = [i for i in range(n_atoms) if not i in atom_list]

    new_coefficients = _set_zero_to_coefficients(electronic_structure['basis'],
                                                 electronic_structure['coefficients'],
                                                 complementary_list)

    new_electronic_structure = deepcopy(electronic_structure)
    new_electronic_structure['coefficients'] = new_coefficients
    return new_electronic_structure


def get_plane(coords, direction=None):
    """
    Returns the center and normal vector of tha plane formed by a list of atomic coordinates

    :param coords: List of atomic coordinates

    :return: center, normal_vector
    """

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
    # works for closed shell only
    alpha_e = np.sum([int(c) for c in configuration['alpha']])
    beta_e = np.sum([int(c) for c in configuration['beta']])
    hole = 0 if configuration['hole'] == '' else 1
    part = 0 if configuration['part'] == '' else 1

    return (structure.number_of_electrons + structure.charge - (alpha_e + beta_e + part - hole))//2


def get_inertia(structure):
    """
    returns the inertia moments and main axis of inertia (in rows)

    :param structure: Structure object containg the molecule

    :return: eigenvalues, eigenvectors
    """

    coordinates = structure.get_coordinates()
    masses = structure.get_atomic_masses()

    coordinates = np.array(coordinates)

    cm = np.average(coordinates, axis=0, weights=masses)

    inertia_tensor = np.zeros((3, 3))
    for c, m in zip(coordinates, masses):
        inertia_tensor += m * (np.dot(c-cm, c-cm) * np.identity(3) - np.outer(c-cm, c-cm))

    eval, ev = np.linalg.eigh(inertia_tensor)

    return eval.tolist(), ev.T.tolist()


def get_basis_functions_ranges_by_atoms(basis, atoms_range=None):

    num = 0
    functions_range = []
    for atoms in np.array(basis['atoms']):
        ini = num
        for shell in atoms['shells']:
            num += shell['functions']
        fin = num
        functions_range.append((ini, fin))

    if atoms_range is not None:
        functions_range = np.array(functions_range)[atoms_range]

    return functions_range

def classify_diabatic_states_of_fragment(diabatic_states, fragments_atoms, tol=0.1):

    print('     Attach      Detach')

    types = []
    for i, state in enumerate(diabatic_states):

        sum_attach = np.sum([state['mulliken']['attach'][i] for i in fragments_atoms])
        sum_detach = np.sum([state['mulliken']['detach'][i] for i in fragments_atoms])

        type = 'unknown'
        if np.abs(sum_attach) > 1.0 - tol:
            if np.abs(sum_detach) < tol:
                type = 'CT+'
            else:
                type = 'LE'
        elif np.abs(sum_attach) < tol:
            if np.abs(sum_detach) > 1 - tol:
                type = 'CT-'

        print('{} {:10.5f} {:10.5f} {}'.format(i + 1, sum_attach, sum_detach, type))
        types.append(type)

    return types


def get_occupated_list(configuration, structure, total_orbitals):
    import numpy as np
    occupied_orbitals = get_occupied_electrons(configuration, structure)
    n_extra = total_orbitals - occupied_orbitals - len(configuration['alpha'])
    vector_alpha = [1] * occupied_orbitals + [int(c) for c in configuration['alpha']] + [0] * n_extra

    n_extra = total_orbitals - occupied_orbitals - len(configuration['beta'])
    vector_beta = [1] * occupied_orbitals + [int(c) for c in configuration['beta']] + [0] * n_extra

    if configuration['hole'] != '':
        if np.sum(vector_alpha) > np.sum(vector_beta):
            vector_alpha[int(configuration['hole']) - 1] = 0
        else:
            vector_beta[int(configuration['hole']) - 1] = 0

    if configuration['part'] != '':
        if np.sum(vector_alpha) < np.sum(vector_beta):
            vector_alpha[int(configuration['part']) - 1] = 1
        else:
            vector_beta[int(configuration['part']) - 1] = 1

    return {'alpha': vector_alpha, 'beta': vector_beta}



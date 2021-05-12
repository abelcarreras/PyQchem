import numpy as np
import re
from pyqchem.utils import get_occupied_electrons


def read_basic_info(output):
    enum = output.find('Molecular Point Group')
    mpg = output[enum:enum+100].split()[3]
    enum = output.find('Largest Abelian Subgroup')
    las = output[enum:enum+100].split()[3]

    there_vector = [m.start() for m in re.finditer('There are ', output)]
    n_alpha = int(output[there_vector[0]:there_vector[0]+100].split()[2])
    n_beta = int(output[there_vector[0]:there_vector[0]+100].split()[5])

    nshell = int(output[there_vector[1]:there_vector[1]+100].split()[2])
    nbas = int(output[there_vector[1]:there_vector[1]+100].split()[5])

    return {'molecular_point_group': mpg,
            'largest_abelian_subgroup': las,
            'n_alpha': n_alpha,
            'n_beta': n_beta,
            'n_shells': nshell,
            'n_basis_functions': nbas}


def get_cis_occupations_list(number_of_orbitals,
                             alpha_electrons,
                             beta_electrons,
                             alpha_transitions=(),
                             beta_transitions=(),
                             ground_state=None):

    if ground_state is None:
        alpha_occupation = [1] * alpha_electrons + (number_of_orbitals - alpha_electrons) * [0]
        beta_occupation = [1] * beta_electrons + (number_of_orbitals - beta_electrons) * [0]
    else:
        alpha_occupation = ground_state['alpha']
        beta_occupation = ground_state['beta']

    for transition in alpha_transitions:
        alpha_occupation[transition['origin'] - 1] = 0
        alpha_occupation[transition['target'] - 1] = 1

    for transition in beta_transitions:
        beta_occupation[transition['origin'] - 1] = 0
        beta_occupation[transition['target'] - 1] = 1

    return {'alpha': alpha_occupation,
            'beta': beta_occupation}


def get_rasci_occupations_list(configuration, structure, total_orbitals):
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


def search_bars(output, from_position=0, bar_type='---'):
    output = output[from_position:]
    positions = []
    previous = 0
    for m in re.finditer(bar_type, output):
        if m.start() > previous + 1:
            positions.append(m.start() + from_position)
        previous = m.end()

    return positions


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


import re
import numpy as np


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


def search_bars(output, from_position=0):
    output = output[from_position:]
    positions = []
    previous = 0
    for m in re.finditer('---', output):
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

    alpha_int = [int(i) for i in state['amplitudes'][0]['alpha']]
    ref_alpha = ''
    for i in range(len(alpha_int)):
        if i < np.sum(alpha_int):
            ref_alpha += '1'
        else:
            ref_alpha += '0'

    reference = {'hole': '', 'alpha': ref_alpha, 'beta': ref_alpha, 'part': ''}

    p = 0
    for configuration in state['amplitudes']:
        if is_transition(configuration, reference, n_electron, max_jump):
            p += configuration['amplitude']**2

    return p


if __name__ == '__main__':
    state = {'amplitudes': [{'hole': '', 'alpha': '110100', 'beta': '111000', 'part': '', 'amplitude': 0.5}]}
    b = get_ratio_of_condition(state)
    print('test', b)
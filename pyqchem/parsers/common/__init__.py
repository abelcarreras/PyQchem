from pyqchem.utils import get_occupied_electrons
from pyqchem.structure import atom_data
from pyqchem.structure import Structure
import numpy as np
import re
import yaml


def read_symmetry_info(output):
    enum = output.find('Molecular Point Group')
    mpg = output[enum:enum + 100].split()[3]
    enum = output.find('Largest Abelian Subgroup')
    if enum > 0:
        las = output[enum:enum + 100].split()[3]
    else:
        las = None

    return {'molecular_point_group': mpg,
            'largest_abelian_subgroup': las}


def read_scf_info(output):
    # scf_energy

    enum = output.find('Total energy in the final basis set') # 6.0 support
    if enum < 0:
        enum = output.find('Total energy =') # >6.3 support
        scf_energy = float(output[enum:enum+100].split()[3])
    else:
        scf_energy = float(output[enum:enum+100].split()[8])

    return {'scf_energy': scf_energy}

def read_basic_info(output):

    there_vector = [m.start() for m in re.finditer('There are ', output)]
    n_alpha = int(output[there_vector[0]:there_vector[0]+100].split()[2])
    n_beta = int(output[there_vector[0]:there_vector[0]+100].split()[5])

    nshell = int(output[there_vector[1]:there_vector[1]+100].split()[2])
    nbas = int(output[there_vector[1]:there_vector[1]+100].split()[5])

    return {'n_alpha': n_alpha,
            'n_beta': n_beta,
            'n_shells': nshell,
            'n_basis_functions': nbas}

def read_input_structure(output):

    enum = output.find('Standard Nuclear Orientation')
    end_section = search_bars(output, from_position=enum, bar_type=r'-----')
    section_structure = output[end_section[0]: end_section[1]].split('\n')

    symbols = []
    coordinates = []
    for line in section_structure[1:-1]:
        symbols.append(line.split()[1])
        coordinates.append(line.split()[2:5])

    coordinates = np.array(coordinates, dtype=float).tolist()

    # basic info
    enum = output.find('Nuclear Repulsion Energy')
    basic_data = read_basic_info(output[enum:enum + 5000])

    n_nucleus = 0
    for s in symbols:
        for i, row in enumerate(atom_data):
            if row[1] == s:
                n_nucleus += i

    charge = n_nucleus - (basic_data['n_alpha'] + basic_data['n_beta'])
    multiplicity = abs(basic_data['n_alpha'] - basic_data['n_beta']) + 1

    return Structure(coordinates=coordinates,
                     symbols=symbols,
                     charge=charge,
                     multiplicity=multiplicity)


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


def get_rasci_occupations_list(configuration, occupied_orbitals, total_orbitals):
    # occupied_orbitals = get_occupied_electrons(configuration, structure)
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


# handle asterisks and convert them to nan
class float_asterisk(float):
    def __new__(cls, value):
        if isinstance(value, str) and '**' in value:
            return super(float_asterisk, cls).__new__(cls, 'nan')
        value = str(float(value)+0)
        return super(float_asterisk, cls).__new__(cls, value)

    @classmethod
    def to_yaml(cls, dumper, data):
        return dumper.represent_float(float(data))


# Representer register in yaml
yaml.add_representer(float_asterisk, float_asterisk.to_yaml)




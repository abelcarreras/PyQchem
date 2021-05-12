__author__ = 'Abel Carreras'
AU_TO_EV = 27.21138

from pyqchem.structure import Structure
from pyqchem.errors import ParserError
from pyqchem.parsers.common import search_bars, standardize_vector
from pyqchem.parsers.common import read_basic_info, get_cis_occupations_list
import numpy as np
import re


def _list_to_complex(list):
    real = list[0].replace('(', '').replace(')', '')
    if real[:2] == '--':
        real = real[2:]

    imag = list[2].split('i')[0].replace('(', '')
    if imag[:2] == '--':
        imag = imag[2:]

    opera = float(list[1] + '1')
    return float(real) + opera * float(imag) * 1.j


def basic_cis(output):
    """
    Parser for CIS/TD-DFT calculations

    :param output:
    :return:
    """
    data_dict = {}

    # Molecule
    n = output.find('$molecule')
    n2 = output[n:].find('$end')

    molecule_region = output[n:n+n2-1].replace('\t', ' ').split('\n')[1:]
    charge, multiplicity = [int(num) for num in molecule_region[0].split()]
    coordinates = [[float(l) for l in line.split()[1:4]] for line in molecule_region[1:]]
    symbols = [line.split()[0].capitalize() for line in molecule_region[1:]]
    n_atoms = len(symbols)

    # structure
    structure_input = Structure(coordinates=coordinates,
                                symbols=symbols,
                                charge=charge,
                                multiplicity=multiplicity)

    enum = output.find('Standard Nuclear Orientation')
    section_structure = output[enum:enum + 200*structure_input.get_number_of_atoms()].split('\n')
    section_structure = section_structure[3:structure_input.get_number_of_atoms()+3]
    coordinates = [[float(num) for num in s.split()[2:]] for s in section_structure]

    data_dict['structure'] = Structure(coordinates=coordinates,
                                       symbols=symbols,
                                       charge=charge,
                                       multiplicity=multiplicity)

    # scf_energy
    enum = output.find('Total energy in the final basis set')
    try:
        data_dict['scf_energy'] = float(output[enum:enum+100].split()[8])
    except IndexError:
        pass

    enum = output.find('Molecular Point Group')
    basic_data = read_basic_info(output[enum:enum + 5000])

    # CIS excited states
    # enum = output.find('CIS Excitation Energies')
    try:
        enum = list(re.finditer('CIS Excitation Energies', output))[-1].end()
    except IndexError:
        enum = list(re.finditer('TDDFT/TDA Excitation Energies', output))[-1].end()

    excited_states = []
    if enum > 0:
        bars = search_bars(output, from_position=enum)

        output_cis = output[bars[0]:bars[1]]

        for m in re.finditer('Excited state ', output_cis):
            state_cis_section = output_cis[m.end():]
            state_cis_lines = state_cis_section.split('\n')

            exc_energy = float(state_cis_lines[0].split()[5])
            exc_energy_units = state_cis_lines[0].split()[3][1:-1]
            tot_energy = float(state_cis_lines[1].split()[5])

            try:
                tot_energy_units = state_cis_lines[1].split()[6]
                mul = state_cis_lines[2].split()[-1]

                trans_mom = [float(mom) for mom in [state_cis_lines[3].split()[2],
                                                    state_cis_lines[3].split()[4],
                                                    state_cis_lines[3].split()[6]]]
                strength = float(state_cis_lines[4].split()[2])
            except ValueError:
                # old version of qchem (< 5.01)
                state_cis_words = output_cis[m.end():].split()
                tot_energy_units = 'au'
                mul = state_cis_words[13]
                trans_mom = [float(mom) for mom in [state_cis_words[16],
                                                    state_cis_words[18],
                                                    state_cis_words[20]]]
                strength = float(state_cis_words[24])

            transitions = []
            for line in state_cis_lines[5:]:
                if line.find('-->') > 0:
                    origin = int(line.split('>')[0].split('(')[1].split(')')[0])
                    target = int(line.split('>')[1].split('(')[1].split(')')[0])
                    amplitude = float(line.split('=')[1])

                    alpha_transitions = []
                    beta_transitions = []
                    try:
                        spin = line[21:].split()[3]
                        if spin == 'alpha':
                            alpha_transitions.append({'origin': origin, 'target': target + basic_data['n_alpha']})
                        elif spin == 'beta':
                            beta_transitions.append({'origin': origin, 'target': target + basic_data['n_beta']})
                        else:
                            raise ParserError('basic_cis', 'Error reading configurations')

                        transitions.append({'origin': origin,
                                            'target': target,
                                            'amplitude': amplitude,
                                            'occupations': get_cis_occupations_list(basic_data['n_basis_functions'],
                                                                                    basic_data['n_alpha'],
                                                                                    basic_data['n_beta'],
                                                                                    alpha_transitions=alpha_transitions,
                                                                                    beta_transitions=beta_transitions)})

                    except (IndexError, ParserError):
                        # This supposes single electron transition
                        alpha_transitions.append({'origin': origin, 'target': target + basic_data['n_alpha']})

                        transitions.append({'origin': origin,
                                            'target': target,
                                            'amplitude': amplitude/np.sqrt(2),
                                            'occupations': get_cis_occupations_list(basic_data['n_basis_functions'],
                                                                                    basic_data['n_alpha'],
                                                                                    basic_data['n_beta'],
                                                                                    alpha_transitions=alpha_transitions,
                                                                                    beta_transitions=beta_transitions)})

                        transitions.append({'origin': origin,
                                            'target': target,
                                            'amplitude': amplitude/np.sqrt(2) if mul == 'Singlet' else -amplitude/np.sqrt(2),
                                            'occupations': get_cis_occupations_list(basic_data['n_basis_functions'],
                                                                                    basic_data['n_alpha'],
                                                                                    basic_data['n_beta'],
                                                                                    alpha_transitions=beta_transitions,
                                                                                    beta_transitions=alpha_transitions)})

                if len(line) < 5:
                    break

            excited_states.append({'total_energy': tot_energy,
                                   'total_energy_units': tot_energy_units,
                                   'excitation_energy': exc_energy,
                                   'excitation_energy_units': exc_energy_units,
                                   'multiplicity': mul,
                                   'transition_moment': standardize_vector(trans_mom),
                                   'strength': strength,
                                   'configurations': transitions})

    data_dict['excited_states'] = excited_states

    # Spin-Orbit coupling
    initial = output.find('*********SPIN-ORBIT COUPLING JOB BEGINS HERE*********')
    final = output.find('*********SOC CODE ENDS HERE*********')

    data_interstate = {}
    if initial > 0:
        soc_section = output[initial:final]

        def label_states(excited_states):
            labels = []
            ns = 1
            nt = 1
            for state in excited_states:
                if state['multiplicity'].lower() == 'singlet':
                    labels.append('S{}'.format(ns))
                    ns += 1
                elif state['multiplicity'].lower() == 'triplet':
                    labels.append('T{}'.format(nt))
                    nt += 1
                else:
                    try:
                        m = float(state['multiplicity'])
                        if abs(m - 1) < 0.1:
                            labels.append('S{}'.format(ns))
                            ns += 1
                        if abs(m - 3) < 0.1:
                            labels.append('T{}'.format(nt))
                            nt += 1
                        state['multiplicity'] = m

                    except ValueError:
                        raise ParserError('basic_cis', 'State multiplicity error')

            return labels, nt-1, ns-1

        labels, n_triplet, n_singlet = label_states(excited_states)

        for i, label in enumerate(labels):
            data_interstate[(i+1, 0)] = {'1e_soc_mat': [0j, 0j, 0j], 'soc_units': 'cm-1'}
            data_interstate[(0, i+1)] = {'1e_soc_mat': [0j, 0j, 0j], 'soc_units': 'cm-1'}
            for j, label2 in enumerate(labels):
                if (label[0] == 'S' or label2[0] == 'S') and (label[0] != label2[0]):
                    data_interstate[(i+1, j+1)] = {'1e_soc_mat': [[0j, 0j, 0j]], 'soc_units': 'cm-1'}
                elif label[0] == 'T' and label2[0] == 'T':
                    data_interstate[(i+1, j+1)] = {'1e_soc_mat': [[0j, 0j, 0j], [0j, 0j, 0j], [0j, 0j, 0j]], 'soc_units': 'cm-1'}
                elif label[0] == 'S' and label2[0] == 'S':
                    data_interstate[(i+1, j+1)] = {'1e_soc_mat': [[0j]], 'soc_units': 'cm-1'}
                else:
                    raise ParserError('basic_cis', 'State multiplicity error')

        for i, label in enumerate(labels):
            for k2, ms2 in enumerate([-1, 0, 1]):
                for j, label2 in enumerate(labels):
                    if label[0] == 'T':
                        for k, ms in enumerate([-1, 0, 1]):
                            enum = soc_section.find('SOC between the {} (ms={}) state and excited triplet states (ms={})'.format(label, ms2, ms))
                            for line in soc_section[enum:enum+50*(n_triplet+1)].split('\n'):
                                if len(line.split()) == 0:
                                    break
                                if line.split()[0] == '{}(ms={})'.format(label2, ms):
                                    data_interstate[(i+1, j+1)]['1e_soc_mat'][k2][k] = _list_to_complex(line.split()[1:4])
                                    data_interstate[(i+1, j+1)]['1e_soc_mat'][k][k2] = _list_to_complex(line.split()[1:4])
                                    data_interstate[(j+1, i+1)]['1e_soc_mat'][k2][k] = _list_to_complex(line.split()[1:4])
                                    data_interstate[(j+1, i+1)]['1e_soc_mat'][k][k2] = _list_to_complex(line.split()[1:4])
                                    break

                    elif label[0] == 'S':
                        for k, ms in enumerate([-1, 0, 1]):
                            enum = soc_section.find('SOC between the {} state and excited triplet states (ms={})'.format(label, ms))
                            for line in soc_section[enum:enum+50*(n_triplet+1)].split('\n'):
                                if len(line.split()) == 0:
                                    break
                                if line.split()[0] == '{}(ms={})'.format(label2, ms):
                                    data_interstate[(i+1, j+1)]['1e_soc_mat'][0][k] = _list_to_complex(line.split()[1:4])
                                    data_interstate[(j+1, i+1)]['1e_soc_mat'][0][k] = _list_to_complex(line.split()[1:4])
                                    break
                    else:
                        raise ParserError('basic_cis', 'SOC reading error')

                enum = soc_section.find('SOC between the singlet ground state and excited triplet states (ms={})'.format(ms2))
                for line in soc_section[enum:enum+50*(n_triplet+1)].split('\n'):
                    if len(line.split()) == 0:
                        break
                    if line.split()[0] == '{}(ms={})'.format(label, ms2):
                        data_interstate[(i+1, 0)]['1e_soc_mat'][k2] = _list_to_complex(line.split()[1:4])
                        data_interstate[(0, i+1)]['1e_soc_mat'][k2] = _list_to_complex(line.split()[1:4])
                        break

        data_dict['interstate_properties'] = data_interstate

    # diabatization
    initial = output.find('Localization Code for CIS excited states')
    if initial > 0:

        bars = search_bars(output, from_position=initial)
        diabat_section = output[initial: bars[0]]

        def read_diabatization_matrix(label):
            matrix = []
            for m in re.finditer(label, diabat_section):
                line = diabat_section[m.end(): m.end() + 50].split('\n')[0]
                matrix.append(float(line.split('=')[1]))

            diabat_dim = int(np.sqrt(len(matrix)))
            return np.array(matrix).reshape(diabat_dim, diabat_dim).T

        rot_matrix = read_diabatization_matrix('showmatrix adiabatic R-Matrix')
        adiabatic_matrix = read_diabatization_matrix('showmatrix adiabatH') * AU_TO_EV
        diabatic_matrix = read_diabatization_matrix('showmatrix diabatH') * AU_TO_EV

        diabat_data = {'rot_matrix': rot_matrix,
                       'adiabatic_matrix': adiabatic_matrix.tolist(),
                       'diabatic_matrix': diabatic_matrix.tolist()}

        if diabat_section.find('showmatrix Total_Decomposed_H_diabatic'):

            tot_decomp_matrix = read_diabatization_matrix('showmatrix Total_Decomposed_H_diabatic') * AU_TO_EV
            decomp_one_matrix = read_diabatization_matrix('showmatrix Decomposed_One_diabatic') * AU_TO_EV
            decomp_j_matrix = read_diabatization_matrix('showmatrix Decomposed_J_diabatic') * AU_TO_EV
            decomp_k_matrix = read_diabatization_matrix('showmatrix Decomposed_K_diabatic') * AU_TO_EV

            diabat_data.update({'tot_decomp_matrix': tot_decomp_matrix,
                                               'decomp_one_matrix': decomp_one_matrix.tolist(),
                                               'decomp_j_matrix': decomp_j_matrix.tolist(),
                                               'decomp_k_matrix': decomp_k_matrix.tolist()})

        mulliken_diabatic = []

        enum = output.find('Mulliken & Loewdin analysis of')
        for m in re.finditer('Mulliken analysis of TDA State', output[enum:]):
            section_mulliken = output[m.end() + enum: m.end() + 10000 + enum]  # 10000: assumed to max of section
            section_mulliken = section_mulliken[:section_mulliken.find('Natural Orbitals stored in FCHK')]
            section_attachment = section_mulliken.split('\n')[10 + n_atoms: 10 + n_atoms * 2]

            mulliken_diabatic.append({'attach': [float(l.split()[1]) for l in section_attachment],
                                      'detach': [float(l.split()[2]) for l in section_attachment],
                                      'total': [float(l.split()[3]) for l in section_attachment]})

        diabatic_states = []
        for i in range(len(rot_matrix)):
            diabat_states_data = {'excitation_energy': diabatic_matrix[i][i],
                                  'excitation_energy_units': 'eV',
                                  'transition_moment': [],
                                  'dipole_moment_units': 'ua'}
            if len(mulliken_diabatic) > 0:
                diabat_states_data['mulliken'] = mulliken_diabatic[i]

            diabatic_states.append(diabat_states_data)
        diabat_data['diabatic_states'] = diabatic_states

        data_dict['diabatization'] = diabat_data

    return data_dict


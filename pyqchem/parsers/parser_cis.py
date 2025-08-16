__author__ = 'Abel Carreras'
AU_TO_EV = 27.21138

from pyqchem.structure import Structure
from pyqchem.errors import ParserError
from pyqchem.parsers.common import search_bars, standardize_vector, read_scf_info
from pyqchem.parsers.common import read_basic_info, get_cis_occupations_list, read_symmetry_info, read_input_structure
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

    # input echo data

    # get diabatic states numbers (starting by zero)
    re_ini = re.search('\$localized_diabatization', output, re.IGNORECASE)
    list_diabat = []
    if re_ini is not None:
        enum_i = re_ini.end()
        enum_f = re.search('\$end', output[enum_i: enum_i+1000], re.IGNORECASE).start()
        diabat_input_section = output[enum_i: enum_i+enum_f]
        enum = re.search('adiabatic states', diabat_input_section, re.IGNORECASE).end()
        list_diabat = [int(val) for val in diabat_input_section[enum:].split('\n')[1].split()]

    # Molecule
    data_dict['structure'] = read_input_structure(output)
    n_atoms = data_dict['structure'].get_number_of_atoms()
    
    # info
    data_dict.update(read_scf_info(output))

    enum = output.find('Molecular Point Group')
    if enum > 0:
        symmetry_data = read_symmetry_info(output[enum:enum + 1000])

    enum = output.find('Nuclear Repulsion Energy')
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

            try:
                tot_energy = float(state_cis_lines[1].split()[5])
            except ValueError:
                tot_energy = float(state_cis_lines[1].split()[4])

            tot_energy_units = 'au'

            try:
                mul = state_cis_lines[2].split()[-1]

                trans_mom = [float(mom) for mom in [state_cis_lines[3].split()[2],
                                                    state_cis_lines[3].split()[4],
                                                    state_cis_lines[3].split()[6]]]
                strength = float(state_cis_lines[4].split()[2])
            except ValueError:
                # old version of qchem (< 5.01)
                state_cis_words = output_cis[m.end():].split()
                mul = state_cis_words[13]
                print(mul)
                trans_mom = [float(mom) for mom in [state_cis_words[16],
                                                    state_cis_words[18],
                                                    state_cis_words[20]]]
                strength = float(state_cis_words[24])

            transitions = []
            for line in state_cis_lines[5:]:
                if line.find('-->') > 0:
                    origin = int(line.split('>')[0].split('(')[1].split(')')[0])
                    target = int(line.split('>')[1].split('(')[1].split(')')[0])
                    amplitude = float(line.split('=')[1].split()[0])

                    alpha_transitions = []
                    beta_transitions = []
                    try:
                        spin = line[21:].split()[-1]
                        if spin == 'alpha':
                            target = target + basic_data['n_alpha']
                            alpha_transitions.append({'origin': origin, 'target': target})
                        elif spin == 'beta':
                            target = target + basic_data['n_beta']
                            beta_transitions.append({'origin': origin, 'target': target})
                        else:
                            raise ParserError('basic_cis', 'Error reading configurations', output)

                        transitions.append({'spin': spin,
                                            'origin': origin,
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

                        transitions.append({'spin': 'alpha',
                                            'origin': origin,
                                            'target': target + basic_data['n_alpha'],
                                            'amplitude': amplitude/np.sqrt(2),
                                            'occupations': get_cis_occupations_list(basic_data['n_basis_functions'],
                                                                                    basic_data['n_alpha'],
                                                                                    basic_data['n_beta'],
                                                                                    alpha_transitions=alpha_transitions,
                                                                                    beta_transitions=beta_transitions)})

                        transitions.append({'spin': 'beta',
                                            'origin': origin,
                                            'target': target + basic_data['n_beta'],
                                            'amplitude': amplitude/np.sqrt(2) if mul == 'Singlet' else -amplitude/np.sqrt(2),
                                            'occupations': get_cis_occupations_list(basic_data['n_basis_functions'],
                                                                                    basic_data['n_alpha'],
                                                                                    basic_data['n_beta'],
                                                                                    alpha_transitions=beta_transitions,
                                                                                    beta_transitions=alpha_transitions)})

                if len(line) < 5:
                    break

            # sort transitions by amplitude
            sorted_transitions = sorted(transitions, key=lambda x: abs(x['amplitude']), reverse=True)

            excited_states.append({'total_energy': tot_energy,
                                   'total_energy_units': tot_energy_units,
                                   'excitation_energy': exc_energy,
                                   'excitation_energy_units': exc_energy_units,
                                   'multiplicity': mul,
                                   'transition_moment': standardize_vector(trans_mom),
                                   'oscillator_strength': strength,
                                   'configurations': sorted_transitions})

    data_dict['excited_states'] = excited_states

    # Multiplicity mapping for singlets, triplets... where multiplicity is indicated as "S1", "T1" in SOC section
    
    try: 
        test = float(data_dict['excited_states'][0]['multiplicity'])
    except ValueError:
        multiplicity_mapping = {'0': 0}
        count_multiplicity = {'Singlet': 0, 'Triplet': 0, 'Quintet': 0, 'Heptet': 0}
        for i, state in enumerate(data_dict['excited_states']):
            count_multiplicity[state['multiplicity']] += 1
            multiplicity_mapping.update({state['multiplicity'][0]+str(count_multiplicity[state['multiplicity']]): i+1})

    # SPIN-ORBIT COUPLINGS

    done_interstate_soc = bool(output.find('SPIN-ORBIT COUPLING JOB USING WIGNER–ECKART THEOREM BEGINS HERE')+1)

    if done_interstate_soc:
        ini_section = output.find('SPIN-ORBIT COUPLING JOB USING WIGNER–ECKART THEOREM BEGINS HERE')
        end_section = output.find('SPIN-ORBIT COUPLING JOB ENDS HERE')
        soc_section = output[ini_section: end_section]
        section_sizes = search_bars(soc_section, bar_type=r'\*'*50+'\n')

        interstate_dict = {}

        for i, m in enumerate(re.finditer('State A:', soc_section)):
            section_length = section_sizes[i+1] - section_sizes[i]
            section_pair = soc_section[m.start():m.start() + section_length]
            
            lines = section_pair.split('\n')
            
            try: # In doublets, where states can be taken as integers (0, 1..)
                state_a = int(lines[0].split()[-1])
                state_b = int(lines[1].split()[-1])

            except ValueError: # In singlets, where states cannot be taken as integers (S1, T1...)
                if 'Ground state' in lines[0]:
                    state_a = 0
                else:
                    state_a = multiplicity_mapping[lines[0].split()[-1]]
                state_b = multiplicity_mapping[lines[1].split()[-1]]

            pair_dict = {}
            for j, line in enumerate(lines):
                if 'Ket state' in line and state_a == 0 and 'scf_multiplicity' not in data_dict:
                        data_dict['scf_multiplicity'] = int(float(line.split()[5]))
                
                if 'WARNING! Clebsh-Gordon coefficient is too small' in line:
                    pair_dict['1e_socc'] = 0
                    pair_dict['1e_soc_mat'] = [0]
                    pair_dict['mf_socc'] = 0
                    pair_dict['mf_soc_mat'] = [0]

                if 'One-electron SO (cm-1)' in line: 
                    soc_keyword = '1e_socc'
                    soc_matrix_keyword = '1e_soc_mat'
                elif 'Mean-field SO (cm-1)' in line:
                    soc_keyword = 'mf_socc'
                    soc_matrix_keyword = 'mf_soc_mat'                        

                if 'SOCC' in line:
                    pair_dict[soc_keyword] = float(line.split()[2])

                if 'Actual matrix elements:' in line:
                    soc_matrix = []
                    for soc_line in lines[j:]:

                        if '<Sz=' in soc_line:
                            matches = re.findall(r"\(([^)]+)\)", soc_line)
                            soc_matrix.append([list(map(float, match.split(','))) for match in matches])

                        if r'_'*30 in soc_line:
                            break
                    
                    pair_dict[soc_matrix_keyword] = soc_matrix

            interstate_dict[(state_a, state_b)] = pair_dict
    
        data_dict['interstate_socs'] = interstate_dict

    # ANGULAR MOMENTUMS

    done_interstate_angmom = bool(output.find(' STATE-TO-STATE TRANSITION ANGULAR MOMENTS')+1)

    if done_interstate_angmom:
        ini_section = output.find('STATE-TO-STATE TRANSITION ANGULAR MOMENTS')
        end_section = output.find('END OF TRANSITION MOMENT CALCULATION')
        angmom_section = output[ini_section: end_section]
        section_sizes = search_bars(angmom_section, bar_type='Transition Angular Moments') # r'-'*40)

        interstate_dict = {}

        for i, m in enumerate(re.finditer('Transition Angular Moments', angmom_section)):
            try: 
                section_length = section_sizes[i+1] - section_sizes[i]
            except IndexError:
                section_length = len(angmom_section) - section_sizes[i]
            section_pair = angmom_section[m.start():m.start() + section_length]
            
            lines = section_pair.split('\n')

            for line in lines[4:]:
                if r'-'*40 in line:
                    break

                pair_dict = {}
                state_a = int(line.split()[0])
                state_b = int(line.split()[1])
                
                pair_dict['angular_momentum'] = [float(j) for j in line.split()[2:-1]]
                pair_dict['length'] = float(line.split()[-1])

                interstate_dict[(state_a, state_b)] = pair_dict
    
        data_dict['interstate_angmom'] = interstate_dict

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

        r_matrix = read_diabatization_matrix('showmatrix adiabatic R-Matrix')
        rot_matrix = read_diabatization_matrix('showmatrix final adiabatic -> diabatic RotMatrix')
        diabatic_matrix = read_diabatization_matrix('showmatrix diabatH') * AU_TO_EV
        if len(rot_matrix) == 0:
            # Boys diabatization
            rot_matrix = read_diabatization_matrix('showmatrix Boys adiabatic->diabatic RotMatrix')
            diabatic_matrix = read_diabatization_matrix('showmatrix Boys diabatH') * AU_TO_EV

        adiabatic_matrix = read_diabatization_matrix('showmatrix adiabatH') * AU_TO_EV

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

        # Mulliken analysis
        mulliken_analysis = []
        enum = output.find('Mulliken & Loewdin analysis of')
        for m in re.finditer('Mulliken analysis of TDA State', output[enum:]):
            section_mulliken = output[m.end() + enum: m.end() + 10000 + enum]  # 10000: assumed to max of section

            enum_i = search_bars(section_mulliken)[2]
            section_mulliken = section_mulliken[:enum_i]
            section_attachment = section_mulliken.split('\n')[4: 4 + n_atoms]

            mulliken_analysis.append({'attach': [float(l.split()[1]) for l in section_attachment],
                                      'detach': [float(l.split()[2]) for l in section_attachment],
                                      'total': [float(l.split()[3]) for l in section_attachment]})

        diabatic_states = []
        for i in range(len(rot_matrix)):
            diabat_states_data = {'excitation_energy': diabatic_matrix[i][i],
                                  'excitation_energy_units': 'eV',
                                  'transition_moment': [],
                                  'dipole_moment_units': 'ua'}
            if len(mulliken_analysis) > 0:
                diabat_states_data['mulliken'] = mulliken_analysis[list_diabat[i]-1]

            diabatic_states.append(diabat_states_data)
        diabat_data['diabatic_states'] = diabatic_states

        data_dict['diabatization'] = diabat_data

    return data_dict


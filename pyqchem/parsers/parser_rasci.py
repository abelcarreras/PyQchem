__author__ = 'Abel Carreras'

import re
import operator
import warnings

import numpy as np
from pyqchem.parsers.common import read_basic_info, get_rasci_occupations_list, read_scf_info
from pyqchem.parsers.common import search_bars, standardize_vector, read_input_structure
from pyqchem.parsers.common import float_asterisk as float


def _read_simple_matrix(header, output, maxchar=10000, foot='-------'):
    matrix_list = []
    for m in re.finditer(header, output):
        section_state = output[m.end():m.end() + maxchar]  # 10000: assumed to max of section
        section_state = section_state[:section_state.find(foot)]
        dim = len(section_state.split('\n')[1].split())
        matrix = section_state.split('\n')[1:dim + 1]
        matrix = [[float(n) for n in l.split()] for l in matrix]
        matrix_list.append(matrix)

    return matrix_list


def _read_soc_matrix(lines, dimensions):

    col_per_line = 5
    matrix = []
    for ib in range(dimensions[0]):
        real = []
        complex = []
        for j in range((dimensions[1] - 1) // col_per_line + 1):
            real += lines[j*dimensions[0] + 1 * (j+1) + ib][11:].split()[0::2]
            complex += lines[j*dimensions[0] + 1 * (j+1) +ib][11:].split()[1::2]

        row = [float(r) + float(c[:-1]) * 1j for r, c in zip(real, complex)]
        matrix.append(row)

    return matrix


def _complete_interstate_pairs(interstate_dict):
    additional_items = {}
    for key, value in interstate_dict.items():
        if not key[::-1] in interstate_dict:
            dict_entry = {}
            for k, v in interstate_dict[key].items():
                if k == 'state_a':
                    dict_entry.update({k: key[1]})
                elif k == 'state_b':
                    dict_entry.update({k: key[0]})
                elif k == 'angular_momentum':
                    dict_entry.update({k: np.conjugate(v).tolist()})
                elif k == '1e_soc_mat':
                    dict_entry.update({k: np.conjugate(v).T.tolist()})
                elif k == 'hso_l-':
                    dict_entry.update({k: (-np.array(v)).tolist()})
                elif k == 'hso_l+':
                    dict_entry.update({k: (-np.array(v)).tolist()})
                elif k == '2e_soc_mat':
                    dict_entry.update({k: np.conjugate(v).T.tolist()})
                elif k == 'total_soc_mat':
                    dict_entry.update({k: np.conjugate(v).T.tolist()})
                else:
                    dict_entry.update({k: v})
            additional_items.update({key[::-1]: dict_entry})

    interstate_dict.update(additional_items)


def parser_rasci(output):
    """
    Parser for RAS-CI calculations
    Include:
    - Diabatization scheme data
    - Structure
    - Adiabatic states
    - SOC

    :param output:
    :return:
    """

    data_dict = {}

    # Molecule
    data_dict['structure'] = read_input_structure(output)
    n_atoms = data_dict['structure'].get_number_of_atoms()

    # basic info
    enum = output.find('Nuclear Repulsion Energy')
    basic_data = read_basic_info(output[enum:enum + 5000])

    # scf_info
    data_dict.update(read_scf_info(output))

    # total energy
    # enum = output.find('Total energy in the final basis set')
    # total_energy = float(output[enum:enum+100].split()[8])

    # RASCI dimensions
    ini_section = output.find('RAS-CI Dimensions')
    end_section = search_bars(output, from_position=ini_section, bar_type=r'\*\*\*')[0]
    dimension_section = output[ini_section: end_section]

    enum = dimension_section.find('Active Elec')
    active_elec = int(dimension_section[enum: enum+50].split()[2])
    enum = dimension_section.find('Active Orb')
    active_orb = int(dimension_section[enum: enum+50].split()[2])
    enum = dimension_section.find('Doubly Occ')
    doubly_occ = int(dimension_section[enum: enum+50].split()[3])
    enum = dimension_section.find('Doubly Vir')
    doubly_vir = int(dimension_section[enum: enum+50].split()[2])
    enum = dimension_section.find('Frozen Occ')
    frozen_occ = int(dimension_section[enum: enum+50].split()[3])
    enum = dimension_section.find('Frozen Vir')
    frozen_vir = int(dimension_section[enum: enum+50].split()[2])

    enum = dimension_section.find('Total CI configurations')
    total_conf = int(dimension_section[enum: enum+50].split()[3])
    enum = dimension_section.find('Active configurations')
    active_conf = int(dimension_section[enum: enum+50].split()[2])
    enum = dimension_section.find('Hole configurations')
    hole_conf = int(dimension_section[enum: enum+50].split()[2])
    enum = dimension_section.find('Particle configurations')
    particle_conf = int(dimension_section[enum: enum+50].split()[2])

    rasci_dimensions = {'active_elec': active_elec,
                        'active_orb': active_orb,
                        'doubly_occupied': doubly_occ,
                        'doubly_virtual': doubly_vir,
                        'frozen_occupied': frozen_occ,
                        'frozen_virtual': frozen_vir,
                        'total_configurations': total_conf,
                        'active_configurations': active_conf,
                        'hole_configurations': hole_conf,
                        'particle_configurations': particle_conf}

    data_dict.update({'rasci_dimensions': rasci_dimensions})

    # Diabatization scheme
    done_diabat = bool(output.find('RASCI DIABATIZATION')+1)
    if done_diabat:
        rot_matrix = _read_simple_matrix('showmatrix final adiabatic -> diabatic', output)[-1]
        adiabatic_matrix = _read_simple_matrix('showing H in adiabatic representation: NO coupling elements', output)[-1]
        diabatic_matrix = _read_simple_matrix('showing H in diabatic representation: WITH coupling elements', output)[-1]

        mulliken_adiabatic = []
        enum = output.find('Mulliken analysis of Adiabatic State')
        for m in re.finditer('Mulliken analysis of Adiabatic State', output[enum:]):
            end_section = search_bars(output, from_position=m.end(), bar_type=r'\-\-\-\-\-\-')[3]
            section_mulliken = output[m.end() + enum: m.end() + enum + end_section]

            #section_mulliken = output[m.end() + enum: m.end() + 20000 + enum]  # 10000: assumed to max of section
            section_mulliken = section_mulliken[:section_mulliken.find('Natural Orbitals stored in FCHK')]
            section_attachment = section_mulliken.split('\n')[9+n_atoms:9+n_atoms*2]

            mulliken_adiabatic.append({'attach': [float(l.split()[1]) for l in section_attachment],
                                       'detach': [float(l.split()[2]) for l in section_attachment],
                                       'total': [float(l.split()[3]) for l in section_attachment]})

        mulliken_diabatic = []
        enum = output.find('showing H in diabatic representation')
        for m in re.finditer('Mulliken Analysis of Diabatic State', output[enum:]):

            end_section_list = search_bars(output[enum:], from_position=m.end(), bar_type=r'\-\-\-\-\-\-')

            dict_mulliken_diabat = {}
            section_mulliken_header = output[enum + end_section_list[0]: enum + end_section_list[1]]
            if 'spin' in section_mulliken_header.lower():
                section_mulliken = output[enum + end_section_list[1]: enum + end_section_list[2]]
                section_mulliken = section_mulliken.split('\n')[1:n_atoms + 1]
                dict_mulliken_diabat.update({'charge': [float(l.split()[2]) for l in section_mulliken],
                                             'spin': [float(l.split()[3]) for l in section_mulliken]})

            section_mulliken_header = output[enum + end_section_list[4]: enum + end_section_list[5]]
            if 'h+' in section_mulliken_header.lower():
                section_mulliken = output[enum + end_section_list[5]: enum + end_section_list[6]]
                section_mulliken = section_mulliken.split('\n')[1:n_atoms + 1]
                dict_mulliken_diabat.update({'e-': [float(l.split()[1]) for l in section_mulliken],
                                             'h+': [float(l.split()[2]) for l in section_mulliken],
                                             'tot_td': [float(l.split()[3]) for l in section_mulliken]})

            section_mulliken_header = output[enum + end_section_list[7]: enum + end_section_list[8]]
            if 'attach' in section_mulliken_header.lower():
                section_mulliken = output[enum + end_section_list[8]: enum + end_section_list[9]]
                section_mulliken = section_mulliken.split('\n')[1:n_atoms + 1]

                dict_mulliken_diabat.update({'attach': [float(l.split()[1]) for l in section_mulliken],
                                             'detach': [float(l.split()[2]) for l in section_mulliken],
                                             'total': [float(l.split()[3]) for l in section_mulliken]})

            mulliken_diabatic.append(dict_mulliken_diabat)

        enum = output.find('Transition dipole moment - diabatic states')

        n_diab_states = len(rot_matrix)
        tdm_section = output[enum: enum + 100 * (n_diab_states + 1) ]

        diabatic_tdm = []
        for i in range(n_diab_states):
            diabatic_tdm.append([float(n) for n in tdm_section.split('\n')[i+2].split()[4:7]])

        diabatic_states = []
        for i, tdm in enumerate(diabatic_tdm):
            diabatic_states.append({'excitation_energy': diabatic_matrix[i][i],
                                    'excitation_energy_units': 'eV',
                                    'transition_moment': tdm,
                                    'dipole_moment_units': 'ua',
                                    'mulliken': mulliken_diabatic[i]})

        data_dict['diabatization'] = {'rot_matrix': rot_matrix,
                                      'adiabatic_matrix': adiabatic_matrix,
                                      'diabatic_matrix': diabatic_matrix,
                                      'diabatic_states': diabatic_states,
                                      'mulliken_adiabatic': mulliken_adiabatic}

    # excited states data
    excited_states = []
    for m in re.finditer('RAS-CI total energy for state', output):
        end_section = search_bars(output, from_position=m.start(), bar_type=r'\*'*20)[0]
        section_state = output[m.start():end_section]

        # energies
        enum = section_state.find('total energy for state')
        tot_energy = float(section_state[enum: enum + 50].split()[5])
        enum = section_state.find('Excitation energy')
        exc_energy_units = section_state[enum: enum + 30].split()[2].strip('(').strip(')')
        exc_energy = float(section_state[enum: enum + 80].split()[4])

        # symmetry
        enum = section_state.find('State symmetry')
        state_symmetry = section_state[enum: enum + 30].split()[2]

        # multiplicity
        n_multi = section_state.find('<S^2>')
        multi_data = section_state[n_multi:n_multi + 30].split(':')[1]
        state_multiplicity = float(multi_data.split()[0])

        # dipole moment
        enum = section_state.find('Dipole Moment')
        dipole_mom = [float(section_state[enum:].split()[2]) + 0.0,
                      float(section_state[enum:].split()[4]) + 0.0,
                      float(section_state[enum:].split()[6]) + 0.0]

        # Transition moment
        enum = section_state.find('Trans. Moment')
        trans_mom = strength = None
        if enum > -1:
            trans_mom = [float(section_state[enum:].split()[2]) + 0.0,
                         float(section_state[enum:].split()[4]) + 0.0,
                         float(section_state[enum:].split()[6]) + 0.0]
            trans_mom = standardize_vector(trans_mom)
            strength = float(section_state[enum:].split()[10])

        # Mulliken population analysis
        mulliken_population = None
        enum = section_state.find('Mulliken population analysis')
        if enum > -1:
            max = search_bars(section_state, from_position=enum, bar_type=r'\-'*30)
            pop_charges = []
            pop_spin = []
            for line in section_state[max[1]: max[2]].split('\n')[1:n_atoms+1]:
                c, s = line.split()[2:4]
                pop_charges.append(float(c))
                pop_spin.append(float(s))

            mulliken_population = {'charge': pop_charges, 'spin': pop_spin, 'units': 'au'}

        # Natural orbitals
        nato_occ = None
        enum = section_state.find('NATURAL OCCUPATION NUMBERS')
        if enum > -1:
            lines = []
            for line in section_state[enum:].split('\n')[2::2]:
                if len(line) == 0:
                    break
                if line.split()[0].isnumeric():
                    lines += line.split()[1:]
                else:
                    break
            nato_occ = [float(num) for num in lines]

        # configurations table
        enum = section_state.find('AMPLITUDE')
        enum2 = section_state.find('Contributions')
        section_table = section_state[enum: enum2].split('\n')[2:-2]

        # ' HOLE  | ALPHA | BETA  | PART | AMPLITUDE'
        table = []
        for row in section_table:
            table.append({'hole': row.split('|')[1].strip(),
                          'alpha': row.split('|')[2].strip(),
                          'beta': row.split('|')[3].strip(),
                          'part': row.split('|')[4].strip(),
                          'amplitude': float(row.split('|')[5]) + 0.0})
            table[-1]['occupations'] = get_rasci_occupations_list(table[-1],
                                                                  doubly_occ,
                                                                  basic_data['n_basis_functions'])

        table = sorted(table, key=operator.itemgetter('hole', 'alpha', 'beta', 'part'))
        table = sorted(table, key=lambda x: abs(x['amplitude']), reverse=True)

        # Contributions RASCI wfn
        contributions_section = section_state[enum2:]
        contributions = {'active' : float(contributions_section.split()[4]),
                         'hole': float(contributions_section.split()[6]),
                         'part': float(contributions_section.split()[8])}

        # complete dictionary
        tot_energy_units = 'au'

        state_dict = {'total_energy': tot_energy,
                      'total_energy_units': tot_energy_units,
                      'excitation_energy': exc_energy,
                      'excitation_energy_units': exc_energy_units,
                      'symmetry': state_symmetry, 
                      'multiplicity': state_multiplicity,
                      'dipole_moment': dipole_mom,
                      'transition_moment': trans_mom,
                      'dipole_moment_units': 'ua',
                      'oscillator_strength': strength,
                      'configurations': table,
                      'contributions_wfn': contributions}

        if nato_occ is not None:
            state_dict.update({'natural_occupation_numbers': nato_occ})

        if mulliken_population is not None:
            state_dict.update({'mulliken_population': mulliken_population})

        excited_states.append(state_dict)

    data_dict.update({'excited_states': excited_states})

    # Interstate transition properties
    done_interstate = bool(output.find('Interstate Transition Properties')+1)
    if done_interstate:
        ini_section = output.find('Interstate Transition Properties')
        end_section = search_bars(output, from_position=ini_section)[1]
        interstate_section = output[ini_section: end_section]
        section_sizes = search_bars(interstate_section, bar_type=r'\*' * 20)

        interstate_dict = {}
        for i, m in enumerate(re.finditer('State A: Root', interstate_section)):
            section_pair = interstate_section[m.start():m.start() + section_sizes[i]]
            end_section = search_bars(section_pair, bar_type=r'\*'*20)[0]
            section_pair = section_pair[:end_section]

            lines = section_pair.split('\n')

            state_a = int(lines[0].lower().split('root')[-1])
            state_b = int(lines[1].lower().split('root')[-1])

            pair_dict = {'state_a': state_a,
                         'state_b': state_b}

            s_a = s_b = 0
            for i, line in enumerate(lines):
                # RAS-CI SOC section
                if 'Angular momentum components' in line:
                    pair_dict['angular_momentum'] = [complex(lines[i+1+k].split()[-1].replace('i', 'j')) for k in range(3)]

                if '||gamma^AB||_total' in line:
                    pair_dict['gamma_total'] = float(lines[i+0].split()[-1])
                    pair_dict['gamma_sym'] = float(lines[i+1].split()[-1])
                    pair_dict['gamma_anti_sym'] = float(lines[i+2].split()[-1])

                if "KET: S',Sz'" in line:
                    s_a = float(lines[i].split('=')[1].split()[0])
                    s_b = float(lines[i+1].split('=')[1].split()[0])

                na = int(2 * s_a + 1)
                nb = int(2 * s_b + 1)

                if 'Spin Matrices' in line:
                    warnings.warn('Spin Matrices parsing is deprecated')
                    spinmat_x = _read_soc_matrix(lines[i + 2:], [nb, na])
                    spinmat_y = _read_soc_matrix(lines[i + 4 + nb:], [nb, na])
                    spinmat_z = _read_soc_matrix(lines[i + 6 + 2*nb:], [nb, na])
                    pair_dict['spin_matrices'] = [spinmat_x, spinmat_y, spinmat_z]

                if 'Spin matrices Sx, Sy and Sz for states' in line:
                    warnings.warn('Spin Matrices parsing is deprecated')
                    pair_dict['spin_matrices'] = [np.zeros((nb, na)).tolist()]*3

                if '1-elec SOC matrix (cm-1)' in line:
                    pair_dict['1e_soc_mat'] = _read_soc_matrix(lines[i+1:], [nb, na])

                if '1-elec SOCC' in line:
                    pair_dict['1e_socc'] = float(line.split()[3])

                if '2e-SOMF Reduced matrix elements (cm-1)' in line:
                    r, c = lines[i+1].replace('i', '').split()[-2:]
                    pair_dict['hso_l-'] = float(r) + float(c) * 1j
                    r, c = lines[i+2].replace('i', '').split()[-2:]
                    pair_dict['hso_l0'] = float(r) + float(c) * 1j
                    r, c = lines[i+3].replace('i', '').split()[-2:]
                    pair_dict['hso_l+'] = float(r) + float(c) * 1j

                if '2-elec mean-field SOC matrix (cm-1)' in line:
                    pair_dict['2e_soc_mat'] = _read_soc_matrix(lines[i + 1:], [nb, na])
                if 'Total mean-field SOC matrix (cm-1)' in line:
                    pair_dict['total_soc_mat'] = _read_soc_matrix(lines[i + 1:], [nb, na])
                if 'Mean-Field SOCC' in line:
                    pair_dict['mf_socc'] = float(line.split()[-2])
                    pair_dict['units'] = line.split()[-1]

                if 'Skipping SOCs between states' in line:
                    pair_dict['1e_soc_mat'] = np.zeros((nb, na)).tolist()
                    pair_dict['1e_socc'] = 0.0

                    pair_dict['hso_l-'] = complex(0.0)
                    pair_dict['hso_l0'] = complex(0.0)
                    pair_dict['hso_l+'] = complex(0.0)

                    pair_dict['2e_soc_mat'] = np.zeros((nb, na)).tolist()
                    pair_dict['total_soc_mat'] = np.zeros((nb, na)).tolist()

                    pair_dict['mf_socc'] = 0.0
                    pair_dict['units'] = 'cm-1'

            interstate_dict[(state_a, state_b)] = pair_dict

        _complete_interstate_pairs(interstate_dict)
        data_dict.update({'interstate_properties': interstate_dict})

    # magnetic hyperfine structure
    done_hyperfine = bool(output.find('RAS-CI: MAGNETIC HYPERFINE STRUCTURE')+1)
    if done_hyperfine:
        ini_section = output.find('RAS-CI: MAGNETIC HYPERFINE STRUCTURE')
        end_section = output.find('RAS-CI timing summary (seconds)')
        hyperfine_section = output[ini_section: end_section]

        for i_state, m in enumerate(re.finditer('RAS-CI root: ', hyperfine_section)):
            state_section = hyperfine_section[m.start():m.end()+(n_atoms*1123)]
            # end_section = search_bars(state_section, bar_type='----')[0]

            fermi_contact = []
            spin_dipole = []

            for l, i in enumerate(re.finditer('Nucleus:', state_section)):
                nucleus_section = state_section[i.start():i.end() + 1123]
                end_nucleus_section = search_bars(nucleus_section, bar_type='----')[0]
                nucleus_section = nucleus_section[:end_nucleus_section]
                lines = nucleus_section.split('\n')

                for k, line in enumerate(lines):

                    if 'fermi contact contribution' in line.lower():
                        fermi_contact.append(float(lines[k+1].split()[0]))

                    if 'spin dipole contribution' in line.lower():
                        dipole_x = float(lines[k+1].split()[0])
                        dipole_y = float(lines[k+1].split()[1])
                        dipole_z = float(lines[k+1].split()[2])

                        spin_dipole.append([dipole_x, dipole_y, dipole_z])

            data_dict['excited_states'][i_state].update({'hyperfine_couplings': {'fermi_contact': fermi_contact,
                                                                                 'spin_dipole': spin_dipole}})

    return data_dict

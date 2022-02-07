__author__ = 'Abel Carreras'

import re
import operator
from pyqchem.structure import Structure
from pyqchem.parsers.common import read_basic_info, get_rasci_occupations_list, search_bars, standardize_vector


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
    # for line in lines:
    #     print(line)

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

    # basic info
    enum = output.find('Molecular Point Group')
    basic_data = read_basic_info(output[enum:enum + 5000])

    # scf_energy
    enum = output.find('SCF   energy in the final basis set')
    scf_energy = float(output[enum:enum+100].split()[8])

    data_dict['scf_energy'] = scf_energy
    # total energy
    # enum = output.find('Total energy in the final basis set')
    # total_energy = float(output[enum:enum+100].split()[8])

    # RASCI dimensions
    ini_section = output.find('RAS-CI Dimensions')
    end_section = search_bars(output, from_position=enum, bar_type='\*\*\*')[1]
    dimension_section = output[ini_section: end_section]

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

    rasci_dimensions = {'doubly_occupied': doubly_occ,
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
            section_mulliken = output[m.end() + enum: m.end() + 10000 + enum]  # 10000: assumed to max of section
            section_mulliken = section_mulliken[:section_mulliken.find('Natural Orbitals stored in FCHK')]
            section_attachment = section_mulliken.split('\n')[9+n_atoms:9+n_atoms*2]

            mulliken_adiabatic.append({'attach': [float(l.split()[1]) for l in section_attachment],
                                       'detach': [float(l.split()[2]) for l in section_attachment],
                                       'total': [float(l.split()[3]) for l in section_attachment]})

        mulliken_diabatic = []
        enum = output.find('showing H in diabatic representation')
        for m in re.finditer('Mulliken Analysis of Diabatic State', output[enum:]):
            section_mulliken = output[m.end() + enum: m.end() + 10000 + enum]  # 10000: assumed to max of section
            section_mulliken = section_mulliken[:section_mulliken.find('Natural Orbitals stored in FCHK')]
            section_attachment = section_mulliken.split('\n')[9+n_atoms:9+n_atoms*2]

            mulliken_diabatic.append({'attach': [float(l.split()[1]) for l in section_attachment],
                                      'detach': [float(l.split()[2]) for l in section_attachment],
                                      'total': [float(l.split()[3]) for l in section_attachment]})

        enum = output.find('Transition dipole moment - diabatic states')

        tdm_section = output[enum: enum + 70 * len(rot_matrix)]

        diabatic_tdm = []
        for m in re.finditer('TDM', tdm_section):
            diabatic_tdm.append([float(n) for n in tdm_section[m.end(): m.end()+70][14:].split()[:3]])

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
        # print('ll found', m.start(), m.end())

        section_state = output[m.end():m.end() + 10000]  # 10000: assumed to max of section
        section_state = section_state[:section_state.find('********')]

        enum = section_state.find('RAS-CI total energy for state')
        section_state = section_state[:enum]

        # energies
        tot_energy = float(section_state.split()[1])
        exc_energy_units = section_state.split()[4][1:-1]
        exc_energy = float(section_state.split()[6])
        state_multiplicity = section_state.split()[8] if section_state.split()[8] != ':' else section_state.split()[9]

        # dipole moment
        enum = section_state.find('Dipole Moment')
        dipole_mom = [float(section_state[enum:].split()[2]) + 0.0,
                      float(section_state[enum:].split()[4]) + 0.0,
                      float(section_state[enum:].split()[6]) + 0.0]

        # Transition moment
        enum = section_state.find('Trans. Moment')
        if enum > -1:
            trans_mom = [float(section_state[enum:].split()[2]) + 0.0,
                         float(section_state[enum:].split()[4]) + 0.0,
                         float(section_state[enum:].split()[6]) + 0.0]
            trans_mom = standardize_vector(trans_mom)
            strength = float(section_state[enum:].split()[10])
        else:
            trans_mom = None
            strength = None

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
                                                                  data_dict['structure'],
                                                                  basic_data['n_basis_functions'])

        table = sorted(table, key=operator.itemgetter('hole', 'alpha', 'beta', 'part'))

        # Contributions RASCI wfn
        contributions_section = section_state[enum2:]
        contributions = {'active' : float(contributions_section.split()[4]),
                         'hole': float(contributions_section.split()[6]),
                         'part': float(contributions_section.split()[8])}

        # complete dictionary
        tot_energy_units = 'au'
        excited_states.append({'total_energy': tot_energy,
                               'total_energy_units': tot_energy_units,
                               'excitation_energy': exc_energy,
                               'excitation_energy_units': exc_energy_units,
                               'multiplicity': state_multiplicity,
                               'dipole_moment': dipole_mom,
                               'transition_moment': trans_mom,
                               'dipole_moment_units': 'ua',
                               'oscillator_strength': strength,
                               'configurations': table,
                               'contributions_fwn': contributions})

    data_dict.update({'excited_states': excited_states})

    # Interstate transition properties
    done_interstate = bool(output.find('Interstate Transition Properties')+1)
    if done_interstate:
        ini_section = output.find('Interstate Transition Properties')
        end_section = search_bars(output, from_position=ini_section)[1]
        interstate_section = output[ini_section: end_section]

        interstate_dict = {}
        for m in re.finditer('State A: Root', interstate_section):
            section_pair = interstate_section[m.start():m.start() + 10000]
            section_pair = section_pair[:section_pair.find('********')]

            lines = section_pair.split('\n')

            state_a = int(lines[0].split()[-1])
            state_b = int(lines[1].split()[-1])

            pair_dict = {'state_a': state_a,
                         'state_b': state_b}

            s_a = s_b = 0
            for i, line in enumerate(lines):
                # RAS-CI SOC section
                if '||gamma^AB||_total' in line:
                    pair_dict['gamma_total'] = float(lines[i+0].split()[-1])
                    pair_dict['gamma_sym'] = float(lines[i+1].split()[-1])
                    pair_dict['gamma_anti_sym'] = float(lines[i+2].split()[-1])

                if "KET: S',Sz'" in line:
                    s_a = float(lines[i].split('=')[1].split()[0])
                    s_b = float(lines[i+1].split('=')[1].split()[0])
                if '1-elec SOC matrix (cm-1)' in line:
                    pair_dict['1e_soc_mat'] = _read_soc_matrix(lines[i+1:], [int(2*s_b + 1), int(2*s_a+1)])
                if '2e-SOMF Reduced matrix elements (cm-1)' in line:
                    r, c = lines[i+1].split()[-2:]
                    pair_dict['hso_l-'] = float(r) + float(c) * 1j
                    r, c = lines[i+2].split()[-2:]
                    pair_dict['hso_l0'] = float(r) + float(c) * 1j
                    r, c = lines[i+3].split()[-2:]
                    pair_dict['hso_l+'] = float(r) + float(c) * 1j

                if '2-elec mean-field SOC matrix (cm-1)' in line:
                    pair_dict['2e_soc_mat'] = _read_soc_matrix(lines[i + 1:], [int(2 * s_b + 1), int(2 * s_a + 1)])
                if 'Total mean-field SOC matrix (cm-1)' in line:
                    pair_dict['total_soc_mat'] = _read_soc_matrix(lines[i + 1:], [int(2 * s_b + 1), int(2 * s_a + 1)])
                if 'Mean-Field SOCC' in line:
                    pair_dict['mf_socc'] = float(line.split()[-2])
                    pair_dict['units'] = line.split()[-1]

            interstate_dict[(state_a, state_b)] = pair_dict
        data_dict.update({'interstate_properties': interstate_dict})

    return data_dict

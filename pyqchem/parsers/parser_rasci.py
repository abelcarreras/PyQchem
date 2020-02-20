__author__ = 'Abel Carreras'

import re
import operator
from pyqchem.utils import standardize_vector
from pyqchem.structure import Structure
from pyqchem.utils import search_bars

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
    matrix = []
    for ib in range(dimensions[0]):
        real = []
        complex = []
        for j in range((dimensions[1] + 1) // 2):
            real += lines[j*3 + 2 + ib].split()[1:][0::2]
            complex += lines[j*3 +2+ib].split()[1:][1::2]

        row = [float(r) + float(c) * 1j for r, c in zip(real, complex)]
        matrix.append(row)

    return matrix


def rasci(output):
    """
    Parser for RAS-CI calculations
    Include:
    - Diabatization scheme data
    - Structure
    - Adiabatic states

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
                                atomic_elements=symbols,
                                charge=charge,
                                multiplicity=multiplicity)

    enum = output.find('Standard Nuclear Orientation')
    section_structure = output[enum:enum + 100*structure_input.get_number_of_atoms()].split('\n')
    section_structure = section_structure[3:structure_input.get_number_of_atoms()+3]
    coordinates = [ [float(num) for num in s.split()[2:]] for s in section_structure]

    data_dict['structure'] = Structure(coordinates=coordinates,
                                       atomic_elements=symbols,
                                       charge=charge,
                                       multiplicity=multiplicity)

    # scf energy
    enum = output.find('SCF   energy in the final basis set')
    scf_energy = float(output[enum:enum+100].split()[8])

    data_dict['scf energy'] = scf_energy
    # total energy
    # enum = output.find('Total energy in the final basis set')
    # total_energy = float(output[enum:enum+100].split()[8])

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

        data_dict['diabatization'] = {'rot_matrix': rot_matrix,
                                      'adiabatic_matrix': adiabatic_matrix,
                                      'diabatic_matrix': diabatic_matrix,
                                      'mulliken_analysis': mulliken_diabatic,
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
        mul = section_state.split()[8]

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
        table = sorted(table, key=operator.itemgetter('hole', 'alpha', 'beta', 'part'))

        # Contributions RASCI wfn
        contributions_section = section_state[enum2:]
        contributions = {'active' : float(contributions_section.split()[4]),
                         'hole': float(contributions_section.split()[6]),
                         'part': float(contributions_section.split()[8])}

        # complete dictionary
        tot_energy_units = 'au'
        excited_states.append({'total_energy': tot_energy,
                               'total energy units': tot_energy_units,
                               'excitation_energy': exc_energy,
                               'excitation energy units': exc_energy_units,
                               'multiplicity': mul,
                               'dipole_moment': dipole_mom,
                               'transition_moment': trans_mom,
                               'oscillator_strength': strength,
                               'configurations': table,
                               'contributions_fwn': contributions})

    data_dict.update({'excited states rasci': excited_states})

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
            # print(section_pair)
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

                if 'KET: S,Sz' in line:
                    s_a = float(lines[i].split()[-2])
                    s_b = float(lines[i+1].split()[-2])
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
                if 'Total mean-field SOC matrix (cm-1' in line:
                    pair_dict['total_soc_mat'] = _read_soc_matrix(lines[i + 1:], [int(2 * s_b + 1), int(2 * s_a + 1)])
                if 'Mean-Field SOCC' in line:
                    pair_dict['mf_socc'] = float(line.split()[-2])
                    pair_dict['units'] = line.split()[-1]

            interstate_dict[(state_a, state_b)] = pair_dict
        data_dict.update({'interstate_properties': interstate_dict})

    return data_dict

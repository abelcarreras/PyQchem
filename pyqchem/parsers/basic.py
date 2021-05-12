__author__ = 'Abel Carreras'
from pyqchem.parsers.common import search_bars


def _get_orbital_energies(orbitals_section):
    # print(orbitals_section)
    occupied = orbitals_section.find('-- Occupied --')
    virtual = orbitals_section.find('-- Virtual --')

    occupied_section = orbitals_section[occupied:virtual]
    virtual_section = orbitals_section[virtual:]

    occupied_energies = [float(energy) if energy != '********' else None for energy in ' '.join(occupied_section.split('\n')[1::2]).split()]
    virtual_energies = [float(energy) if energy != '********' else None for energy in ' '.join(virtual_section.split('\n')[1::2]).split()]

    return occupied_energies + virtual_energies


def basic_parser_qchem(output):
    """
    This showcases the format of  Q-Chem version parser compatibility.
    Just by creating a docstring with the following line:

    compatibility: 5.1, 5.2+

    will activate the check for parser-qchem version compatibility.
    If they are not compatible a warning will rise.

    """
    data_dict = {}

    # scf_energy
    enum = output.find('Total energy in the final basis set')
    data_dict['scf_energy'] = float(output[enum:enum+100].split()[8])

    # Orbitals energy
    enum = output.find('Orbital Energies (a.u.)')
    bars = search_bars(output, from_position=enum, bar_type='----')
    orbitals_section = output[enum:bars[1]]

    alpha_mos = orbitals_section.find('Alpha MOs')
    beta_mos = orbitals_section.find('Beta MOs')

    # print(orbitals_section[alpha_mos:beta_mos])
    # print(orbitals_section[alpha_mos:])
    if beta_mos > 0:
        alpha_energies = _get_orbital_energies(orbitals_section[alpha_mos:beta_mos])
        beta_energies = _get_orbital_energies(orbitals_section[beta_mos:])
    else:
        alpha_energies = _get_orbital_energies(orbitals_section[alpha_mos:beta_mos])
        beta_energies = alpha_energies

    data_dict['orbital_energies'] = {'alpha': alpha_energies, 'beta': beta_energies, 'units': 'au'}

    # Mulliken Net Atomic Charges
    enum = output.find('Ground-State Mulliken Net Atomic Charges')
    bars = search_bars(output, from_position=enum, bar_type='----')
    mulliken_section = output[bars[0]:bars[1]]
    data_dict['mulliken_charges'] = [float(line.split()[2]) for line in mulliken_section.split('\n')[1:-1]]

    # Multipole Moments
    enum = output.find('Cartesian Multipole Moments')
    bars = search_bars(output, from_position=enum, bar_type='----')
    multipole_section = output[bars[0]: bars[1]]
    multipole_lines =  multipole_section.split('\n')[1:-1]

    multipole_dict = {}

    multipole_dict['charge'] = float(multipole_lines[1])
    multipole_dict['charge_units'] = 'ESU x 10^10'

    multipole_dict['dipole_moment'] = [float(val) for val in multipole_lines[3].split()[1::2]]
    multipole_dict['dipole_units'] = 'Debye'

    quadrupole = [float(val) for val in multipole_lines[6].split()[1::2]] + \
                 [float(val) for val in multipole_lines[7].split()[1::2]]

    # create quadrupole array
    multipole_dict['quadrupole_moment'] = [[quadrupole[0], quadrupole[1], quadrupole[2]],
                                           [quadrupole[1], quadrupole[3], quadrupole[4]],
                                           [quadrupole[2], quadrupole[4], quadrupole[5]]]

    # multipole_dict['quadrupole_moment'] = [float(val) for val in multipole_lines[6].split()[1::2]] + \
    #                                       [float(val) for val in multipole_lines[7].split()[1::2]]


    multipole_dict['quadrupole_units'] = 'Debye-Ang'

    octopole = [float(val) for val in multipole_lines[9].split()[1::2]] + \
               [float(val) for val in multipole_lines[10].split()[1::2]] + \
               [float(val) for val in multipole_lines[11].split()[1::2]] + \
               [float(val) for val in multipole_lines[12].split()[1::2]]

    # create octopole array
    multipole_dict['octopole_moment'] = [
        [[octopole[0], octopole[1], octopole[4]],
         [octopole[1], octopole[2], octopole[5]],
         [octopole[4], octopole[5], octopole[7]]],

        [[octopole[1], octopole[2], octopole[5]],
         [octopole[2], octopole[3], octopole[6]],
         [octopole[5], octopole[6], octopole[8]]],

        [[octopole[4], octopole[5], octopole[7]],
         [octopole[5], octopole[6], octopole[8]],
         [octopole[7], octopole[8], octopole[9]]],
    ]

    # multipole_dict['octopole_moment'] = [float(val) for val in multipole_lines[9].split()[1::2]] + \
    #                                     [float(val) for val in multipole_lines[10].split()[1::2]] + \
    #                                     [float(val) for val in multipole_lines[11].split()[1::2]] + \
    #                                     [float(val) for val in multipole_lines[12].split()[1::2]]


    multipole_dict['octopole_units'] = 'Debye-Ang^2'

    data_dict['multipole'] = multipole_dict

    return data_dict



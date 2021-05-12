__author__ = 'Abel Carreras'

import re
import operator
from pyqchem.parsers.common import standardize_vector


def basic_rasci(output):
    """
    Basic Parser for RAS-CI calculations

    :param output:
    :return:
    """

    # scf_energy
    enum = output.find('SCF   energy in the final basis set')
    scf_energy = float(output[enum:enum+100].split()[8])

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

        # amplitudes table
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
                               'total_energy_units': tot_energy_units,
                               'excitation_energy': exc_energy,
                               'excitation_energy_units': exc_energy_units,
                               'multiplicity': mul,
                               'dipole_moment': dipole_mom,
                               'transition_moment': trans_mom,
                               # 'oscillator_strength': strength,
                               'configurations': table,
                               'contributions_fwn': contributions})

    return {'scf_energy': scf_energy,
            'excited_states': excited_states}

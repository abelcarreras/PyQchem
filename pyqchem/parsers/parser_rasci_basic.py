__author__ = 'Abel Carreras'

import re


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


def basic_rasci(output):

    # scf energy
    enum = output.find('SCF   energy in the final basis set')
    scf_energy = float(output[enum:enum+100].split()[8])

    # total energy
    # enum = output.find('Total energy in the final basis set')
    # total_energy = float(output[enum:enum+100].split()[8])

    # excited states data
    excited_states = []
    for m in re.finditer('RAS-CI total energy for state', output):
        # print('ll found', m.start(), m.end())

        section_state = output[m.end():m.end() + 10000]
        enum = section_state.find('RAS-CI total energy for state')
        section_state = section_state[:enum]

        # energies
        tot_energy = float(section_state.split()[1])
        exc_energy = float(section_state.split()[6])
        mul = section_state.split()[8]

        # dipole moment
        enum = section_state.find('Dipole Moment')
        dipole_mom = [float(section_state[enum:enum+100].split()[2]) + 0.0,
                      float(section_state[enum:enum + 100].split()[4]) + 0.0,
                      float(section_state[enum:enum + 100].split()[6]) + 0.0]

        # Transition moment
        enum = section_state.find('Trans. Moment')
        if enum > -1:
            trans_mom = [float(section_state[enum:enum + 100].split()[2]) + 0.0,
                         float(section_state[enum:enum + 100].split()[4]) + 0.0,
                         float(section_state[enum:enum + 100].split()[6]) + 0.0]
            trans_mom = standardize_vector(trans_mom)
        else:
            trans_mom = None

        # amplitudes table
        enum = section_state.find('AMPLITUDE')
        enum2 = section_state.find('Contributions')
        section_table = section_state[enum:enum2].split('\n')[2:-2]

        # ' HOLE  | ALPHA | BETA  | PART | AMPLITUDE'

        table = []
        for row in section_table:
            table.append({'hole': row.split('|')[1].strip(),
                          'alpha': row.split('|')[2].strip(),
                          'beta': row.split('|')[3].strip(),
                          'part': row.split('|')[4].strip(),
                          'amplitude': float(row.split('|')[5]) + 0.0})

        # Contributions RASCI wfn
        enum3 = section_state[enum2:enum2+300].find('********') + enum2
        contributions_section = section_state[enum2: enum3]
        contributions = {'active' : float(contributions_section.split()[4]),
                         'hole': float(contributions_section.split()[6]),
                         'part': float(contributions_section.split()[8])}

        # complete dictionary
        excited_states.append({'total_energy': tot_energy,
                               'excitation_energy': exc_energy,
                               'multiplicity': mul,
                               'dipole_moment': dipole_mom,
                               'transition_moment': trans_mom,
                               'amplitudes': table,
                               'contributions_fwn': contributions})

    return {'scf energy': scf_energy,
            'excited states rasci': excited_states}

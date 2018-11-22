__author__ = 'Abel Carreras'

import re


def search_bars(output, from_position=0):
    output = output[from_position:]
    positions = []
    previous = 0
    for m in re.finditer('---', output):
        if m.start() > previous + 1:
            positions.append(m.start() + from_position)
        previous = m.end()

    return positions


def basic_parser_qchem(output):

    # scf energy
    enum = output.find('Total energy in the final basis set')
    scf_energy = float(output[enum:enum+100].split()[8])

    # excited states data
    excited_states = []
    for m in re.finditer('RAS-CI total energy for state', output):
        # print('ll found', m.start(), m.end())
        tot_energy = float(output[m.end():m.end() + 100].split()[1])
        exc_energy = float(output[m.end():m.end() + 100].split()[6])
        mul = output[m.end():m.end() + 100].split()[8]

        excited_states.append({'total energy': tot_energy,
                               'excitation energy': exc_energy,
                               'multiplicity': mul})

    # CIS Excitation energies
    enum = output.find('CIS Excitation Energies')
    excited_states_cis = []
    if enum > 0:
        bars = search_bars(output, from_position=enum)

        output_cis = output[bars[0]:bars[1]]

        for m in re.finditer('Excited state ', output_cis):
            state_cis_words = output_cis[m.end():].split()

            exc_energy = float(state_cis_words[5])
            exc_energy_units = state_cis_words[3][1:-1]
            tot_energy = float(state_cis_words[11])
            try:
                tot_energy_units = state_cis_words[12]
                mul = state_cis_words[14]
                trans_mom = [float(mom) for mom in [state_cis_words[17],
                                                    state_cis_words[19],
                                                    state_cis_words[21]]]
                strength = float(state_cis_words[25])
            except ValueError:
                # old version of qchem (< 5.01)
                tot_energy_units = 'au'
                mul = state_cis_words[13]
                trans_mom = [float(mom) for mom in [state_cis_words[16],
                                                    state_cis_words[18],
                                                    state_cis_words[20]]]
                strength = float(state_cis_words[24])

            transitions = []
            for i, word in enumerate(state_cis_words):
                if word == '-->':
                    origin = int(state_cis_words[i-1][:-1])
                    target = int(state_cis_words[i+2][:-1])
                    amplitude = float(state_cis_words[i+5])

                    transitions.append({'origin': origin,
                                        'target': target,
                                        'amplitude': amplitude})
                if word == 'Excited':
                    break

            excited_states_cis.append({'total energy': tot_energy,
                                       'total energy units': tot_energy_units,
                                       'excitation energy': exc_energy,
                                       'excitation energy units': exc_energy_units,
                                       'multiplicity': mul,
                                       'transition moment': trans_mom,
                                       'strength': strength,
                                       'transitions': transitions})

    return {'scf energy': scf_energy,
            'excited states': excited_states,
            'excited states cis': excited_states_cis}


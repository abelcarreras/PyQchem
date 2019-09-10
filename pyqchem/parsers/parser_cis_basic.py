__author__ = 'Abel Carreras'

import re
from pyqchem.utils import search_bars, standardize_vector


def basic_cis(output):
    """
    Parser for CIS calculations

    :param output:
    :return:
    """

    # scf energy
    enum = output.find('Total energy in the final basis set')
    scf_energy = float(output[enum:enum+100].split()[8])

    # CIS excited states
    # enum = output.find('CIS Excitation Energies')
    enum = list(re.finditer('CIS Excitation Energies', output))[-1].end()
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
                                       'transition moment': standardize_vector(trans_mom),
                                       'strength': strength,
                                       'transitions': transitions})

    return {'scf energy': scf_energy,
            'excited states cis': excited_states_cis}


__author__ = 'Abel Carreras'

import re
from pyqchem.utils import search_bars, standardize_vector
from pyqchem.errors import ParserError


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

    # scf_energy
    enum = output.find('Total energy in the final basis set')
    data_dict['scf_energy'] = float(output[enum:enum+100].split()[8])

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

            excited_states.append({'total_energy': tot_energy,
                                   'total_energy_units': tot_energy_units,
                                   'excitation_energy': exc_energy,
                                   'excitation_energy_units': exc_energy_units,
                                   'multiplicity': mul,
                                   'transition_moment': standardize_vector(trans_mom),
                                   'strength': strength,
                                   'transitions': transitions})

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
                    raise ParserError('basic_cis', 'State multiplicity error')
            return labels, nt-1, ns-1

        labels, n_triplet, n_singlet = label_states(excited_states)

        for i, label in enumerate(labels):
            data_interstate[(i+1, 0)] = {'soc': [0j, 0j, 0j], 'soc_units': 'cm-1'}
            data_interstate[(0, i+1)] = {'soc': [0j, 0j, 0j], 'soc_units': 'cm-1'}
            for j, label2 in enumerate(labels):
                if (label[0] == 'S' or label2[0] == 'S') and (label[0] != label2[0]):
                    data_interstate[(i+1, j+1)] = {'soc': [0j, 0j, 0j], 'soc_units': 'cm-1'}
                elif label[0] == 'T' and label2[0] == 'T':
                    data_interstate[(i+1, j+1)] = {'soc': [[0j, 0j, 0j], [0j, 0j, 0j], [0j, 0j, 0j]], 'soc_units': 'cm-1'}
                elif label[0] == 'S' and label2[0] == 'S':
                    data_interstate[(i+1, j+1)] = {'soc': [0j], 'soc_units': 'cm-1'}
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
                                    data_interstate[(i+1, j+1)]['soc'][k2][k] = _list_to_complex(line.split()[1:4])
                                    data_interstate[(i+1, j+1)]['soc'][k][k2] = _list_to_complex(line.split()[1:4])
                                    data_interstate[(j+1, i+1)]['soc'][k2][k] = _list_to_complex(line.split()[1:4])
                                    data_interstate[(j+1, i+1)]['soc'][k][k2] = _list_to_complex(line.split()[1:4])
                                    break

                    elif label[0] == 'S':
                        for k, ms in enumerate([-1, 0, 1]):
                            enum = soc_section.find('SOC between the {} state and excited triplet states (ms={})'.format(label, ms))
                            for line in soc_section[enum:enum+50*(n_triplet+1)].split('\n'):
                                if len(line.split()) == 0:
                                    break
                                if line.split()[0] == '{}(ms={})'.format(label2, ms):
                                    data_interstate[(i+1, j+1)]['soc'][k] = _list_to_complex(line.split()[1:4])
                                    data_interstate[(j+1, i+1)]['soc'][k] = _list_to_complex(line.split()[1:4])
                                    break
                    else:
                        raise ParserError('basic_cis', 'SOC reading error')

                enum = soc_section.find('SOC between the singlet ground state and excited triplet states (ms={})'.format(ms2))
                for line in soc_section[enum:enum+50*(n_triplet+1)].split('\n'):
                    if len(line.split()) == 0:
                        break
                    if line.split()[0] == '{}(ms={})'.format(label, ms2):
                        data_interstate[(i+1, 0)]['soc'][k2] = _list_to_complex(line.split()[1:4])
                        data_interstate[(0, i+1)]['soc'][k2] = _list_to_complex(line.split()[1:4])
                        break

        data_dict['interstate_properties'] = data_interstate

    return data_dict


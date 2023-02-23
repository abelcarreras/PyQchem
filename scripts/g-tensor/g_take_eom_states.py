import numpy as np

from g_plots import get_bar_chart

def get_eom_type(eom_input):
    """
    Define the type of equation of motions eom_input
    :param: eom_input
    :return: eom_type
    """
    searches = ['Solving for']

    with open(eom_input) as data:
        for line in data:
            if any(i in line for i in searches):

                if ('EOMIP' in line): eom_type = 'eomIP'
                if ('EOMEA' in line): eom_type = 'eomEA'
    return eom_type


def get_SCF_energy(eom_input):
    """
    Define the type of equation of motions eom_input
    :param: eom_input
    :return: SCF_energy
    """
    searches = ['SCF   energy in the final basis set']

    with open(eom_input, encoding="utf8") as data:
        for line in data:
            if any(i in line for i in searches):
                element = line[39:]
                SCF_energy = float(element)
    return SCF_energy


def get_irreps_energies(eom_input):
    """
    Take the irreducible representations and the transitions in each irrep.,
    as well as their total and excitation energies
    :param: eom_input
    :return:  irreps, energies, excit_energies
    """
    searches = [' EOMIP transition', 'EOMEA transition']
    irreps = []
    energy_list = []
    excitation_energy_list_eV = []

    with open(eom_input) as data:
        for line in data:
            if any(i in line for i in searches):
                elements = []

                line = line.replace('EOMIP', '')
                line = line.replace('EOMEA', '')
                line = line.replace('transition', '')
                line = line.replace('/', ' ')

                # Take eom transition number and irreps
                element = line
                irreps.append(element.split())

                # Take total and excitation energies
                next_line = next(data)
                elements.append(next_line.split())

                element = elements[0][3]
                energy_list.append(element)
                element = elements[0][8]
                excitation_energy_list_eV.append(element)

    total_energies = np.array(energy_list, dtype=float)
    excitation_energies_respect_reference_eV = np.array(excitation_energy_list_eV, dtype=float)

    # Excitation energies are expressed as difference with minimum energy (usually ground state)
    excitation_energies = excitation_energies_respect_reference_eV
    excitation_energies[:] = (excitation_energies_respect_reference_eV[:] - min(
        excitation_energies_respect_reference_eV)) \
                             / 27.211399  # always in a.u.

    search_2 = [' State A: ']
    elements = []

    with open(eom_input) as data:
        for line in data:
            if any(i in line for i in search_2):
                line = line.replace('/', ' ')
                elements.append(line.split())
                element = elements[0][4:]
                break

    optimized_state_index = irreps.index(element)

    return irreps, optimized_state_index, total_energies, excitation_energies


def get_eom_SOCC_values(eom_input, optimized_state_index):
    """
    Take the contribution of hole configurations for each state
    :param: lines_eom_input
    :return: hole_contributions
    """
    search_1 = ' Mean-field SO (cm-1)'
    search_2 = 'SOCC '
    elements = []
    SOCC_list = []

    i = 0

    with open(eom_input, encoding="utf8") as file:
        for line in file:

            if search_1 in line:

                next_line = next(file)

                if search_2 in next_line:
                    if (i == optimized_state_index): SOCC_list.append('0.000000')
                    # SOCC in main state, in which interstate properties are calculated, is zero

                    elements = []
                    elements.append(next_line.split())
                    SOCC_list.append(elements[0][2])

                    i += 1

    return SOCC_list


def get_maximum_amplitude_orbitals(eom_input, eom_type):
    """
     Gives orbitals related with maximum amplitudes transitions
     :param: eom_input
     :return: significant_orbitals
     """

    def get_transitions_between_orbitals(line, eom_type):
        """
         Used in "get_significant_orbitals". It takes the
         "Transitions between orbitals" lines without the amplitudes
         :param: line
         :return: transition_orbitals
         """
        transition_orbitals = []

        if (eom_type == 'eomIP'):
            line_split = []

            if ('infty' in line):
                line_split.append(line.split())
                element = line_split[0][1] + ' ' + line_split[0][2]
                transition_orbitals.append(element)
                element = '-'
                transition_orbitals.append(element)
                element = '-'
                transition_orbitals.append(element)

                # element = line[13:22]
                # transition_orbitals.append(element)
                # element = line[28:36]
                # transition_orbitals.append(element)
                # element = line[46:55]
                # transition_orbitals.append(element)
                # print(transition_orbitals)

            elif ('infty' not in line):
                line_split.append(line.split())
                element = line_split[0][1] + ' ' + line_split[0][2]
                transition_orbitals.append(element)
                element = line_split[0][4] + ' ' + line_split[0][5]
                transition_orbitals.append(element)
                element = line_split[0][8] + ' ' + line_split[0][9]
                transition_orbitals.append(element)

        elif (eom_type == 'eomEA'):
            line_split = []

            if ('infty' in line):
                line_split.append(line.split())
                element = '-'
                transition_orbitals.append(element)
                element = '-'
                transition_orbitals.append(element)
                element = line_split[0][3] + ' ' + line_split[0][4]
                transition_orbitals.append(element)

                # transition_orbitals = []
                # element = line[13:22]
                # transition_orbitals.append(element)
                # element = line[32:40]
                # transition_orbitals.append(element)
                # element = line[46:54]
                # transition_orbitals.append(element)

            elif ('infty' not in line):
                line_split.append(line.split())
                element = line_split[0][1] + ' ' + line_split[0][2]
                transition_orbitals.append(element)
                element = line_split[0][5] + ' ' + line_split[0][6]
                transition_orbitals.append(element)
                element = line_split[0][8] + ' ' + line_split[0][9]
                transition_orbitals.append(element)

        return transition_orbitals

    def get_summary_significant_orbitals(max_amplitude_transition_line, transition_orbitals, next_line, data, orbitals):
        """
         Used in "get_significant_orbitals". It takes the orbital number from
         "Summary of significant orbitals" section, that has the symmetry of those irreps
         in "Transitions between orbitals"
         :param: max_amplitude_transition_line, transition_orbitals, next_line, data, orbitals
         :return: orbitals
         """
        for i in range(0, len(transition_orbitals)):
            if ('-' in transition_orbitals[i]):
                orbitals.append('-')

            elif ('-' not in transition_orbitals[i]):
                while transition_orbitals[i] not in next_line: next_line = next(data)

                line_orbitals = []
                line_orbitals.append(next_line.split())
                element = line_orbitals[0][0]
                orbitals.append(element)
        return orbitals

    searches = ['Transitions between orbitals']
    significant_orbitals = []

    with open(eom_input) as data:
        for line in data:
            if any(i in line for i in searches):

                amplitude_list = []
                all_transition_lines = []
                next_line = next(data)

                while '->' in next_line:
                    # Take amplitudes without their sign
                    next_line = next_line.replace('-', '')

                    line_split = []
                    line_split.append(next_line.split())
                    element = line_split[0][0]
                    amplitude_list.append(element)

                    all_transition_lines.append(next_line)
                    next_line = next(data)

                # Obtain line with the highest amplitude in all transitions
                max_index = amplitude_list.index(max(amplitude_list))
                max_amplitude_transition_line = all_transition_lines[max_index]

                transition_orbitals = get_transitions_between_orbitals(max_amplitude_transition_line, eom_type)
                significant_orbitals = get_summary_significant_orbitals(max_amplitude_transition_line,transition_orbitals,next_line, data, significant_orbitals)

    return significant_orbitals


def prepare_presentation_list(eom_version):
    presentation_list = []

    if (eom_version == 'eomIP'):
        presentation_list.append(['Symmetry', 'Transition', 'Energy (au)', 'Excitation energy (eV)', 'SOCC (cm-1)','Occ. 1', 'Occ. 2', 'Virtual 1'])
    elif (eom_version == 'eomEA'):
        presentation_list.append(['Symmetry', 'Transition', 'Energy (au)', 'Excitation energy (eV)', 'SOCC (cm-1)','Occ. 1', 'Virtual 1', 'Virtual 2'])
    return presentation_list


def eom_results(eom_presentation_list, irreps, total_energies, excitation_energies, eom_soc_constants, orbitals):
    def second_smallest_number(numbers):
        m1 = m2 = float('inf')
        for x in numbers:
            if x <= m1:
                m1, m2 = x, m1
            elif x < m2:
                m2 = x
        return m2

    eom_state = 0
    transition = 0
    selected_eom_excitation_energies_list = []
    selected_eom_soc_constants_list = []

    threshold_excitation_energy = 8  # in eV, energy difference with respect to ground state

    for i in range(0, len(irreps)):
        state_irreps = irreps[i][1]
        excit_energy = np.round(float(total_energies[i]), 3)
        excitation_energy = np.round(float(excitation_energies[i] * 27.211399), 3)
        SOCC = np.round(float(eom_soc_constants[i]), 2)
        # SOCC = 0

        if (excitation_energy <= threshold_excitation_energy):
            selected_eom_excitation_energies_list.append(excitation_energies[i])
            # selected_eom_soc_constants_list.append(eom_soc_constants[i])
            selected_eom_soc_constants_list.append(SOCC)

            transition += 1

            eom_presentation_list.append([state_irreps, transition, excit_energy, excitation_energy, SOCC,orbitals[i * 3], orbitals[i * 3 + 1], orbitals[i * 3 + 2]])
            eom_state += 1

        if (i < len(irreps) - 1) and (irreps[i][1] != irreps[i + 1][1]):
            eom_presentation_list.append(['---'])

    second_smallest_excit_energy = second_smallest_number(excitation_energies)
    if (second_smallest_excit_energy > threshold_excitation_energy):
        print('Increase the threshold_excitation_energy')
        exit()

    selected_excitation_energies_eom = np.array(selected_eom_excitation_energies_list, dtype=float)
    selected_eom_soc_constants = np.array(selected_eom_soc_constants_list, dtype=float)

    return eom_presentation_list, selected_excitation_energies_eom, selected_eom_soc_constants


def get_eom_transitions_analysis(eom_input):
    """"
    Obtain transitions analysis in Equation of motion output from Q-Chem.
    :param: eom_input
    :return: eom_presentation_list
    """
    from g_take_eom_states import get_eom_type, get_SCF_energy, get_irreps_energies, prepare_presentation_list, \
        get_maximum_amplitude_orbitals, get_eom_SOCC_values, eom_results

    eom_version = get_eom_type(eom_input)

    SCF_reference_energy = get_SCF_energy(eom_input)

    irreps, optimized_state_index, total_energies_eom, all_excitation_energies_eom = get_irreps_energies(eom_input)

    eom_SOCC = get_eom_SOCC_values(eom_input, optimized_state_index)

    orbitals = get_maximum_amplitude_orbitals(eom_input, eom_version)

    eom_presentation_list = prepare_presentation_list(eom_version)

    eom_presentation_list, selected_excitation_energies_eom, selected_eom_soc_constants = eom_results(
        eom_presentation_list, irreps, total_energies_eom, all_excitation_energies_eom, eom_SOCC, orbitals)

    print("")
    print("------------------------")
    print("    eom-CC ANALYSIS ")
    print("------------------------")
    print("SCF   reference energy: ", SCF_reference_energy)

    print('\n'.join([''.join(['{:^20}'.format(item) for item in row]) \
                     for row in (eom_presentation_list)]))

"""
    EXCITED STATES MOST IMPORTANT CONFIGURATIONS
 Analysis of the excited states and print of those
 orbitals involved in the configurations with the
 highest amplitudes in the excited states
"""
import numpy as np
import sys

from g_read import get_number_of_states, get_eigenenergies, get_selected_states, \
    get_symmetry_states, get_hole_part_contributions, get_SOCC_values, get_ground_state_orbital_momentum

def get_RAS2_and_RAS_OCC_and_HOMO(RAS_input):
    """
    Get the active space selected in input (meaning SCF order),
    the number of orbitals in RAS1 and the HOMO (or SOMO) position
    :param: lines_input
    :return: active_orbitals, RAS_OCC, HOMO_orbital
    """
    with open(RAS_input, encoding="utf-8") as file:
        data = file.readlines()

    searches = ['RAS_ACT_ORB']
    elements = []

    for line in data:
        if any(i in line for i in searches):
            line = line.replace('[', '')
            line = line.replace(']', '')
            line = line.replace(',', ' ')
            line = line.replace('!user-selected RAS2 orbitals', '')

            element = line[12:]
            elements.append(element.split())
            break
    active_orbitals = np.array(elements[0], dtype=int)
    active_orbitals = np.transpose(active_orbitals)

    searches = ['RAS_OCC']
    RAS_OCC = 0
    for line in data:
        if any(i in line for i in searches):
            split_line = line.split()
            RAS_OCC = int(split_line[1])
            break

    searches = ['beta electrons']
    HOMO_orbital = 0
    for line in data:
        if any(i in line for i in searches):
            element = line[13:20]
            HOMO_orbital = int(element)
            break

    return active_orbitals, RAS_OCC, HOMO_orbital

def get_highest_amplitudes(file):

    def sort_amplitudes(amplitudes, input_orbitals):
        amplitudes = abs(amplitudes)
        change_order = np.zeros(len(amplitudes))
        change_list = []

        for i in range(0, len(amplitudes)):
            for j in range(i, len(amplitudes)):
                if (amplitudes[j] > amplitudes[i]):
                    change_order[0] = amplitudes[j]
                    amplitudes[j] = amplitudes[i]
                    amplitudes[i] = change_order[0]

                    change_list = input_orbitals[j]
                    input_orbitals[j] = input_orbitals[i]
                    input_orbitals[i] = change_list
        return amplitudes, input_orbitals

    def other_important_orbitals(amplitudes, amplitude_cutoff):
        cut_off = (1 - amplitude_cutoff) * amplitudes[0]

        indexes_list = []
        indexes_list.append('0')

        for i in range(1, len(amplitudes)):
            # Value in 0 is ignored since amplitudes are sorted
            if (amplitudes[i] > cut_off):
                indexes_list.append(i)

        main_indexes = np.array(indexes_list, dtype=int)
        return main_indexes

    next_line = next(file)
    next_line = next(file)

    elements_amplitude = []
    orbitals_each_state = []

    while next_line.startswith(" |"):  # Take all the line
        next_line = next_line.replace('|       |', '-1')  # "-1" to identify where there is no hole-part orbital
        next_line = next_line.replace('|', '')

        element = next_line.split()[-1]  # Take the amplitudes value
        elements_amplitude.append(element)

        element = next_line.split()[0:4]  # Take the HOLE, alpha, beta, PART sections
        orbitals_each_state.append(element)

        next_line = next(file)
    amplitudes = np.array(elements_amplitude, dtype=float)

    # Take the line data where amplitude is maximum and where difference
    # between maximum amplitude and amplitude is less than a cut-off
    amplitudes, orbitals_each_state = sort_amplitudes(amplitudes, orbitals_each_state)
    index_max_amplitudes = other_important_orbitals(amplitudes, amplitude_cutoff=0.3)  # Cut-off value between 1st and 2nd highest amplitudes

    return index_max_amplitudes, orbitals_each_state

def get_orbital(HOMO_orbital, configuration_data, initial_active_orbitals):
    """
    Take orbitals depending on hole, active or particle configurations
     """

    def from_RAS_to_SCF_order(HOMO_orbital, initial_SCF_space, RAS_orbital):
        """
        Change from RAS order to SCF order.
        Example: if in (1,2,3,4) orbitals the active space selected is [1,3],
        RAS_order is (1,2,3,4) while SCF_order is (2,1,3,4), since orbitals are
        order by RAS1-RAS2-RAS3.
        :param: HOMO_orbital, initial_SCF_space, final_RAS_space
        :return: final_SCF_space
        """
        initial_SCF_space = initial_SCF_space.tolist()
        SCF_orbitals_order = []

        # SCF_orbitals_order: orbitals in HF energetic order
        for orbital in range(1, HOMO_orbital + 1):  # RAS1
            if (orbital not in initial_SCF_space): SCF_orbitals_order.append(orbital)

        for orbital in initial_SCF_space:  # RAS2
            SCF_orbitals_order.append(orbital)

        for orbital in range(HOMO_orbital + 1, HOMO_orbital + 50):  # RAS3
            if (orbital not in initial_SCF_space): SCF_orbitals_order.append(orbital)

        # RAS_orbitals_order: orbitals in RAS energetic order, meaning RAS1 - RAS2 - RAS3
        RAS_orbitals_order = list(range(1, HOMO_orbital + 50))

        # print('RAS_orbitals_order', 'SCF_orbitals_order')
        # for i in range(0, len(RAS_orbitals_order)):
        #     print(RAS_orbitals_order[i], '--', SCF_orbitals_order[i])
        # exit()

        indx = RAS_orbitals_order.index(RAS_orbital)
        final_SCF_orbital = SCF_orbitals_order[indx]
        # print(RAS_orbital, final_SCF_orbital)
        return final_SCF_orbital

    # When orbital is in active space
    if (configuration_data[0] == '-1') and (configuration_data[3] ==  '-1'):
        alpha_config = str(configuration_data[1])
        beta_config = str(configuration_data[2])

        for active_orbital_index in range(0, len(alpha_config)):
            if alpha_config[active_orbital_index] != beta_config[active_orbital_index]: break
        new_orbital = initial_active_orbitals[active_orbital_index]

    # When orbital is in hole
    elif (configuration_data[0] !=  '-1') and (configuration_data[3] ==  '-1'):
        RAS_orbital = int( configuration_data[0] )
        new_orbital = from_RAS_to_SCF_order(HOMO_orbital, initial_active_orbitals, RAS_orbital)

    # When orbital is in particle
    elif (configuration_data[0] ==  '-1') and (configuration_data[3] !=  '-1'):
        RAS_orbital = int( configuration_data[3] )
        new_orbital = from_RAS_to_SCF_order(HOMO_orbital, initial_active_orbitals, RAS_orbital)

    new_orbital = int(new_orbital)

    return new_orbital

def print_excited_states(presentation_list,n_states,hole_contributions,part_contributions,SOCC_values,excitation_energies_eV,state_symmetries,new_orbital,orbital_momentum):
    """
    Prepare te presentation list with the values of each excited state
    :param: presentation_list,n_states,hole_contributions,part_contributions,SOCC_values,excitation_energies_eV,state_symmetries,new_orbital
    :return: presentation_list, SOC
     """
    state = np.round(int(n_states), 0) + 1
    symmetry = state_symmetries[n_states]

    hole = np.round(float(hole_contributions[n_states]), 2)
    part = np.round(float(part_contributions[n_states]), 2)

    SOC = np.round(float(SOCC_values[n_states]), 0)
    excit_energy = np.round(float(excitation_energies_eV[n_states]), 3)

    L_ground_state = np.round(float(orbital_momentum[n_states]), 3)

    presentation_list.append([state, symmetry, hole, part, excit_energy, new_orbital, SOC, L_ground_state])

    return presentation_list, SOC

def get_new_active_space_electrons(new_active_space, HOMO_orbital):
    """
    Electrons in the new active space considering if they are occupied
    or unoccupied observing the HOMO position
     """
    electrons = 0
    for i in range(0, len(new_active_space)):
        if new_active_space[i] <= HOMO_orbital:
            electrons += 2
        else:
            pass

    if electrons <= 0:
        electrons = 0
    else:
        electrons = electrons - 1

    return electrons

def get_excited_states_analysis(input):
    """
    Obtention of a matrix with several data for each excited state. Cut-off determines the amplitude
    difference between the 1st and 2nd (or more) important configurations in each state: if it is lower,
    the state will appear another time with the 2nd configuration.
    :param: input
    :return: excited_states_presentation_matrix
    """
    totalstates = get_number_of_states(input)

    state_symmetries, ordered_state_symmetries = get_symmetry_states(input, totalstates)

    states_ras = get_selected_states(input, totalstates, states_ras=0, selected_states = 1, symmetry_selection='None')

    eigenenergies_ras, excitation_energies_ras = get_eigenenergies(input, totalstates, states_ras)

    hole_contributions, part_contributions = get_hole_part_contributions(input, totalstates)

    SOCC_values = get_SOCC_values(input, totalstates)

    orbital_momentum = get_ground_state_orbital_momentum(input, totalstates)

    initial_active_orbitals, RAS_OCC, HOMO_orbital = get_RAS2_and_RAS_OCC_and_HOMO(input)

    excited_states_presentation_list = []
    excited_states_presentation_list.append(['State', 'Symmetry', 'Hole', 'Part', 'Excitation energy (eV)','Orbitals',
                                             'SOCC (cm-1)', 'Orbital momentum'])

    searches = ' | HOLE  | '
    n_states = 0

    with open(input, encoding="utf-8") as file:
        for line in file:
            if searches in line: # Go to configurations line

                index_max_amplitudes, state_orbitals = get_highest_amplitudes(file)

                for i in index_max_amplitudes:
                    new_orbital = get_orbital(HOMO_orbital, state_orbitals[i], initial_active_orbitals)

                    excited_states_presentation_list, SOC = print_excited_states(excited_states_presentation_list,
                                                                                 n_states,hole_contributions,
                                                                                 part_contributions,SOCC_values,
                                                                                 excitation_energies_ras * 27.211399,
                                                                                 ordered_state_symmetries,new_orbital,
                                                                                 orbital_momentum)

                n_states += 1

    excited_states_presentation_matrix = np.array(excited_states_presentation_list, dtype=object)

    print("------------------------")
    print(" EXCITED STATES ANALYSIS")
    print("------------------------")
    print('Most important settings for each state (amplitude_cutoff: 0.3) :')
    print('\n'.join( [''.join(['{:^20}'.format(item) for item in row])\
                     for row in (excited_states_presentation_matrix[:,:]) ] ) )

    print(" ")


def improved_active_space(input):
    from g_read import get_number_of_states, get_eigenenergies, get_selected_states, get_symmetry_states, \
        get_hole_part_contributions, get_SOCC_values

    from g_excited_states_analysis import get_RAS2_and_RAS_OCC_and_HOMO, get_highest_amplitudes, get_orbital, \
        get_new_active_space_electrons, print_excited_states

    totalstates = get_number_of_states(input)

    state_symmetries, ordered_state_symmetries = get_symmetry_states(input, totalstates)

    states_ras = get_selected_states(input, totalstates, states_ras=0, selected_states=1, symmetry_selection='None')

    eigenenergies_ras, excitation_energies_ras = get_eigenenergies(input, totalstates, states_ras)

    hole_contributions, part_contributions = get_hole_part_contributions(input, totalstates)

    orbital_momentum = get_ground_state_orbital_momentum(input, totalstates)

    SOCC_values = get_SOCC_values(input, totalstates)

    initial_active_orbitals, RAS_OCC, HOMO_orbital = get_RAS2_and_RAS_OCC_and_HOMO(input)

    # TAKE STATES WITH HIGHEST CONTRIBUTION
    new_space_list = []
    final_SCF_ordered_space = []
    excited_states_presentation_list = []
    excited_states_presentation_list.append(['State', 'Symmetry', 'Hole', 'Part', 'Excitation energy (eV)','Orbitals',
                                             'SOCC (cm-1)','Orbital momentum'])

    searches = ' | HOLE  | '
    n_states = 0

    with open(input, encoding="utf-8") as file:
        for line in file:
            if searches in line:  # Go to configurations line

                index_max_amplitudes, state_orbitals = get_highest_amplitudes(file)

                for i in index_max_amplitudes:
                    new_orbital = get_orbital(HOMO_orbital, state_orbitals[i], initial_active_orbitals)

                    if new_orbital not in new_space_list:
                        new_space_list.append(new_orbital)
                        excited_states_presentation_list, SOC = print_excited_states(excited_states_presentation_list,
                                                                                     n_states, hole_contributions,
                                                                                     part_contributions,SOCC_values,
                                                                                     excitation_energies_ras * 27.211399,
                                                                                     state_symmetries, new_orbital,
                                                                                     orbital_momentum)

                        if SOC != 0 or i == 0:
                            final_SCF_ordered_space.append(new_orbital)

                n_states += 1

    new_active_space = np.array(new_space_list, dtype=int)

    print("------------------------")
    print(" IMPROVED ACTIVE SPACE ")
    print("------------------------")

    print('Initial active space (HOMO =', HOMO_orbital, '):')
    initial_active_orbitals = initial_active_orbitals.tolist()
    electrons = get_new_active_space_electrons(initial_active_orbitals, HOMO_orbital)
    print('[',electrons, ',', len(initial_active_orbitals), '] ;', initial_active_orbitals)

    # print('')
    # print('New active space (SOCC not zero, HOMO singly occupied):')
    # electrons = get_new_active_space_electrons(final_SCF_ordered_space, HOMO_orbital)
    # final_SCF_ordered_space.sort()
    # print('[',electrons, ',', len(final_SCF_ordered_space), '] ;', final_SCF_ordered_space)

    print('')
    print('Initial active space + New active space:')

    initial = set(initial_active_orbitals)
    final = set(final_SCF_ordered_space)
    in_final_but_not_in_initial = final - initial
    initial_plus_final_space = initial_active_orbitals + list(in_final_but_not_in_initial)

    electrons = get_new_active_space_electrons(initial_plus_final_space, HOMO_orbital)
    initial_plus_final_space.sort()
    print('[',electrons, ',', len(initial_plus_final_space), '] ;', initial_plus_final_space)

import numpy as np


def change_RAS_with_eom_energies(RAS_states_to_change, eom_states_to_change, excitation_energies_ras,
                                 eom_excitation_energies):
    """
    Change the RAS energies with the eom energies of the states selected.
    :param RAS_states_to_change,eom_states_to_change,excitation_energies_ras,eom_excitation_energies
    :return SOC_matrix: SOC matrix
    """
    changed_excitation_energies = np.array(excitation_energies_ras, dtype=float)

    for i in range(0, len(RAS_states_to_change)):
        ras_index = RAS_states_to_change[i] - 1
        eom_index = eom_states_to_change[i] - 1
        changed_excitation_energies[ras_index] = eom_excitation_energies[eom_index]
    return changed_excitation_energies


def get_comparison_results(RAS_states_to_change, eom_states_to_change, ras_energies, eom_energies, soc_constant_ras,
                           selected_eom_soc_constants):
    comparison_presentation_list = []

    comparison_presentation_list.append(
        ['RAS State', 'eom transition', 'RAS excitation energy (eV)', 'eom-CC excitation energy (eV)',
         'RAS SOCC (cm-1)', 'eom-CC SOCC (cm-1)'])

    # comparison_presentation_list.append(['State', 'RAS excitation energy (eV)', 'eom-CC excitation energy (eV)'])

    for i in range(0, len(RAS_states_to_change)):
        ras_index = RAS_states_to_change[i]
        eom_index = eom_states_to_change[i]

        ras_excit_ener = np.round(float(ras_energies[ras_index - 1]) * 27.211399, 4)
        eom_excit_ener = np.round(float(eom_energies[eom_index - 1]) * 27.211399, 4)

        ras_soc = np.round(float(soc_constant_ras[ras_index - 1]), 0)
        eom_soc = np.round(float(selected_eom_soc_constants[eom_index - 1]), 0)

        comparison_presentation_list.append([ras_index, eom_index, ras_excit_ener, eom_excit_ener, ras_soc, eom_soc])
    return comparison_presentation_list


def ras_and_eom_energy_exchange(eom_input, input, states_ras, RAS_states_to_change, eom_states_to_change):
    from g_read import get_number_of_states, get_eigenenergies, get_selected_states, get_spin_orbit_couplings, \
        get_SOCC_values

    from g_operations import from_energies_SOC_to_g_values

    from g_take_eom_states import get_eom_type, get_SCF_energy, get_irreps_energies, prepare_presentation_list, \
        get_maximum_amplitude_orbitals, get_eom_SOCC_values, eom_results

    from g_ras_eom_change import change_RAS_with_eom_energies, get_comparison_results

    eom_version = get_eom_type(eom_input)

    RAS_SCF_reference_energy = get_SCF_energy(input)

    eom_SCF_reference_energy = get_SCF_energy(eom_input)

    irreps, optimized_state_index, total_energies_eom, all_excitation_energies_eom = get_irreps_energies(eom_input)

    eom_SOCC = get_eom_SOCC_values(eom_input, optimized_state_index)

    orbitals = get_maximum_amplitude_orbitals(eom_input, eom_version)

    eom_presentation_list = prepare_presentation_list(eom_version)

    eom_presentation_list, selected_excitation_energies_eom, selected_eom_soc_constants = eom_results(
        eom_presentation_list, irreps, total_energies_eom, all_excitation_energies_eom, eom_SOCC, orbitals)

    # G-tensor calculation section

    totalstates = get_number_of_states(input)

    states_ras = get_selected_states(input, totalstates, states_ras, selected_states=1, symmetry_selection='None')

    eigenenergies_ras, excitation_energies_ras = get_eigenenergies(input, totalstates, states_ras)

    changed_excitation_energies = change_RAS_with_eom_energies(RAS_states_to_change, eom_states_to_change,
                                                               excitation_energies_ras,
                                                               selected_excitation_energies_eom)

    SOC_ras = get_spin_orbit_couplings(input, totalstates, states_ras, selected_SOC=0)

    soc_constant_ras = get_SOCC_values(input, totalstates)

    G_tensor, G_tensor_results = from_energies_SOC_to_g_values(input, states_ras, totalstates,
                                                               changed_excitation_energies, SOC_ras)

    comparison_presentation_list = get_comparison_results(RAS_states_to_change, eom_states_to_change,
                                                          excitation_energies_ras, selected_excitation_energies_eom,
                                                          soc_constant_ras, selected_eom_soc_constants)

    print("")
    print("--------------------------------------")
    print("    eom-CC + RAS change in energy ")
    print("--------------------------------------")

    print('RAS SCF reference energy:', RAS_SCF_reference_energy)
    print('eom SCF reference energy:', eom_SCF_reference_energy)
    print('')

    print('RAS states to be changed:', RAS_states_to_change)
    print('eom states to be changed:', eom_states_to_change)
    print('')

    return comparison_presentation_list, G_tensor, G_tensor_results

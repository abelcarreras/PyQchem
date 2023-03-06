#####################################
#          MODULES SELECTION
#####################################
# module load Python/3.7.4-Anaconda3-2019.10
import numpy as np
import sys

from g_read import get_number_of_states, get_eigenenergies, get_selected_states, get_spin_orbit_couplings, \
    get_SOCC_values, get_ground_state_orbital_momentum, get_symmetry_states, get_spin_orbit_couplings_pyqchem

from g_operations import from_energies_SOC_to_g_values, print_g_calculation

from g_excited_states_analysis import get_excited_states_analysis, improved_active_space

# from g_take_eom_states import get_eom_transitions_analysis
#
# from g_ras_eom_change import ras_and_eom_energy_exchange

from g_plots import get_bar_chart, sos_analysis_and_plot

#####################################
#            INPUT VALUES
#####################################

# G-TENSOR CALCULATION
g_calculation = 1
ras_input = '../../../../PycharmProjects/my_g-tensor/' \
'RASCI_results/cucl4_2-/cucl4_2-_def2tzvp_11_6.out' # str(sys.argv[1])

selected_states = 1 # 0: use "state_ras" ; 1: use all states ; 2: use states by selected symmetry
states_ras = [1,2,3,4,5,6,7,8,9,10,11] # States to be included when "selected_states = 0"
symmetry_selection = 'B2u' # Symmetry selected states
selected_SOC = 0 # 0: Total mean-field SOC matrix; 1: 1-elec SOC matrix; 2: 2-elec mean-field SOC matrix

# EXCITED STATES ANALYSIS IN ras
excited_states_analysis = 0
new_active_space = 0
sos_analysis = 0
bar_plots = 0

# eom ANALYSIS AND ras-eom ENERGIES EXCHANGE
eom_information = 0
eom_input = 'Sven_eom-MP2/ag_cucl4_2-_D2h_eomip_soc_def2-TZVP.in.out'

eom_change_energies = 0
ras_states_to_change = [2,3,4]
eom_states_to_change = [4,7,9]

# OUTPUT
write_file = 0 # 0: write results directly; 1: write in output file
output_file = ras_input + '-gvalues.txt'
if (write_file == 1):
    sys.stdout = open(output_file, "w")

#####################################
#      G-VALUE CALCULATION
#####################################
if (g_calculation == 1):

    totalstates = get_number_of_states(ras_input)

    states_ras = get_selected_states(ras_input, totalstates, states_ras, selected_states, symmetry_selection)

    eigenenergies_ras, excitation_energies_ras = get_eigenenergies(ras_input, totalstates, states_ras)

    soc_ras = get_spin_orbit_couplings_pyqchem(ras_input, totalstates, states_ras, selected_SOC)

    ras_G_matrix, ras_g_values, eigenvalues, eigenvector = from_energies_SOC_to_g_values(ras_input, states_ras,
                                                                                         totalstates,
                                                                                         excitation_energies_ras,
                                                                                         soc_ras)

    print_g_calculation(ras_input, totalstates, selected_states, symmetry_selection, states_ras, ras_g_values)

    # for i in range(0, len(eigenvalues),2):
    #     energy_diff = (eigenvalues[i] - excitation_energies_ras[i//2]) * 27.211399
    #     print('State:', i//2, '; Relativ. - non-relativ. energy diff (eV): ', np.round(energy_diff,4), '; c**2: ',
    #           np.round(np.absolute(eigenvector[i,i]),3) )


#####################################
#        EXCITED STATE ANALYSIS
#####################################
if (excited_states_analysis == 1):
    get_excited_states_analysis(ras_input)

if (new_active_space == 1):
    improved_active_space(ras_input)

#####################################
#        PLOT ANALYSIS
#####################################

if (sos_analysis == 1):
    sos_analysis_and_plot(ras_input, save_picture=0)

if (bar_plots == 1):

    totalstates = get_number_of_states(ras_input)

    states_ras = get_selected_states(ras_input, totalstates, states_ras, selected_states=1, symmetry_selection='None')

    eigenenergies_ras, excitation_energies_ras = get_eigenenergies(ras_input, totalstates, states_ras)

    soc_ras = get_spin_orbit_couplings(ras_input, totalstates, states_ras, selected_SOC)

    ras_G_matrix, ras_g_values, eigenvalues, eigenvector = from_energies_SOC_to_g_values(ras_input, states_ras, totalstates,
                                                               excitation_energies_ras, soc_ras)

    # Printing excitation energies versus orbital symmetries:
    state_symmetries, ordered_state_symmetries = get_symmetry_states(ras_input, totalstates)

    get_bar_chart(ras_input, ordered_state_symmetries, excitation_energies_ras * 27.211399, 'State', 'Energy (eV)',
                  'Excitation energies')

    # Printing SOCC versus orbital symmetries:
    SOCC_values = get_SOCC_values(ras_input, totalstates)
    get_bar_chart(ras_input, ordered_state_symmetries, SOCC_values / 8065.540107, 'State', 'Energy (eV)',
                  'Mean-field spin-orbit coupling constants')

    # Printing orbital angular momentum versus orbital symmetries:
    orbital_momentum = get_ground_state_orbital_momentum(ras_input, totalstates)
    get_bar_chart(ras_input, ordered_state_symmetries, orbital_momentum, 'State', 'Orbital angular momentum','')

#####################################
#      eom COMPARISON
#####################################
# if eom_information == 1:
#     get_eom_transitions_analysis(eom_input)
#
# if (eom_change_energies == 1):
#     comparison_presentation_list, G_tensor, G_tensor_results = ras_and_eom_energy_exchange(eom_input, ras_input, states_ras, ras_states_to_change, eom_states_to_change)
#
#     print('g-factor (x y z dimensions) with ras energies:')
#     print(np.round(ras_g_values.real[0], 3), np.round(ras_g_values.real[1], 3),np.round(ras_g_values.real[2], 3))
#     print('')
#
#     print('g-factor (x y z dimensions) with eom eigenenergies:')
#     print(np.round(G_tensor_results.real[0], 3), np.round(G_tensor_results.real[1], 3),np.round(G_tensor_results.real[2], 3))
#     print('')
#
#     print("ras file selected: ", ras_input)
#     print("eom-CC file selected: ", eom_input)
#     print('\n'.join([''.join(['{:^30}'.format(item) for item in row]) \
#                      for row in (comparison_presentation_list)]))

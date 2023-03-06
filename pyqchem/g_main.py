"""
g-tensor calculation by USING ONLY GROUND STATE DOUBLETS
"""

from g_read import get_number_of_states, get_eigenenergies, get_selected_states, get_spin_orbit_couplings_pyqchem

from g_operations import from_energies_SOC_to_g_values, print_g_calculation

#####################################
#            INPUT VALUES
#####################################

ras_input = '../scripts/h2o.out'

selected_states = 1 # 0: use "state_ras" ; 1: use all states ; 2: use states by selected symmetry
states_ras = [1,2,3,4,5] # States to be included when "selected_states = 0"
symmetry_selection = 'B2u' # Symmetry selected states

#####################################
#      G-VALUE CALCULATION
#####################################

totalstates = get_number_of_states(ras_input)

states_ras = get_selected_states(ras_input, totalstates, states_ras, selected_states, symmetry_selection)

eigenenergies_ras, excitation_energies_ras = get_eigenenergies(ras_input, totalstates, states_ras)

soc_ras = get_spin_orbit_couplings_pyqchem(ras_input, totalstates, states_ras)

ras_G_matrix, ras_g_values, eigenvalues, eigenvector = from_energies_SOC_to_g_values(ras_input, states_ras,
                                                                                     totalstates,
                                                                                     excitation_energies_ras,
                                                                                     soc_ras)

print_g_calculation(ras_input, totalstates, selected_states, symmetry_selection, states_ras, ras_g_values)

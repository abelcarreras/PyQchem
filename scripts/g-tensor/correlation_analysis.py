"""
PROGRAM TO CALCULATE THE CORRELATION BETWEEN ENERGIES-SOC
AND THE G-TENSOR VALUES
"""
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import stats

from g_read import get_number_of_states, get_eigenenergies, get_selected_states, get_spin_orbit_couplings, \
    get_SOCC_values, get_spin_matrices, get_orbital_matrices

from g_operations import get_Hamiltonian_construction, Hamiltonian_diagonalization, angular_matrixes_obtention, \
    angular_matrixes_obtention, g_factor_calculation, from_energies_SOC_to_g_values, print_g_calculation

from g_excited_states_analysis import get_excited_states_analysis, improved_active_space

from g_take_eom_states import get_eom_transitions_analysis

from g_ras_eom_change import ras_and_eom_energy_exchange

from g_plots import get_bar_chart, sos_analysis_and_plot, plot_g_tensor_vs_states

#####

ras_input = 'g-tensor_final_results/no2/no2_def2-TZVP_11_9.out'  # str(sys.argv[1])
selected_states = 0  # 0: use "state_ras" ; 1: use all states ; 2: use states by selected symmetry
states_ras = [1, 2]  # States to be included when "selected_states = 0"
symmetry_selection = 'B1'  # Symmetry selected states
selected_SOC = 0  # 0: Total mean-field SOC matrix; 1: 1-elec SOC matrix; 2: 2-elec mean-field SOC matrix

totalstates = get_number_of_states(ras_input)

states_ras = get_selected_states(ras_input, totalstates, states_ras, selected_states, symmetry_selection)

eigenenergies_ras, excitation_energies_ras = get_eigenenergies(ras_input, totalstates, states_ras)

soc_ras = get_spin_orbit_couplings(ras_input, totalstates, states_ras, selected_SOC)

SOCC_values = get_SOCC_values(ras_input, totalstates)

ras_G_matrix, ras_g_values = from_energies_SOC_to_g_values(ras_input, states_ras, totalstates, excitation_energies_ras,
                                                           soc_ras)

print_g_calculation(ras_input, totalstates, selected_states, symmetry_selection, states_ras, ras_g_values)


def soc_gvalues_correlation(ras_input, states_ras, totalstates, excitation_energies_ras, soc_ras,SOCC_values):
    """"
    Correlation analysis between SOC and gvalues. Change is done only for imaginary part of SOC
    between ground and first excited state.
    :param: ras_input, states_ras, totalstates, excitation_energies_ras, soc_ras,SOCC_values
    :return: soc_gvalues_matrix, r_square of fit
    """
    presentation_tuple = []

    max_socc = (max(SOCC_values)) * 3
    min_socc = (min(SOCC_values)) / 20

    for soc_value in np.linspace(min_socc, max_socc , 50):
        soc_value = 0 + soc_value * 1j / 219474.63068 # from cm-1 to au

        soc_ras[0, 3] = soc_value
        soc_ras[1,2] = soc_value
        soc_ras[3,0] = np.conj(soc_ras[0, 3])
        soc_ras[2,1] = np.conj(soc_ras[1,2])

        ras_G_matrix, ras_g_values = from_energies_SOC_to_g_values(ras_input, states_ras, totalstates, excitation_energies_ras, soc_ras)

        presentation_tuple.append([abs(soc_value), np.round(ras_g_values.real[0], 3), np.round(ras_g_values.real[1], 3), np.round(ras_g_values.real[2], 3)])

    presentation_matrix = np.array(presentation_tuple, dtype=object)
    presentation_matrix[:, 0] = presentation_matrix[:,0] * 27.211399 # from au to eV

    plot_g_tensor_vs_states(presentation_matrix, x_title='SOC (eV)', y_title='g-values (ppt)',
                            main_title='g-factor evolution with SOC', save_picture=0)
    plt.show()
    plt.close()

    # R-squares obtention
    r_square = []

    for i in range(1, len(presentation_matrix[0]) ):
        res = stats.linregress(presentation_matrix[:,0].astype(float), presentation_matrix[:,i].astype(float))
        r_square.append(res.rvalue ** 2)

    return presentation_matrix, r_square


def energy_gvalues_correlation(ras_input, states_ras, totalstates, excitation_energies_ras, soc_ras):
    """"
    Correlation analysis between energies and gvalues. Change is done only for energy
    between ground and first excited state.
    :param: ras_input, states_ras, totalstates, excitation_energies_ras, soc_ras,SOCC_values
    :return: soc_gvalues_matrix, r_square of fit
    """
    presentation_tuple = []

    max_ener = (max(excitation_energies_ras)) * 5
    min_ener = (max(excitation_energies_ras)) / 20

    for ener_value in np.linspace(min_ener, max_ener, 50):
        excitation_energies_ras[1] = ener_value

        ras_G_matrix, ras_g_values = from_energies_SOC_to_g_values(ras_input, states_ras, totalstates,
                                                                   excitation_energies_ras, soc_ras)

        presentation_tuple.append([ener_value * 27.211399, np.round(ras_g_values.real[0], 3), np.round(ras_g_values.real[1], 3),
                                   np.round(ras_g_values.real[2], 3)])
        print(ener_value * 27.211399, np.round(ras_g_values.real[0], 3), np.round(ras_g_values.real[1], 3),
                                   np.round(ras_g_values.real[2], 3))

    presentation_matrix = np.array(presentation_tuple, dtype=object)

    plot_g_tensor_vs_states(presentation_matrix, x_title='Energy (eV)', y_title='g-values (ppt)',
                            main_title='g-factor evolution with energy', save_picture=0)

    plt.show()
    plt.close()

    return presentation_matrix


def orbital_momentum_gvalues_correlation(input, states_ras, totalstates, excitation_energies_ras, soc_ras):
    """"
    Correlation analysis between energies and gvalues. Change is done only for energy
    between ground and first excited state.
    :param: ras_input, states_ras, totalstates, excitation_energies_ras, soc_ras,SOCC_values
    :return: soc_gvalues_matrix, r_square of fit
    """
    Hamiltonian_ras = get_Hamiltonian_construction(states_ras, excitation_energies_ras, soc_ras)

    eigenvalues, eigenvector, kramers_states = Hamiltonian_diagonalization(Hamiltonian_ras)

    spin_matrix = get_spin_matrices(input, states_ras)

    sigma_matrix = angular_matrixes_obtention(eigenvalues, eigenvector, kramers_states, spin_matrix)

    l_matrix = get_orbital_matrices(input, totalstates, states_ras)

    # print('Angular momentums (x,y,z):')
    # for k in range(0,3):
    #    print('Dimension: ', k)
    #    print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
    #                     for row in np.round((l_matrix[:,:,k]),5)]))
    #    print(" ")
    # exit()

    # Change L values

    presentation_tuple = []

    max_momentum = 1.5
    min_momentum = 0.01

    for l_value in np.linspace(min_momentum, max_momentum, 100):
        l_value = 0 + l_value * 1j

        l_matrix[2, 0, 0] = l_value
        l_matrix[3, 1, 0] = l_value
        l_matrix[0, 2, 0] = np.conj(l_matrix[2, 0, 0])
        l_matrix[1, 3, 0] = np.conj(l_matrix[3, 1, 0])

        lambda_matrix = angular_matrixes_obtention(eigenvalues, eigenvector, kramers_states, l_matrix)

        ras_G_matrix, ras_g_values = g_factor_calculation(lambda_matrix, sigma_matrix)

        presentation_tuple.append([abs(l_value), np.round(ras_g_values.real[0], 3), np.round(ras_g_values.real[1], 3),
                                   np.round(ras_g_values.real[2], 3)])

    presentation_matrix = np.array(presentation_tuple, dtype=object)

    plot_g_tensor_vs_states(presentation_matrix, x_title='Orbital Angular Momentum', y_title='g-values (ppt)',
                            main_title='g-factor evolution with orbital angular momentum', save_picture=0)

    plt.show()
    plt.close()

    # R-squares obtention
    r_square = []

    for i in range(1, len(presentation_matrix[0]) ):
        res = stats.linregress(presentation_matrix[:,0].astype(float), presentation_matrix[:,i].astype(float))
        r_square.append(res.rvalue ** 2)

    return presentation_matrix, r_square

# soc_gvalues_matrix, r2 = soc_gvalues_correlation(ras_input,states_ras,totalstates,excitation_energies_ras,soc_ras,SOCC_values)
# print('R-square in three dimmensions: ', np.round(r2[0], 2), np.round(r2[1], 3), np.round(r2[2], 3))
#
# ener_gvalues_matrix = energy_gvalues_correlation(ras_input,states_ras,totalstates,excitation_energies_ras,soc_ras)

momentum_gvalues_matrix, r2_momentum = orbital_momentum_gvalues_correlation(ras_input, states_ras, totalstates, excitation_energies_ras, soc_ras)
print('R-square in three dimmensions: ', np.round(r2_momentum[0], 2), np.round(r2_momentum[1], 3), np.round(r2_momentum[2], 3))

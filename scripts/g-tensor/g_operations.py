"""
    G-TENSOR CALCULATION WITH ORBITAL AND SPIN
    ANGULAR MOMENTUM AND SOC BETWEEN STATES
 Analysis of the excited states and print of those
 orbitals involved in the configurations with the
 highest amplitudes in the excited states
"""
import numpy as np
from numpy import linalg, sqrt


def get_Hamiltonian_construction(selected_states, eigenenergies, spin_orbit_coupling):
    """
    Hamiltonian is written 'bra' in rows and 'ket' in columns,
    with spin order -1/2 , +1/2.
    :param: selected_states, eigenenergies, spin_orbit_coupling
    :return: Hamiltonian
    """
    Hamiltonian = np.zeros((len(selected_states) * 2, len(selected_states) * 2), dtype=complex)
    # spin_orbit_coupling[:,:] = 0

    for i in range(0, len(selected_states) * 2):
        for j in range(0, len(selected_states) * 2):
            if (i == j):
                Hamiltonian[i, i] = eigenenergies[i // 2]
            else:
                Hamiltonian[i, j] = spin_orbit_coupling[i, j]

    # print('Hamiltonian (SOC in cm-1, energies in a.u):')
    # print('\n'.join([''.join(['{:^20}'.format(item) for item in row])\
    #                  for row in np.round((Hamiltonian[:,:]),5)]))
    # print(" ")
    # exit()

    return Hamiltonian


def Hamiltonian_diagonalization(Hamiltonian):
    """
    1) Hamiltonian is diagonalized
    2) eigenvectors-eigenvalues are ordered by weight coefficients
    3) Doublet with the lowest energy is set as the new basis (Kramer doublet)
    :param: Hamiltonian
    :return eigenvalues, eigenvectors, kramer_st: kramer_st is the index
    of the state set as Kramer doublet * 2
    """
    # def normalize_vector(vector):
    #     """
    #     Normalize a complex vector
    #     :param vector:
    #     :return:
    #     """
    #     sum_of_coefficients = 0
    #     for i in range(0, len(vector)):
    #         coefficient = vector.real[i]**2 + vector.imag[i]**2
    #         sum_of_coefficients += coefficient
    #
    #     module = sqrt( sum_of_coefficients )
    #     vector[:] = vector[:] / module
    #     return vector, module

    eigenvalues, eigenvectors = linalg.eigh(Hamiltonian)

    # Reorder eigenvectors (and eigenenergies) by weight coefficients
    change_order = np.zeros(len(eigenvectors), dtype=complex)
    for v_1 in range(0, len(eigenvectors)):
        for v_2 in range(v_1, len(eigenvectors)):

            if (abs(eigenvectors[v_1, v_2]) > abs(eigenvectors[v_1, v_1])):
                change_order[:] = eigenvectors[:, v_1]
                eigenvectors[:, v_1] = eigenvectors[:, v_2]
                eigenvectors[:, v_2] = change_order[:]

                change_order.real[0] = eigenvalues[v_1]
                eigenvalues[v_1] = eigenvalues[v_2]
                eigenvalues[v_2] = change_order.real[0]

    # Normalize the eigenvectors
    # for i in range(0, len(eigenvalues)):
    #     print('Before eigenvector:', eigenvectors[:, i])
    #     module = sqrt(sum( ( abs(eigenvectors[:, i]) )** 2) )
    #     eigenvectors[:, i] = eigenvectors[:, i] / module
    #     print('After eigenvector:', eigenvectors[:, i])
    #     print('module:', module)
    #     print()
    # print('-------------------------------------')
    # for i in range(0, len(eigenvalues)):
    #     eigenvectors[:, i],module = normalize_vector(eigenvectors[:,i])
    #     print('eigenvalue:', eigenvalues[i])
    #     print('eigenvector:', eigenvectors[:, i])
    #     print('module:', module)
    #     print()

    # Kramer doublets selection:
    minimum_energy = min(eigenvalues)
    eigenvalues_list = list(eigenvalues)
    kramer_st = eigenvalues_list.index(minimum_energy)

    # The index of the selected state must be even since
    # Kramer doublets are [kramer_st, kramer_st+1]
    if (kramer_st % 2) != 0:
        kramer_st = kramer_st - 1

    # for i in range(0, len(eigenvalues)):
    #     print('eigenvalue:', eigenvalues[i])
    #     print('eigenvector:', eigenvectors[:, i])
    #     print()

    return eigenvalues, eigenvectors, kramer_st


def angular_matrixes_obtention(eigenvalues, eigenvectors, kramer_st, input_angular_matrix):
    """
    Spin or orbital angular matrix calculation using:
    1) coeff_bra, coeff_ket: coefficients of the lineal combination of non-relativistic states,
    that come from Kramer doublet states eigenvectors
    2) angular_value: angular momentum between states. Depending on the
    column of the final matrix, it takes real (col 0), imaginary (col 1) or
    both parts (col 2).

    :param eigenvalues, eigenvectors, kramer_st, input_angular_matrix
    :return angular_matrix: contains the spin value < B(S,Sz) | Sx | A(S',Sz') >,
    meaning < B(S,Sz)| corresponds to state "i" (in rows) and | A(S',Sz') > to
    state "j" (in columns)
    """
    # Matrices calculation:
    angular_matrix = np.zeros((3, 3), dtype=complex)

    for row in range(0, 3):  # dimension x,y,z
        for column in range(0, 3):  # dimension x,y,z

            for bra in range(0, len(eigenvalues)):  # state <B|
                for ket in range(0, len(eigenvalues)):  # state |A>

                    # coeff_ket for 1st and 2nd column ("x", "y" dir)
                    # coeff_ket_2 for 3rd column ("z" direction)
                    coeff_bra = np.conj(eigenvectors[bra, kramer_st + 1])
                    coeff_ket = (eigenvectors[ket, kramer_st])
                    coeff_ket_2 = (eigenvectors[ket, kramer_st + 1])

                    angular_value = (input_angular_matrix[bra, ket, row])

                    if column == 0:
                        element = coeff_bra * coeff_ket * angular_value
                        angular_matrix[row, column] += 2 * element.real

                    elif column == 1:
                        element = coeff_bra * coeff_ket * angular_value
                        angular_matrix[row, column] += 2 * element.imag

                    elif column == 2:
                        element = coeff_bra * coeff_ket_2 * angular_value
                        angular_matrix[row, column] += 2 * element

    # print('SIGMA matrix with all spin angular momentums:')
    # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
    #                  for row in np.round((angular_matrix[:,:]),8)]))
    # print(" ")

    return angular_matrix


def g_factor_calculation(lambda_matrix, sigma_matrix):
    """
    Calculation of the G-tensor with lambda and sigma matrices. Then, g-factors
    are calculated as square roots of the eigenvalues of the G-tensor.
    """
    # G-tensor matrix obtention:
    lande_factor = 2.002319304363

    sigma_plus_lambda = lande_factor * sigma_matrix + lambda_matrix

    G_matrix = np.matmul((sigma_plus_lambda), np.transpose(sigma_plus_lambda))

    # Diagonalize and reorder by weight coefficients:
    G_matrix_diagonal, rotation_matrix = linalg.eigh(G_matrix)

    change_order = np.zeros(len(G_matrix_diagonal), dtype=complex)
    for i in range(0, 3):
        for j in range(i, 3):
            if (abs(rotation_matrix[i, j]) > abs(rotation_matrix[i, i])):
                change_order[:] = rotation_matrix[:, j]
                rotation_matrix[:, j] = rotation_matrix[:, i]
                rotation_matrix[:, i] = change_order[:]

                change_order.real[0] = G_matrix_diagonal[j]
                G_matrix_diagonal[j] = G_matrix_diagonal[i]
                G_matrix_diagonal[i] = change_order.real[0]

    # Obtain the g-factor: g = O^(1r) * X * sqrt(Gdiag) * O
    # X_mat = np.identity(3)
    # g_value = (np.transpose(rotation_matrix)).dot(X_mat).dot(sqrt(G_matrix_diagonal)).dot(rotation_matrix)
    # g_value = ( (g_value) - lande_factor ) * 1000

    g_tensor_values = np.zeros(3, dtype=complex)
    for i in range(0, 3):
        g_tensor_values[i] = (sqrt(G_matrix_diagonal[i]) - lande_factor) * 1000

    # print('G_tensor matrix:')
    # print(''.join(['{:^15}'.format(item) for item in ['x','y','z']]) )
    # print('\n'.join([''.join(['{:^15}'.format(item) for item in row])\
    #                  for row in np.round((G_matrix[:,:]),5)]))
    # print('')
    # print(np.round(g_tensor_values.real[0], 3), np.round(g_tensor_values.real[1], 3), \
    #       np.round(g_tensor_values.real[2], 3))
    # print('')

    return G_matrix, g_tensor_values


def from_energies_SOC_to_g_values(input, states_ras, totalstates, excitation_energies_ras, SOC_ras):
    """"
    Obtention of the g-values from the eigenenergies and the SOCs.
    :param: input, states, totalstates, excitation_energies, SOC
    :return: G_matrix, G_tensor_results
    """
    from g_read import get_spin_matrices, get_orbital_matrices

    Hamiltonian_ras = get_Hamiltonian_construction(states_ras, excitation_energies_ras, SOC_ras)

    eigenvalues, eigenvector, kramers_states = Hamiltonian_diagonalization(Hamiltonian_ras)

    spin_matrix = get_spin_matrices(input, states_ras)

    l_matrix = get_orbital_matrices(input, totalstates, states_ras)

    sigma_matrix = angular_matrixes_obtention(eigenvalues, eigenvector, kramers_states, spin_matrix)

    lambda_matrix = angular_matrixes_obtention(eigenvalues, eigenvector, kramers_states, l_matrix)

    G_matrix, g_values = g_factor_calculation(lambda_matrix, sigma_matrix)

    return G_matrix, g_values, eigenvalues, eigenvector


def print_g_calculation(input, totalstates, selected_states, symmetry_selection, states_ras, G_tensor_results_ras):
    print("--------------------------------------")
    print("     INPUT SECTION")
    print("--------------------------------------")
    print("File selected: ", input)
    print("Number of states: ", totalstates)
    if (selected_states == 2):
        print("Symmetry: ", symmetry_selection)
        print("Selected states: ", states_ras)
    else:
        print("Selected states: ", states_ras)

    print(" ")
    print("------------------------")
    print(" ras-CI OUTPUT SECTION")
    print("------------------------")
    print('g-factor (x y z dimensions):')
    print(np.round(G_tensor_results_ras.real[0], 3), np.round(G_tensor_results_ras.real[1], 3), \
          np.round(G_tensor_results_ras.real[2], 3))
    print('')


def perturbative_method(totalstates, eigenenergies, spin_orbit_coupling):
    """
    The coefficients that previously were on the eigenvectors
    are now calculated using equation 27: they are called "eta". Similar
    process to the one used in Hamiltonian construction, BUT considering
    also the eigenvalues substraction in the denominator: en[a] - en[b],
    where the "a" corresponds to the row number and the "b" to the column number
    """
    coef = np.zeros((len(states_ras) * 2, len(states_ras) * 2), dtype=complex)
    k1 = 0
    k2 = 0
    for y in range(0, len(totalstates) * 2):
        k1 = 0
        for x in range(0, len(totalstates) * 2):
            if (y == x) or (x % 2 != 0 and x - 1 == y):
                coef[x][y] = 0
                k1 = k1 + 1
            elif (y % 2 == 0):
                coef[x][y] = spin_orbit_coupling[x - k1 + k2][0] / (eigenenergies[x] - eigenenergies[y])
                coef[x][y + 1] = spin_orbit_coupling[x - k1 + k2][1] / (eigenenergies[x] - eigenenergies[y + 1])
        k2 = k2 + len(states_ras) - 1

    # The procedure is the same than for the Hamiltonian case. Coefficients matrix
    # obtained by perturbation expressions is symmetric, meaning multiplication does not
    # change if it is done by columns or by rows
    # gpt_mat_ras = np.zeros((len(states_ras) * 2, len(states_ras) * 2, 3), dtype=complex)
    # for z in range(0, 3):
    #     for y in range(0, len(states_ras) * 2, 2):
    #         for x in range(0, len(states_ras) * 2):
    #             if (y == x) or (x % 2 != 0 and x - 1 == y):
    #                 gpt_mat_ras[y, x, z] = 0
    #                 gpt_mat_ras[x, y, z] = 0
    #             elif (x % 2 == 0):
    #                 gpt_mat_ras[x, y, z] = -4000 * etacoeff_ras[y, x] * l_matrix[x // 2, y // 2, z]
    #                 gpt_mat_ras[x + 1, y, z] = -4000 * etacoeff_ras[y, x + 1] * l_matrix[x // 2, y // 2, z]
    #                 gpt_mat_ras[x, y + 1, z] = -4000 * etacoeff_ras[y + 1, x] * l_matrix[x // 2, y // 2, z]
    #                 gpt_mat_ras[x + 1, y + 1, z] = -4000 * etacoeff_ras[y + 1, x + 1] * l_matrix[x // 2, y // 2, z]

    return

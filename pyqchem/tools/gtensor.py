__author__ = 'Antonio Cebreiro-Gallardo'

from pyqchem.tools.spin import s2_to_s, spin_matrices
from numpy import linalg, sqrt
from scipy import constants
import numpy as np

lande_factor = 2.002319304363


def take_selected_states_values(input_list, selected_states):
    """
    From "input_list", take only the elements in "selected_states" positions and put them in "output_list"
    :param: input_list, selected_states
    :return: output_list
    """
    output_list = []
    for i in selected_states:
        output_list.append(input_list[i - 1])
    return output_list


def get_eigenenergies(selected_states, data):
    """
    Get energies of RAS-CI states selected in Q-Chem output. Energies in a.u., excitation energies in eV.
    :param: file, totalstates, selected_states
    :return: eigenenergies, excitation_energies
    """

    element_energy = [state['total_energy'] for state  in data['excited_states']]
    elements_excitenergy = [state['excitation_energy'] for state  in data['excited_states']]

    energies_selected = take_selected_states_values(element_energy, selected_states)
    excited_energies_selected = take_selected_states_values(elements_excitenergy, selected_states)

    eigenenergies = np.array(energies_selected, dtype=float)
    excitation_energies_ev = np.array(excited_energies_selected, dtype=float)
    excitation_energies = excitation_energies_ev / constants.physical_constants['Hartree energy in eV'][0]
    return eigenenergies, excitation_energies


def get_states_sz(qchem_file, states_selected):
    """
    Get a list of: i) s2 of all the states ii) [-sz,..,+sz] of the highest s2 state iii) [-sz,..,+sz] of ground
    state.
    :param: qchem_file, states_selected
    :return: all_multip, all_sz, ground_sz
    """
    def from_s2_to_sz_list(s2):
        # Obtain s values
        s = 0.5 * (-1 + np.sqrt(1 + 4 * s2))
        # Making multiplicity list from -n to +n in 1/2 intervals
        s = 0.5 * np.round(s / 0.5)
        sz = list(np.arange(-s, s + 1, 1))
        return sz

    read_multip = []
    for i, state in enumerate(qchem_file['excited_states']):
        read_multip.append(np.round(state['multiplicity'], 2))

    all_multip = []
    for i in range(0, len(read_multip)):
        # Making multiplicity numbers multiples of doublets (s2=0.75).
        element = read_multip[i]
        n_times = np.round(element / 0.75)
        new_multip = 0.75 * n_times
        all_multip.append(new_multip)

    ordered_multip = []
    for i in range(0, len(states_selected)):
        selected_mult = states_selected[i] - 1
        ordered_multip.append(all_multip[selected_mult])

    all_sz = from_s2_to_sz_list(max(all_multip))
    ground_sz = from_s2_to_sz_list(ordered_multip[0])

    if len(ground_sz) == 1:
        raise ValueError("Warning! It is not allowed the calculation of the g-tensor in a singlet ground state. "
                         "Ground state corresponds to the first of the included states selected.")
    return all_multip, all_sz, ground_sz


def hermitian_test(matrix, sz_list):
    """
    Check if matrix is Hermitian. If not, "ValueError".
    :param: matrix, sz_list
    """
    for i in range(0, len(matrix)):
        for j in range(i, len(matrix)):
            element_1 = np.round(matrix[i, j], 4)
            element_2 = np.round(np.conjugate(matrix[j, i]), 4)
            if element_1 != element_2:
                state_1 = i // len(sz_list)
                state_2 = j // len(sz_list)

                print('State 1:', state_1, ', State 2:', state_2, ', row:', i, ', column:', j,
                      ', value:', matrix[i, j])
                print('State 2:', state_2, ', State 1:', state_2, ', row:', j, ', column:', i,
                      ', value:', matrix[j, i])
                raise ValueError("Matrix is not Hermitian: see the elements shown above (SOCs in cm-1)")


def reordering_eigenvectors(eigenval, eigenvect):
    """
    Reorder eigenvectors and eigenenergies by eigenvectors weight coefficients.
    :param: eigenvalues, eigenvectors
    :return: eigenvalues, eigenvectors
    """
    change_order = np.zeros(len(eigenvect), dtype=complex)
    for v_1 in range(0, len(eigenvect)):
        for v_2 in range(v_1, len(eigenvect)):

            if abs(eigenvect[v_1, v_2]) > abs(eigenvect[v_1, v_1]):
                change_order[:] = eigenvect[:, v_1]
                eigenvect[:, v_1] = eigenvect[:, v_2]
                eigenvect[:, v_2] = change_order[:]

                change_order.real[0] = eigenval[v_1]
                eigenval[v_1] = eigenval[v_2]
                eigenval[v_2] = change_order.real[0]
    return eigenval, eigenvect


def diagonalization(matrix):
    """
    Diagonalize Hamiltonian. Eigenvectors-eigenvalues are ordered by weight coefficients. Construct the diagonal matrix.
    :param: matrix
    :return: eigenvalues, eigenvectors, diagonal_matrix
    """
    eigenvalues, eigenvectors = linalg.eigh(matrix)
    eigenvalues, eigenvectors = reordering_eigenvectors(eigenvalues, eigenvectors)

    rotation_inverse = np.linalg.inv(eigenvectors)
    diagonal_matrix = np.matmul(np.matmul(rotation_inverse, matrix), eigenvectors)
    return eigenvalues, eigenvectors, diagonal_matrix


def angular_matrices_obtention(hamiltonian, input_angular_matrix, sz_list):
    """
    Angular matrix from non-relativistic states (states from Q-Chem) is expanded in the relativistic states. In
    < B(S,Sz) | Sx | A(S',Sz') >, < B(S,Sz)| corresponds to rows and | A(S',Sz') > to columns.
    :param: eigenvectors, input_angular_matrix, sz_list
    :return: angular_matrix
    """

    eigenvectors = diagonalization(hamiltonian)[1]

    angular_matrix = np.zeros((len(sz_list), len(sz_list), 3), dtype=complex)

    for k in range(0, 3):
        for row in range(0, len(sz_list)):  # state i
            for column in range(0, len(sz_list)):  # state j

                for bra in range(0, len(eigenvectors)):  # state <B|
                    for ket in range(0, len(eigenvectors)):  # state |A>

                        coeff_bra = np.conj(eigenvectors[bra, row])
                        coeff_ket = (eigenvectors[ket, column])
                        angular_value = input_angular_matrix[bra, ket, k]

                        element = coeff_bra * coeff_ket * angular_value
                        angular_matrix[row, column, k] += element
    return angular_matrix


class GTensor:
    def __init__(self,
                 parsed_data,
                 selected_states=None):
        """
        Compute g-tensor for RAS-CI states

        :param: parsed_data: parsed data dictionary from a RAS-CI calculation
        :param: selected_states: custom selected states (if None use all)
        """

        self._parsed_data = parsed_data
        if selected_states is None:
            selected_states = [i + 1 for i in range(len(parsed_data['excited_states']))]

        self._selected_states = selected_states

    def get_spin_orbit_coupling(self):
        """
        Get: i) Spin-orbit matrix with dimensions 'bra' x 'ket', with spin order (-Ms , +Ms) in the order of
        "selected_states", ii) [-sz,..,+sz] of the highest s2 state, iii) [-sz,..,+sz] of ground state

        :return: selected spin-orbit couplings
        """

        def get_all_socs(line, n_states, multiplicities, all_sz, soc_selection):
            """
            Get SOC matrix with spin order (-Ms , +Ms) in the order of "selected_states". If there is no value in that
            |I,S,Ms> state, SOC = 0

            :param: line, n_states, multiplicities, all_sz, soc_selection
            :return: all_soc
            """
            all_soc = np.zeros((n_states * len(all_sz), n_states * len(all_sz)), dtype=complex)

            for i in range(0, n_states):
                for j in range(0, n_states):
                    if i != j:
                        i_multip = -0.5 * (-1 + np.sqrt(1 + 4 * multiplicities[i]))
                        i_multip = 0.5 * np.round(i_multip / 0.5)
                        # Making multiplicity numbers multiples of doublets (s=0.5)
                        j_multip = -0.5 * (-1 + np.sqrt(1 + 4 * multiplicities[j]))
                        j_multip = 0.5 * np.round(j_multip / 0.5)
                        # Making multiplicity numbers multiples of doublets (s=0.5)

                        i_index = all_sz.index(i_multip)
                        j_index = all_sz.index(j_multip)

                        i_position = i_index + len(all_sz) * i
                        j_position = j_index + len(all_sz) * j

                        i_j_soc_matrix = line[(i + 1, j + 1)][soc_selection]

                        for sz_1 in range(0, len(i_j_soc_matrix)):
                            for sz_2 in range(0, len(i_j_soc_matrix[0])):
                                all_soc[j_position + sz_1, i_position + sz_2] = i_j_soc_matrix[sz_1][sz_2]
            return all_soc

        def get_selected_states_socs(n_states, all_sz, socs):
            """
            Get SOC matrix of the selected states "n_states" from all the socs "all_sz"

            :param: n_states, all_sz, socs
            :return: selected_soc
            """
            selected_soc = np.zeros((len(n_states) * len(all_sz), len(n_states) * len(all_sz)), dtype=complex)
            for i, all_i in enumerate(n_states):
                for j, all_j in enumerate(n_states):

                    for sz_1 in range(0, len(all_sz)):
                        for sz_2 in range(0, len(all_sz)):
                            i_index = i * len(all_sz) + sz_1
                            j_index = j * len(all_sz) + sz_2
                            all_i_index = (all_i - 1) * len(all_sz) + sz_1
                            all_j_index = (all_j - 1) * len(all_sz) + sz_2

                            selected_soc[i_index][j_index] = socs[all_i_index][all_j_index]
            return selected_soc

        data_interstate = self._parsed_data['interstate_properties']
        total_states = len(self._parsed_data['excited_states'])

        all_multiplicities, sz_max_allstates, sz_ground_state = get_states_sz(self._parsed_data, self._selected_states)
        all_socs = get_all_socs(data_interstate, total_states, all_multiplicities, sz_max_allstates, 'total_soc_mat')
        selected_socs = get_selected_states_socs(self._selected_states, sz_max_allstates, all_socs)
        selected_socs = selected_socs / (constants.physical_constants['hartree-inverse meter relationship'][0] / 100)

        return selected_socs

    def get_hamiltonian(self):
        """
        Construct Hamiltonian matrix with dimensions 'bra' x 'ket', with spin order (-Ms , +Ms) in the order of
        "selected_states". Make hermitian test to the matrix

        :return: hamiltonian
        """

        eigen_energies, excitation_energies_ras = get_eigenenergies(self._selected_states, self._parsed_data)  # extraction energies
        spin_orbit_couplings = self.get_spin_orbit_coupling()

        sz_list = get_states_sz(self._parsed_data, self._selected_states)[1]

        hamiltonian = np.zeros((len(self._selected_states) * len(sz_list),
                                len(self._selected_states) * len(sz_list)), dtype=complex)

        for i in range(0, len(self._selected_states) * len(sz_list)):
            for j in range(0, len(self._selected_states) * len(sz_list)):
                if i == j:
                    hamiltonian[i, i] = eigen_energies[i // len(sz_list)]
                else:
                    hamiltonian[i, j] = spin_orbit_couplings[i, j]
        hermitian_test(hamiltonian, sz_list)

        return hamiltonian

    def get_combination_spin_matrix(self):

        """
        Get spin matrix with dimensions ['bra' x 'ket' x 3] (x,y,z), with spin order (-Ms , +Ms) in the order of
        "selected_state". Functions "s2_to_s", "s_to_s2" and "spin_matrices" are taken from inside PyQChem

        :return: spin_matr, standard_spin_mat.
        """

        def expand_spin_matrices(s_x, s_y, s_z, max_mult, state_mult):
            """
            Expand (sxx, syy, szz) matrices with dimension of the st multiplicity to (sxx, syy, szz) with dimension of the
            maximum multiplicity of states

            :param: s_x, s_y, s_z, max_mult, state_mult
            :return: long_sx, long_sy, long_sz
            """
            long_sx = np.zeros((max_mult, max_mult), dtype=complex)
            long_sy = np.zeros((max_mult, max_mult), dtype=complex)
            long_sz = np.zeros((max_mult, max_mult), dtype=complex)

            multipl_difference = int((max_mult - state_mult) // 2)

            if multipl_difference != 0:
                for n_row in range(0, len(sx)):
                    for n_column in range(0, len(sx)):
                        iii = n_row + multipl_difference
                        jj = n_column + multipl_difference
                        long_sx[iii, jj] = s_x[n_row, n_column]
                        long_sy[iii, jj] = s_y[n_row, n_column]
                        long_sz[iii, jj] = s_z[n_row, n_column]
            else:
                long_sx = s_x
                long_sy = s_y
                long_sz = s_z
            return long_sx, long_sy, long_sz

        def form_big_spin_matrix(st, selected_state, max_mult, sxx, syy, szz, spin_matr):
            """
            Get big spin matrix with dimensions "(len(selected_state) * max_mult, len(selected_state) * max_mult, 3)"
            with s_x, s_y, s_z

            :param: st, selected_state, max_mult, sxx, syy, szz, spin_matr
            :return: spin_matr
            """
            s_dim = 0
            total_dim = selected_state.index(selected_state[st]) * max_mult
            for row in range(0, max_mult):
                for column in range(0, max_mult):
                    spin_matr[total_dim, total_dim, 0] = sxx[s_dim, s_dim]
                    spin_matr[total_dim, total_dim + column, 0] = sxx[s_dim, s_dim + column]
                    spin_matr[total_dim + row, total_dim, 0] = sxx[s_dim + row, s_dim]
                    spin_matr[total_dim + row, total_dim + column, 0] = sxx[s_dim + row, s_dim + column]

                    spin_matr[total_dim, total_dim, 1] = syy[s_dim, s_dim]
                    spin_matr[total_dim, total_dim + column, 1] = syy[s_dim, s_dim + column]
                    spin_matr[total_dim + row, total_dim, 1] = syy[s_dim + row, s_dim]
                    spin_matr[total_dim + row, total_dim + column, 1] = syy[s_dim + row, s_dim + column]

                    spin_matr[total_dim, total_dim, 2] = szz[s_dim, s_dim]
                    spin_matr[total_dim, total_dim + column, 2] = szz[s_dim, s_dim + column]
                    spin_matr[total_dim + row, total_dim, 2] = szz[s_dim + row, s_dim]
                    spin_matr[total_dim + row, total_dim + column, 2] = szz[s_dim + row, s_dim + column]
            return spin_matr

        def get_standard_spin_matrix(s2_selected_dict, max_mult, spin_mat):
            """
            Construct Standard Spin matrix from the previous Spin Matrix

            :param: s2_selected_dict, max_mult, spin_mat
            :return: standard_spin_mat
            """
            s2_ground = (float(s2_selected_dict[0]['s2']))
            ground_multip = int(2 * s2_to_s(s2_ground) + 1)
            standard_spin_mat = np.zeros((ground_multip, ground_multip, 3), dtype=complex)

            multip_difference = (max_mult - ground_multip) // 2
            for k in range(0, 3):
                for ii in range(0, ground_multip):
                    for jj in range(0, ground_multip):
                        standard_spin_mat[ii, jj, k] = spin_mat[ii + multip_difference, jj + multip_difference, k]
            return standard_spin_mat

        # Take s2 of all states and each of the selected states.
        s2_all_states = [state['multiplicity'] for i, state in enumerate(self._parsed_data['excited_states'])]
        s2_selected_dicts = [{'st': i + 1, 's2': state['multiplicity']}
                             for i, state in enumerate(self._parsed_data['excited_states']) if i + 1 in self._selected_states]

        # Take maximum multiplicity, which determines the dimension of the spin matrix (max_mult x max_mult x 3).
        max_multip = int(2 * s2_to_s(max(s2_all_states)) + 1)
        spin_matrix = np.zeros((len(self._selected_states) * max_multip, len(self._selected_states) * max_multip, 3),
                               dtype=complex)

        for state in range(0, len(self._selected_states)):
            # Form the spin matrix of the state
            s2_state = (float(s2_selected_dicts[state]['s2']))
            s = s2_to_s(s2_state)
            sx, sy, sz = spin_matrices(s)

            # Expand the spin matrix of the st to the dimension of the maximum multiplicity (max_mult)
            state_multip = int(2 * s2_to_s(s2_state) + 1)
            sx, sy, sz = expand_spin_matrices(sx, sy, sz, max_multip, state_multip)

            # Mix (sxx,syy,szz) in one spin matrix:
            spin_matrix = form_big_spin_matrix(state, self._selected_states, max_multip, sx, sy, sz, spin_matrix)

        # Take Standard Spin Matrix from Spin Matrix
        standard_spin_matrix = get_standard_spin_matrix(s2_selected_dicts, max_multip, spin_matrix)

        sz_list = get_states_sz(self._parsed_data, self._selected_states)[1]

        combination_spin_matrix = angular_matrices_obtention(self.get_hamiltonian(), spin_matrix, sz_list)

        return combination_spin_matrix, standard_spin_matrix

    def get_combination_orbital_matrix(self):
        """
        Get orbital angular momentum matrix with dimensions ['bra' x 'ket' x 3] (x,y,z), with spin order (-Ms , +Ms) in the

        :return: combination_orbital_matrix
        """

        def get_all_momentum(line, n_states):
            """
            Get Lk between all the states selected in x,y,z dimensions. | A > in columns, < B | in rows

            :param: line, selected_states
            :return: all_orbital_momentum
            """
            all_orbital_momentum = np.zeros((n_states, n_states, 3), dtype=complex)

            for i in range(0, n_states):
                for j in range(0, n_states):
                    if i != j:
                        element = line[(i + 1, j + 1)]['angular_momentum']

                        for k in range(0, 3):
                            all_orbital_momentum[j, i, k] = element[k]
            return all_orbital_momentum

        def get_selected_states_momentum(n_states, all_momentum):
            """
            Get Lk between the selected states selected in x,y,z dimensions

            :param: selected_states, all_momentum
            :return: selected_momentum
            """
            selected_momentum = np.zeros((len(n_states), len(n_states), 3), dtype=complex)

            for k in range(0, 3):
                for i, all_i in enumerate(n_states):
                    for j, all_j in enumerate(n_states):
                        selected_momentum[i][j][k] = all_momentum[all_i - 1][all_j - 1][k]
            return selected_momentum

        def get_all_multip_momentum(all_momentums, all_sz, selected_states):
            """
            Get Lk between the selected states selected in x,y,z dimensions for doublets,
            i.e. in each row (column) there are state A,-1/2 and A,+1/2

            :param: all_momentums, all_sz
            :return: lk_values
            """
            lk_values = np.zeros((len(all_momentums) * len(all_sz), len(selected_states) * len(all_sz), 3), dtype=complex)

            for k in range(0, 3):
                for i in range(0, len(all_momentums) * len(all_sz), len(all_sz)):
                    for j in range(0, len(all_momentums) * len(all_sz), len(all_sz)):

                        for multip in range(0, len(all_sz)):
                            lk_values[i + multip, j + multip][k] = all_momentums[i // len(all_sz)][j // len(all_sz)][k]
            return lk_values

        sz_list = get_states_sz(self._parsed_data, self._selected_states)[1]

        total_states = len(self._parsed_data['excited_states'])
        all_lk = get_all_momentum(self._parsed_data['interstate_properties'], total_states)
        selected_lk = get_selected_states_momentum(self._selected_states, all_lk)
        orbital_matrix = get_all_multip_momentum(selected_lk, sz_list, self._selected_states)

        return angular_matrices_obtention(self.get_hamiltonian(), orbital_matrix, sz_list)

    def get_number_of_states(self):
        """
        get number of states used to compute the g-tensor

        :return: number of states
        """
        return len(self._selected_states)

    def get_selected_states(self):
        """
        get selected states used to compute the g-tensor

        :return:
        """
        return self._selected_states

    def get_g_tensor(self):
        """
        get diagonal values of g-tensor

        :return: diagonal values of g-tensor
        """

        def j_matrix_formation(spin, orbital, list_sz, sz_ground):
            """
            Get the total angular momentum matrix from the orbital and spin angular momentums. Then, expand it to the
            multiplicity of the ground state multiplicity "sz_ground".
            :param: spin, orbital, list_sz, sz_ground
            :return: j_matr
            """
            j_big_matrix = lande_factor * spin + orbital

            sz_difference = (len(list_sz) - len(sz_ground)) // 2
            j_matr = np.zeros((len(sz_ground), len(sz_ground), 3), dtype=complex)
            for k in range(0, 3):
                for ii in range(0, len(j_matr)):
                    for j in range(0, len(j_matr)):
                        j_matr[ii, j, k] = j_big_matrix[ii + sz_difference, j + sz_difference, k]

                hermitian_test(j_matr[:, :, k], list_sz)
            return j_matr

        def j_diagonalization(matrix):
            """
            J-matrix diagonalization giving i) eigenvectors ii) diagonal matrix

            :param: matrix
            :return: eigenvalues, eigenvectors, diagonal_matrix
            """
            eigenvalues, eigenvectors = linalg.eigh(matrix)
            rotation_inverse = np.linalg.inv(eigenvectors)
            diagonal_matrix = np.matmul(np.matmul(rotation_inverse, matrix), eigenvectors)
            return eigenvectors, diagonal_matrix

        def trace_g_values(total_j, spin_matr):
            """
            Obtaining g_value with projection of the traces

            :param: total_j, spin_matr
            :return: g_value
            """
            a = np.matmul(total_j, spin_matr)
            b = np.matmul(spin_matr, spin_matr)
            g_value = np.trace(a) / np.trace(b)
            return g_value

        combination_spin_matrix, standard_spin_matrix = self.get_combination_spin_matrix()
        combination_orbital_matrix = self.get_combination_orbital_matrix()
        _, sz_list, sz_ground = get_states_sz(self._parsed_data, self._selected_states)

        j_matrix = j_matrix_formation(combination_spin_matrix,
                                      combination_orbital_matrix,
                                      sz_list, sz_ground)

        # PROJECTION TECHNIQUE TO OBTAIN THE TRIANGULAR G-MATRIX:
        # 1) g-value zz
        j_matrix_rotation, j_matrix_diagonal_z = j_diagonalization(j_matrix[:, :, 2])
        g_matrix_triangular = np.zeros((3, 3), dtype=complex)
        g_matrix_triangular[2, 2] = trace_g_values(j_matrix_diagonal_z, standard_spin_matrix[:, :, 2])

        # 2) pseudospin z
        pseudospin_matrix = np.zeros((len(j_matrix[0, :]), len(j_matrix[:, 0]), 3), dtype=complex)
        pseudospin_matrix[:, :, 2] = j_matrix[:, :, 2] / g_matrix_triangular[2, 2]

        # 3) g-value yz
        g_matrix_triangular[1, 2] = trace_g_values(j_matrix[:, :, 1], pseudospin_matrix[:, :, 2])
        residue = j_matrix[:, :, 1] - g_matrix_triangular[1, 2] * pseudospin_matrix[:, :, 2]

        # 4) g-value yy
        j_matrix_rotation, j_matrix_transformed_y = j_diagonalization(j_matrix[:, :, 1])
        g_matrix_triangular[1, 1] = trace_g_values(j_matrix_transformed_y, standard_spin_matrix[:, :, 2])

        # 5) pseudospin y
        pseudospin_matrix[:, :, 1] = residue[:, :] / g_matrix_triangular[1, 1]

        # 6) g-value xy and xz
        g_matrix_triangular[0, 1] = trace_g_values(j_matrix[:, :, 0], pseudospin_matrix[:, :, 1])
        g_matrix_triangular[0, 2] = trace_g_values(j_matrix[:, :, 0], pseudospin_matrix[:, :, 2])

        # 7) g-value xx
        j_matrix_rotation, j_matrix_transformed_x = j_diagonalization(j_matrix[:, :, 0])
        g_matrix_triangular[0, 0] = trace_g_values(j_matrix_transformed_x, standard_spin_matrix[:, :, 2])

        # 8) from g_matrix_triangular to g-shifts
        upper_g_matrix = np.matmul(g_matrix_triangular, np.transpose(g_matrix_triangular))
        upper_g_matrix_eigenvalues, rotation_g_matrix, upper_g_matrix_diagonal = diagonalization(upper_g_matrix)

        tensor_g = [sqrt(upper_g_matrix_diagonal[i, i]) for i in range(3)]

        return tensor_g

    def get_g_shift(self, use_ppm=False):
        """
        get g-tensor shifts in ppt

        :param use_ppm: if True return shifts in ppm
        :return: shifts list
        """

        tensor_g = self.get_g_tensor()
        g_shift = (np.array(tensor_g) - lande_factor) * 1000  # ppt
        if use_ppm:
            g_shift = g_shift * 1000

        return g_shift

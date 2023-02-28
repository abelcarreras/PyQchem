import numpy as np
from pyqchem.tools import rotate_coordinates
from scipy.optimize import minimize
from itertools import product
from pyqchem.structure import atom_data
from pyqchem.units import ANGSTROM_TO_AU, AMU_TO_ELECTRONMASS, BOHR_TO_CM, SPEEDOFLIGHT_AU, AU_TO_EV, KB_EV
import warnings


def do_reflection(coordinates, axis):
    coordinates = np.array(coordinates)
    coordinates[:, axis] = -coordinates[:, axis]
    return coordinates.tolist()


def get_principal_axis_and_moments_of_inertia(cm_vectors, masses):
    """
    Compute the principal axis and moments of inertia
    (compact and vectorized)

    :param cm_vectors: list of distance vectors of atoms respect center of mass
    :param masses: list of atomic masses
    :return: moments of inertia, principal axis of inertia in rows
    """

    # Build inertia tensor
    inertia_tensor = np.zeros((3, 3))
    for m, c in zip(masses, cm_vectors):
        inertia_tensor += m * (np.identity(3) * np.dot(c, c) - np.outer(c, c))

    # Compute eigenvalues and eigenvectors of inertia tensor
    e_values, e_vectors = np.linalg.eigh(inertia_tensor)  # be careful eigenvectors in columns!
    axis_of_inertia = e_vectors.T
    moments_of_inertia = e_values

    return moments_of_inertia, axis_of_inertia


def get_reduced_mass(structure, modes, mass_weighted=True):
    """
    compute the reduced mass of normal modes (AMU)

    :param structure: Structure object
    :param modes: normal modes
    :param mass_weighted: if True, Normal modes are assumed mass weighted, else False
    :return: reduced mass in AMU
    """

    mass_vector = []
    for m in structure.get_atomic_masses():
        for _ in range(3):
            mass_vector.append(m)

    r_mass = []
    for mode in modes:
        mode = np.array(mode).flatten()

        if mass_weighted:
            mode = mode * np.sqrt(mass_vector)
            mode = mode / np.linalg.norm(mode)

        mode = mode / np.sqrt(mass_vector)
        r_mass.append(1 / np.dot(mode, mode))

    return r_mass


class NormalModes:
    def __init__(self, structure, modes, frequencies):
        self._coordinates = np.array(structure.get_coordinates())
        self._modes = modes  # in mass weighted coordinates
        self._symbols = structure.get_symbols()
        self._frequencies = frequencies  # cm-1
        self._reduced_mass = get_reduced_mass(structure, modes)  # AMU

        is_freq_negative = np.array(self._frequencies) < 0
        if is_freq_negative.any():
            warnings.warn('Some modes have negative frequency')

    def __len__(self):
        return len(self._modes)  # mass weighted coordinates

    def get_coordinates(self):
        return self._coordinates  # angstrom

    def get_frequencies(self):
        return self._frequencies  # cm-1

    def get_reduced_mass(self):
        return self._reduced_mass  # AMU

    def get_atomic_masses(self):
        atomic_numbers = [[data[1].upper() for data in atom_data].index(element.upper()) for element in self._symbols]
        return [atom_data[an][3] for an in atomic_numbers]

    def get_displacements(self):

        mass_atoms = np.array([[m] * 3 for m in self.get_atomic_masses()]).flatten()
        modes = np.array([np.array(mode).flatten().tolist() for mode in self._modes])

        rot_modes = []
        for mode, rm in zip(modes, self._reduced_mass):
            m_b = np.sqrt(rm / np.array(mass_atoms))
            rot_modes.append(mode/m_b)

        return np.array(rot_modes).T

    def apply_rotation_angles(self, angles):
        self._coordinates = rotate_coordinates(self._coordinates, angles[0], [1, 0, 0])
        self._coordinates = rotate_coordinates(self._coordinates, angles[1], [0, 1, 0])
        self._coordinates = rotate_coordinates(self._coordinates, angles[2], [0, 0, 1])

        rot_modes = []
        for mode in self._modes:
            mode_temp = rotate_coordinates(mode,      angles[0], [1, 0, 0])
            mode_temp = rotate_coordinates(mode_temp, angles[1], [0, 1, 0])
            mode_temp = rotate_coordinates(mode_temp, angles[2], [0, 0, 1])
            rot_modes.append(mode_temp)

        self._modes = rot_modes

    def apply_rotation(self, angle, axis):
        self._coordinates = rotate_coordinates(self._coordinates, angle, axis)

        rot_modes = []
        for mode in self._modes:
            mode_temp = rotate_coordinates(mode, angle, axis)
            rot_modes.append(mode_temp)

        self._modes = rot_modes

    def apply_reflection(self, axis):
        self._coordinates = np.array(self._coordinates)
        self._coordinates[:, axis] = -self._coordinates[:, axis]
        self._coordinates = self._coordinates.tolist()

        rot_modes = []
        for mode in self._modes:
            mode_new = np.array(mode)
            mode_new[:, axis] = -mode_new[:, axis]
            rot_modes.append(mode_new.tolist())
        self._modes = rot_modes

    def number_of_modes(self):
        return len(self._modes)

    def center_of_mass(self):
        coor = self.get_coordinates()
        mass = self.get_atomic_masses()
        return np.average(coor, axis=0, weights=mass)

    def set_center_of_mass(self):
        coor = self.get_coordinates()
        cm = self.center_of_mass()
        self._coordinates = list(np.array(coor) - cm)

    def align_axis_of_inertia(self, vector_1, vector_2):

        vector_1 = np.array(vector_1)/np.linalg.norm(vector_1)
        vector_2 = np.array(vector_2)/np.linalg.norm(vector_2)

        self.set_center_of_mass()

        mass = self.get_atomic_masses()
        moments, axis = get_principal_axis_and_moments_of_inertia(self._coordinates, mass)

        angle_x = np.arccos(np.dot(axis[2], vector_1)/np.linalg.norm(axis[2]))
        axis_px = np.cross(axis[2], vector_1)

        if not (np.linalg.norm(axis_px) == 0):
            self.apply_rotation(-angle_x, axis_px)

        moments, axis = get_principal_axis_and_moments_of_inertia(self._coordinates, mass)

        angle_y = np.arccos(np.dot(axis[1], vector_2)/np.linalg.norm(axis[1]))
        self.apply_rotation(-angle_y, vector_1)

    def trim_modes(self, modes_list):
        modes = []
        frequencies = []
        reduced_mass = []
        for i in modes_list:
            modes.append(self._modes[i])
            frequencies.append(self._frequencies[i])
            reduced_mass.append(self._reduced_mass[i])

        self._modes = modes
        self._frequencies = frequencies
        self._reduced_mass = reduced_mass

    def trim_negative_frequency_modes(self):

        trim_list = []
        for i, f in enumerate(self._frequencies):
            if f > 0:
                trim_list.append(i)
        self.trim_modes(trim_list)

class VibrationalState:
    def __init__(self, q_index=0, vector_rep=(), frequencies=()):
        """
        Define a vibrational state

        :param q_index: index of electronic state associated to this vibrational state
        :param vector_rep: vector representation of vibrational state
        :param frequencies: list of frequencies associated to the vector_rep

        """
        self.q_index = q_index
        self.vector_rep = vector_rep
        self.frequencies = frequencies

    def __hash__(self):
        return hash(tuple(self.vector_rep))

    def get_label(self):
        """
        get label of the vibrational state

        :return: label string
        """
        label = '{}('.format(self.q_index)
        c = 0
        for i, v in enumerate(self.vector_rep):
            if v != 0:
                if c > 0:
                    label += ','

                label += '{}v{}'.format(v, i)
                c += 1
        if np.sum(self.vector_rep) == 0:
            label += '0'

        label += ')'
        return label

    def total_quanta(self):
        return np.sum(self.vector_rep)

    def get_vib_energy(self):
        return np.sum(self.frequencies * self.vector_rep) * AU_TO_EV


class VibrationalTransition:
    def __init__(self, origin, target, fcf, excitation_energy=0.0, reorganization_energy=0.0):
        """
        Define a transition between two vibrational states

        :param origin: origin VibrationalState
        :param target: target VibrationalState
        :param fcf: Franck-Condon factor
        :param intensity: intensity of spectrum peak
        :param excitation_energy: electronic excitation energy [eV]
        :param reorganization_energy: total reorganization energy [eV]
        """
        self.excitation_energy = excitation_energy
        self.reorganization_energy = reorganization_energy
        self.origin = origin
        self.target = target
        self.fcf = fcf
        #self.intensity = fcf*fcf

        #self.intensity = fcf*fcf * np.exp(-self.origin.get_vib_energy() / (600 * KB_EV))

    @property
    def energy(self):
        """
        energy associated to the vibrational transition (vib + electronic) [eV]

        :return: energy
        """
        return self.excitation_energy + (self.target.get_vib_energy() - self.origin.get_vib_energy())

    @property
    def energy_absorption(self):
        """
        energy associated to the vibrational transition (vib + electronic) [eV]
        corrected with reorganization energy for absorption spectrum

        :return: energy
        """
        delta_vib_energy = self.target.get_vib_energy() - self.origin.get_vib_energy()
        return self.excitation_energy + delta_vib_energy + self.reorganization_energy/2

    @property
    def energy_emission(self):
        """
        energy associated to the vibrational transition (vib + electronic) [eV]
        corrected with reorganization energy for emission spectrum

        :return: energy
        """
        delta_vib_energy = self.target.get_vib_energy() - self.origin.get_vib_energy()
        return self.excitation_energy + delta_vib_energy - self.reorganization_energy/2

    def get_label(self, sign=None):
        """
        get label of the vibrational transition

        :return: label string
        """
        if sign == 0:
            return '{} <- {}'.format(self.origin.get_label(), self.target.get_label())
        elif sign == 1:
            return '{} -> {}'.format(self.origin.get_label(), self.target.get_label())
        else:
            return '{} <-> {}'.format(self.origin.get_label(), self.target.get_label())

    def get_intensity_absorption(self, temperature=300):
        """
        get transition absorption intensity as a function of the temperature

        :param temperature: the temperature in Kelvin
        :return: intensity
        """
        return self.fcf**2 * np.exp(-self.origin.get_vib_energy() / (temperature * KB_EV))

    def get_intensity_emission(self, temperature=300):
        """
        get transition absorption intensity as a function of the temperature

        :param temperature: the temperature in Kelvin
        :return: intensity
        """
        return self.fcf**2 * np.exp(-self.target.get_vib_energy() / (temperature * KB_EV))


class Duschinsky:
    def __init__(self,
                 structure_initial,  # angstrom
                 structure_final,  # angstrom
                 modes_initial,  # mass weighted coordinates
                 modes_final,  # mass weighted coordinates
                 frequencies_initial,  # cm -1
                 frequencies_final,  # cm -1
                 n_max_modes=None,
                 ):

        self._modes_initial = NormalModes(structure_initial,
                                          modes_initial,
                                          frequencies_initial)

        self._modes_final = NormalModes(structure_final,
                                        modes_final,
                                        frequencies_final)

        self._modes_initial.trim_negative_frequency_modes()
        self._modes_final.trim_negative_frequency_modes()

        # make sure duschinsky matrix will be square after trim negative
        if len(self._modes_initial) != len(self._modes_final):
            n = np.min([len(self._modes_initial), len(self._modes_final)])
            self._modes_initial.trim_modes(list(range(n)))
            self._modes_final.trim_modes(list(range(n)))

        # Define restrictions of necessary
        self._modes_origin_indices, self._modes_target_indices = self._get_restricted_modes(n_max_modes)

    def _get_restricted_modes(self, n_max_modes):
        """
        compute the list of restricted modes to use for origin and target
        :param n_max_modes: max nodes to use
        :return:
        """

        if n_max_modes is None:
            n_modes_total = len(self._modes_final)
            return list(range(n_modes_total)), list(range(n_modes_total))

        s = self.get_s_matrix()
        d_ini, d_fin = self.get_d_vector()

        if n_max_modes > len(d_fin):
            warnings.warn('n_max_modes ({}) > n_modes ({}).'.format(n_max_modes, len(d_fin)))
            n_max_modes = len(d_fin)

        q_vector = np.zeros_like(d_fin)
        for i in range(n_max_modes):
            q_vector[i] = 1

        sdot_initial = np.dot(s, q_vector) + d_fin
        sdot_final = np.dot(s.T, q_vector - d_fin)

        indices_origin = np.abs(sdot_final).argsort()[::-1]
        indices_target = np.abs(sdot_initial).argsort()[::-1]

        modes_origin_indices = set()
        modes_target_indices = set()

        for i in range(n_max_modes):
            modes_origin_indices.add(i)
            modes_origin_indices.add(indices_origin[i])

            modes_target_indices.add(i)
            modes_target_indices.add(indices_target[i])

        return modes_origin_indices, modes_target_indices

    def get_restricted_modes(self):
        """
        return restricted modes used for origin and target
        :return: lists of normal modes indices for origin and target
        """
        return self._modes_origin_indices, self._modes_target_indices

    def align_coordinates_pmi(self):
        """
        Align principal axis of inertia of molecules along x,y,z axis

        :return:
        """

        vector_1 = [0, 0, 1]
        vector_2 = [0, 1, 0]

        self._modes_initial.set_center_of_mass()
        self._modes_final.set_center_of_mass()

        self._modes_initial.align_axis_of_inertia(vector_1, vector_2)
        self._modes_final.align_axis_of_inertia(vector_1, vector_2)


    def align_coordinates_rsm(self):
        """
        Align molecules by minimizing RSM between the two geometries

        :return:
        """

        self._modes_final.set_center_of_mass()
        self._modes_initial.set_center_of_mass()

        c_initial = self._modes_initial.get_coordinates()
        c_final = self._modes_final.get_coordinates()

        def optimize_function_1(angle):
            coor = rotate_coordinates(c_final, angle[0], [1, 0, 0])
            coor = rotate_coordinates(coor,    angle[1], [0, 1, 0])
            coor = rotate_coordinates(coor,    angle[2], [0, 0, 1])

            return np.sum(np.linalg.norm(np.subtract(c_initial, coor), axis=1))

        def optimize_function_2(angle):
            coor = do_reflection(c_final, 0)
            coor = rotate_coordinates(coor, angle[0], [1, 0, 0])
            coor = rotate_coordinates(coor, angle[1], [0, 1, 0])
            coor = rotate_coordinates(coor, angle[2], [0, 0, 1])

            return np.sum(np.linalg.norm(np.subtract(c_initial, coor), axis=1))

        res_1 = minimize(optimize_function_1, np.array([0, 0, 0]), method='Nelder-Mead', tol=1e-6)
        res_2 = minimize(optimize_function_2, np.array([0, 0, 0]), method='Nelder-Mead', tol=1e-6)
        if res_1.fun < res_2.fun:
            self._modes_final.apply_rotation_angles(res_1.x)
        else:
            self._modes_final.apply_reflection(0)
            self._modes_final.apply_rotation_angles(res_2.x)

    def get_s_matrix(self):
        """
        compute Duschinsky rotation matrix S

        :return: the matrix
        """
        l_i = np.array(self._modes_initial.get_displacements())
        l_f = np.array(self._modes_final.get_displacements())

        s = np.dot(l_f.T, l_i)

        return s

    def get_d_vector(self):
        """
        compute displacement vector of displacements

        :return: the displacement vector in A.U.
        """

        mass_vector = np.array([[m] * 3 for m in self._modes_initial.get_atomic_masses()]).flatten() * AMU_TO_ELECTRONMASS
        coor_initial = self._modes_initial.get_coordinates()
        coor_final = self._modes_final.get_coordinates()

        l_i = np.array(self._modes_initial.get_displacements())
        l_f = np.array(self._modes_final.get_displacements())

        t = np.diag(mass_vector)

        diff = np.subtract(coor_initial, coor_final).flatten() * ANGSTROM_TO_AU
        d_fin = np.dot(np.dot(l_f.T, np.sqrt(t)), diff)
        d_ini = np.dot(np.dot(l_i.T, np.sqrt(t)), diff)

        return d_ini, d_fin

    def _get_lambda_matrices(self):

        initial_freq, final_freq = self._get_frequencies()

        freq_vector_ini = np.diag(np.sqrt(initial_freq))
        freq_vector_fin = np.diag(np.sqrt(final_freq))

        return freq_vector_ini, freq_vector_fin

    def get_j_matrix(self):
        """
        J = Lmb_f * S * Lmb_i^(-1)

        :return: J matrix
        """
        lmb_ini, lmb_fin = self._get_lambda_matrices()
        s = self.get_s_matrix()

        return np.dot(np.dot(lmb_fin, s), np.linalg.inv(lmb_ini))

    def get_q_matrix(self):
        """
        Q = (1 + J^T * J)^(-1)

        :return: Q matrix
        """
        j_matrix = self.get_j_matrix()
        return np.linalg.inv(np.identity(len(j_matrix)) + np.dot(j_matrix.T, j_matrix))

    def get_p_matrix(self):
        """
        P = J * Q * J^T

        :return: P matrix
        """
        j_matrix = self.get_j_matrix()
        q_matrix = self.get_q_matrix()

        return np.dot(np.dot(j_matrix, q_matrix), j_matrix.T)

    def get_r_matrix(self):
        """
        R = Q * J^T

        :return: R matrix
        """
        j_matrix = self.get_j_matrix()
        q_matrix = self.get_q_matrix()

        return np.dot(q_matrix, j_matrix.T)

    def get_dt_vector(self):
        """
        dt = Lmb_f * d

        :return: dt vector
        """
        lmb_ini, lmb_fin = self._get_lambda_matrices()
        d_ini, d_fin = self.get_d_vector()

        return np.dot(lmb_fin, d_fin)

    def _get_frequencies(self):
        """
        angular frequencies in A.U.

        :return:
        """

        units_factor = BOHR_TO_CM * 2 * np.pi * SPEEDOFLIGHT_AU  # from ordinary frequency in cm-1 to angular in A.U.

        i_freq = np.array(self._modes_initial.get_frequencies()) * units_factor
        f_freq = np.array(self._modes_final.get_frequencies()) * units_factor

        return i_freq, f_freq

    def get_transitions(self,
                        max_vib_origin=2,
                        max_vib_target=2,
                        excitation_energy=0.0,
                        reorganization_energy=0.0,
                        add_hot_bands=False,
                        q_index_origin=0,
                        q_index_target=1
                        ):

        freq_origin, freq_target = self._get_frequencies()
        s = self.get_s_matrix()
        assert np.abs(np.linalg.det(s)) > 1e-2

        q = self.get_q_matrix()
        p = self.get_p_matrix()
        r = self.get_r_matrix()
        dt = self.get_dt_vector()
        n_modes = len(s)

        # compute FCF <0|0>
        freq_rel_prod = np.product(np.divide(freq_origin, freq_target))
        exponential = np.exp(-1 / 2 * np.dot(dt, np.dot(np.identity(n_modes) - p, dt)))
        pre_factor = 2**(n_modes / 2) / np.sqrt(abs(np.linalg.det(s)))
        fcf_00 = pre_factor * freq_rel_prod ** (-1 / 4) * np.sqrt(np.linalg.det(q)) * exponential

        # evaluate pre-factor for FCF matrices
        ompd = np.sqrt(2) * np.dot(np.identity(n_modes) - p, dt)
        tpmo = 2 * self.get_p_matrix() - np.identity(len(self.get_p_matrix()))
        tqmo = 2 * self.get_q_matrix() - np.identity(len(self.get_p_matrix()))
        tr = 2 * self.get_r_matrix()
        rd = np.sqrt(2) * np.dot(r, dt)

        # cache function to accelarate FCF
        def cache_fcf(func):
            cache = {}
            def wrapper(v1, k1, v2, k2):
                key = (tuple(v1), k1, tuple(v2), k2)
                if key not in cache:
                    cache[key] = func(v1, k1, v2, k2)
                return cache[key]

            return wrapper

        @cache_fcf
        def evalSingleFCFpy(origin_vector, k_origin, target_vector, k_target):

            if np.sum(origin_vector) == 0 and np.sum(target_vector) == 0:
                return fcf_00

            if k_origin == 0:
                ksi = 0
                while target_vector[ksi] == 0:
                    ksi = ksi + 1
                target_vector[ksi] -= 1

                fcf = ompd[ksi] * evalSingleFCFpy(origin_vector, 0, target_vector, k_target - 1)

                if k_target > 1:
                    for theta in range(ksi, n_modes):
                        if target_vector[theta] > 0:
                            tmp_dbl = tpmo[ksi, theta] * np.sqrt(target_vector[theta])
                            target_vector[theta] -= 1
                            tmp_dbl *= evalSingleFCFpy(origin_vector, 0, target_vector, k_target - 2)
                            fcf += tmp_dbl
                            target_vector[theta] += 1
                fcf /= np.sqrt(target_vector[ksi] + 1)
                target_vector[ksi] += 1

            else:
                ksi = 0
                while origin_vector[ksi] == 0:
                    ksi = ksi + 1

                origin_vector[ksi] -= 1
                fcf = -rd[ksi] * evalSingleFCFpy(origin_vector, k_origin - 1, target_vector, k_target)

                for theta in range(ksi, n_modes):
                    if origin_vector[theta] > 0:
                        tmp_dbl = tqmo[ksi, theta] * np.sqrt(origin_vector[theta])
                        origin_vector[theta] -= 1
                        tmp_dbl *= evalSingleFCFpy(origin_vector, k_origin - 2, target_vector, k_target)
                        fcf += tmp_dbl
                        origin_vector[theta] += 1

                if k_target > 0:
                    for theta in range(0, n_modes):
                        if target_vector[theta] > 0:
                            tmp_dbl = tr[ksi, theta] * np.sqrt(target_vector[theta])
                            target_vector[theta] -= 1
                            tmp_dbl *= evalSingleFCFpy(origin_vector, k_origin - 1, target_vector, k_target - 1)
                            fcf += tmp_dbl
                            target_vector[theta] += 1

                fcf /= np.sqrt(origin_vector[ksi] + 1)
                origin_vector[ksi] += 1

            return fcf

        def generate_configurations(n, restrict_list=None):
            """
            generate vector configuration of lengh n_modes and filled with total quanta == n
            :param n: total quanta
            :param restrict_list: list of allowed vibrational modes
            :return:
            """

            def transform(conf, n_modes):
                new_conf = np.zeros((n_modes), dtype=int)
                for i, j in enumerate(restrict_list):
                    new_conf[j] = conf[i]
                return new_conf

            if restrict_list is None:
                for conf in product(range(0, n + 1), repeat=n_modes):
                    if np.sum(conf) == n:
                        yield conf
            else:
                for conf in product(range(0, n + 1), repeat=len(restrict_list)):
                    if np.sum(conf) == n:
                        yield transform(conf, n_modes)

        def get_state_list(n, q_index, frequencies=(), restrict_list=None):
            state_list = []
            for conf in generate_configurations(n, restrict_list):
                state_list.append(VibrationalState(q_index=q_index,
                                                   vector_rep=list(conf),
                                                   frequencies=frequencies))

            return state_list

        s0_origin = VibrationalState(q_index=q_index_origin, vector_rep=[0] * n_modes, frequencies=freq_origin)
        s0_target = VibrationalState(q_index=q_index_target, vector_rep=[0] * n_modes, frequencies=freq_target)

        # (0)->(0) transition
        transition_list = [VibrationalTransition(s0_origin, s0_target,
                                                 fcf=fcf_00,
                                                 excitation_energy=excitation_energy,
                                                 reorganization_energy=reorganization_energy)]


        # (0)<->(n) transitions
        for i in range(0, max_vib_target):
            state_list = get_state_list(i + 1, q_index=q_index_target, frequencies=freq_target, restrict_list=self._modes_target_indices)
            for target_state in state_list:

                fcf = evalSingleFCFpy(s0_origin.vector_rep, 0, target_state.vector_rep, i+1)
                transition_list.append(VibrationalTransition(s0_origin, target_state,
                                                             fcf=fcf,
                                                             excitation_energy=excitation_energy,
                                                             reorganization_energy=reorganization_energy))

        # (n)<->(0) transitions
        for i in range(0, max_vib_origin):
            state_list = get_state_list(i + 1, q_index=q_index_origin, frequencies=freq_origin, restrict_list=self._modes_origin_indices)
            for origin_state in state_list:

                fcf = evalSingleFCFpy(origin_state.vector_rep, i+1, s0_target.vector_rep, 0)
                transition_list.append(VibrationalTransition(origin_state, s0_target,
                                                             fcf=fcf,
                                                             excitation_energy=excitation_energy,
                                                             reorganization_energy=reorganization_energy))

        if add_hot_bands:
            # (m)<->(n) transitions [Hot bands]
            state_list_origin = []
            for i in range(0, max_vib_origin):
                state_list_origin += get_state_list(i + 1, q_index=q_index_origin, frequencies=freq_origin, restrict_list=self._modes_origin_indices)

            state_list_target = []
            for i in range(0, max_vib_target):
                state_list_target += get_state_list(i + 1, q_index=q_index_target, frequencies=freq_target, restrict_list=self._modes_target_indices)

            for origin_state in state_list_origin:
                for target_state in state_list_target:
                    # print(origin_state.get_label(), ' - ', target_state.get_label())

                    k_origin = origin_state.total_quanta()
                    k_target = target_state.total_quanta()

                    fcf = evalSingleFCFpy(origin_state.vector_rep, k_origin, target_state.vector_rep, k_target)
                    transition_list.append(VibrationalTransition(origin_state, target_state,
                                                                 fcf=fcf,
                                                                 excitation_energy=excitation_energy,
                                                                 reorganization_energy=reorganization_energy))

        # sort transition list by energy
        transition_list.sort(key=lambda x: x.energy, reverse=False)

        return transition_list

    def get_huang_rys(self):
        freq_origin, freq_target = self._get_frequencies()

        d_ini, d_fin = self.get_d_vector()

        lmb_i = 0.5*freq_origin**2 * d_ini**2
        lmb_f = 0.5*freq_target**2 * d_fin**2

        s_i = lmb_i/(freq_origin)
        s_f = lmb_f/(freq_target)

        return s_i, s_f

    def get_reorganization_energies(self):
        """
        individual reorganization energies for each mode in eV

        :return:
        """
        freq_origin, freq_target = self._get_frequencies()

        d_ini, d_fin = self.get_d_vector()

        lmb_i = 0.5*freq_origin**2 * d_ini**2 * AU_TO_EV
        lmb_f = 0.5*freq_target**2 * d_fin**2 * AU_TO_EV

        return lmb_i , lmb_f


def get_duschinsky(origin_frequency_output, target_frequency_output, n_max_modes=None):
    """
    build Duschinsky instance object from frequency parser dictionary

    :param origin_frequency_output: frequency parsed output of origin state (typically ground state)
    :param target_frequency_output: frequency parsed output of target state (typically excited state)
    :param n_max_modes: restrict calculation to the first n modes
    :return: Duschinsky object
    """

    return  Duschinsky(structure_initial=origin_frequency_output['structure'],
                       structure_final=target_frequency_output['structure'],
                       modes_initial=[mode['displacement'] for mode in origin_frequency_output['modes']],
                       modes_final=[mode['displacement'] for mode in target_frequency_output['modes']],
                       frequencies_initial=[mode['frequency'] for mode in origin_frequency_output['modes']],
                       frequencies_final=[mode['frequency'] for mode in target_frequency_output['modes']],
                       n_max_modes=n_max_modes,
                       )

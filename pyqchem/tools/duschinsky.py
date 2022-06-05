import numpy as np
from pyqchem.tools import rotate_coordinates
from scipy.optimize import minimize
from pyqchem.structure import atom_data
from pyqchem.units import ANGSTROM_TO_AU, AMU_TO_ELECTRONMASS


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



class NormalModes:
    def __init__(self, coordinates, symbols, modes, reduced_mass):
        self._coordinates = np.array(coordinates)
        self._modes = modes
        self._reduced_mass = reduced_mass
        self._symbols = symbols

    def get_coordinates(self):
        return self._coordinates

    def get_atomic_masses(self):
        atomic_numbers = [[data[1].upper() for data in atom_data].index(element.upper()) for element in self._symbols]
        return [atom_data[an][3] for an in atomic_numbers]

    def get_mass_weighted_displacements(self):
        return [np.array(mode).flatten().tolist() for mode in self._modes ]

    def get_displacements(self):

        mass_atoms = np.array([[m] * 3 for m in self.get_atomic_masses()]).flatten()
        modes = self.get_mass_weighted_displacements()

        rot_modes = []
        for mode, rm in zip(modes, self._reduced_mass):
            m_b = np.sqrt(rm / np.array(mass_atoms))
            rot_modes.append(np.array(mode)/m_b)

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

        angle_x = np.arccos(np.dot(axis[0], vector_1)/np.linalg.norm(axis[0]))
        axis_px = np.cross(axis[0], vector_1)

        self.apply_rotation(-angle_x, axis_px)

        moments, axis = get_principal_axis_and_moments_of_inertia(self._coordinates, mass)

        angle_y = np.arccos(np.dot(axis[1], vector_2)/np.linalg.norm(axis[1]))
        self.apply_rotation(-angle_y, axis[0])



class Duschinsky:
    def __init__(self,
                 coordinates_initial,
                 coordinates_final,
                 modes_initial,
                 modes_final,
                 r_mass_initial,
                 r_mass_final,
                 symbols_initial,
                 symbols_final
                 ):

        self._modes_initial = NormalModes(coordinates_initial, symbols_initial, modes_initial, r_mass_initial)
        self._modes_final = NormalModes(coordinates_final, symbols_final, modes_final, r_mass_final)

    def align_coordinates(self):
        """
        Align principal axis of inertia of molecules along x,y,z axis

        :return:
        """

        vector_1 = [1, 0, 0]
        vector_2 = [0, 1, 0]

        self._modes_initial.set_center_of_mass()
        self._modes_final.set_center_of_mass()

        self._modes_initial.align_axis_of_inertia(vector_1, vector_2)
        self._modes_final.align_axis_of_inertia(vector_1, vector_2)


    def align_coordinates_2(self):
        """
        Align molecules by minimizing RSM between the two geometries

        :return:
        """

        self._modes_final.set_center_of_mass()
        self._modes_initial.set_center_of_mass()

        # self._modes_initial.apply_reflection(1)
        # self._modes_final.apply_rotation([-1.23, 0.12, -2.01])

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
        compute vector of displacements

        :return: the vector
        """

        mass_vector = np.array([[m] * 3 for m in self._modes_initial.get_atomic_masses()]).flatten()
        coor_initial = self._modes_initial.get_coordinates()
        coor_final = self._modes_final.get_coordinates()

        l_f = np.array(self._modes_final.get_displacements())
        t = np.diag(mass_vector)


        diff = np.subtract(coor_initial, coor_final).flatten()
        d = np.dot(np.dot(l_f.T, np.sqrt(t)), diff)

        # switch to atomic units
        d *= ANGSTROM_TO_AU * np.sqrt(AMU_TO_ELECTRONMASS);

        return d



def get_duschinsky(gs_structure, sg_structure, gs_output, sg_output):
    """
    build Duschinsky instance object from frequency parser dictionary

    :param gs_structure: Structure of ground state
    :param sg_structure: Structure of target state (excited state)
    :param gs_output: frequency parsed output of ground state
    :param sg_output: frequency parsed output of target state (excited state)
    :return: Duschinsky object
    """

    return Duschinsky(coordinates_initial=gs_structure.get_coordinates(),
                      coordinates_final=sg_structure.get_coordinates(),
                      modes_initial=[mode['displacement'] for mode in gs_output['modes']],
                      modes_final=[mode['displacement'] for mode in sg_output['modes']],
                      r_mass_initial=[mode['reduced_mass'] for mode in gs_output['modes']],
                      r_mass_final=[mode['reduced_mass'] for mode in sg_output['modes']],
                      symbols_initial=sg_structure.get_symbols(),
                      symbols_final=gs_structure.get_symbols()
                      )

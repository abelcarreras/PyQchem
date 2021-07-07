__author__ = 'Abel Carreras'
import numpy as np
import hashlib, json
from pyqchem.errors import StructureError


def int_to_xyz(molecule, no_dummy=True):

    internal = molecule._get_full_z_matrix()
    coordinates = [[0.0, 0.0, 0.0]]

    for line in internal[1:]:
        bi = int(line[0])  # bond index
        B = line[1]        # bond value
        ai = int(line[2])  # Angle index
        A = line[3]        # Angle value
        ci = int(line[4])  # Dihedral index
        C = line[5]        # Dihedral value

        bond = np.array(coordinates[ai-1]) - np.array(coordinates[bi-1])
        if np.linalg.norm(bond) == 0:
            bond = np.array([1, 0, 0])

        bond2 = np.array(coordinates[ci-1]) - np.array(coordinates[ai-1])
        if np.linalg.norm(bond2) == 0:
            bond2 = np.array([0, 1, 0])

        origin = bond/np.linalg.norm(bond)*B
        ref2 = bond
        ref3 = np.cross(bond, bond2)

        # Check case of linear structure
        if np.linalg.norm(ref3) == 0:
            ref3 = [0.0, 0.0, 0.1]

        inter = np.dot(rotation_matrix(ref3, np.deg2rad(A)), origin)
        final = np.dot(rotation_matrix(ref2, np.deg2rad(C)), inter)
        final = final + np.array(coordinates[bi-1])
        coordinates.append(final)

    coordinates = np.array(coordinates)

    if no_dummy:
      #  mask = np.argwhere(molecule.get_atomic_elements_with_dummy()[:,0]  == 'X')
        mask = np.argwhere((molecule.get_symbols_with_dummy()[:, 0] == 'X') |
                           (molecule.get_symbols_with_dummy()[:, 0] == 'x')).flatten()
        coordinates = np.delete(coordinates,mask,axis=0)

    return np.array(coordinates, dtype=float)


def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation around the given axis

    :param axis: rotation axis
    :param theta: rotation angle in radians

    :return: the rotation matrix
    """
    if np.dot(axis, axis) == 0.0:
        print ('Warning, reference rotation axis module is 0')
        exit()

    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/np.sqrt(np.dot(axis, axis))

    a = np.cos(theta/2)
    b, c, d = -axis*np.sin(theta/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


class Structure:
    """
    Structure object containing all the geometric data of the molecule
    """
    def __init__(self,
                 coordinates=None,
                 internal=None,
                 z_matrix=None,
                 int_label=None,
                 atom_types=None,
                 symbols=None,
                 atomic_numbers=None,
                 connectivity=None,
                 file_name=None,
                 charge=0,
                 multiplicity=1,
                 name=None,
                 int_weights=None):
        """
        :param coordinates: List containing the cartesian coordinates of each atom in Angstrom
        :param symbols: Symbols of the atoms within the molecule
        :param atomic_numbers: Atomic numbers of the atoms within the molecule
        :param charge: charge of the molecule
        :param multiplicity: multiplicity of the molecule
        """

        self._coordinates = np.array(coordinates)
        self._internal = internal
        self._z_matrix = z_matrix
        self._int_label = int_label
        self._atom_types = atom_types
        self._atomic_numbers = atomic_numbers
        self._connectivity = connectivity
        self._symbols = symbols
        self._charge = charge
        self._multiplicity = multiplicity
        self._name = name

        self._file_name = file_name
        self._int_weights = int_weights

        self._atomic_masses = None
        self._number_of_atoms = None
        self._number_of_internal = None
        self._energy = {}
        self._modes = None

        self._full_z_matrix = None

        # check input data
        if symbols is not None and coordinates is not None:
            if len(coordinates) != len(symbols):
                raise StructureError('coordinates and symbols do not match')

        if atomic_numbers is not None:
            self._symbols = [atom_data[i][1] for i in atomic_numbers]

    def __str__(self):
        return self.get_xyz()

    def __hash__(self):
        digest = hashlib.md5(json.dumps((self.get_xyz(), self.alpha_electrons, self.beta_electrons),
                                        sort_keys=True).encode()).hexdigest()
        return int(digest, 16)

    def get_coordinates(self, fragment=None):
        """
        gets the cartesian coordinates

        :param fragment: list of atoms that are part of the fragment

        :return: coordinates list
        """
        if self._coordinates is None:
            self._coordinates = int_to_xyz(self)

        if fragment is None:
            return np.array(self._coordinates).tolist()
        else:
            return np.array(self._coordinates)[fragment].tolist()


    def set_coordinates(self, coordinates):
        """
        sets the cartessian coordinates

        :param coordinates: cartesian coordinates matrix
        """

        self._coordinates = np.array(coordinates)
        self._number_of_atoms = None
        self._energy = {}

    def _get_internal(self):
        if self._internal is None:
            print('No internal coordinates available\n Load internal file')
            exit()
        return self._internal.copy()

    def _set_internal(self, internal):
        self._internal = internal
        self._energy = None
        self._coordinates = int_to_xyz(self)
        self._full_z_matrix = None

    def _get_full_z_matrix(self):
        if self._full_z_matrix is None:
            num_z_atoms = self._get_z_matrix().shape[0]
            self._full_z_matrix = np.zeros((num_z_atoms,6))

            for row, i in enumerate(self._get_z_matrix()[1:]):
                    for col, k in enumerate(i[0]):
                        try:
                            self._full_z_matrix[row+1, col] = float(k)
                        except ValueError:
                            self._full_z_matrix[row+1, col] = self._get_int_dict()[k]


        return self._full_z_matrix

    def _get_z_matrix(self):
        if self._z_matrix is None:
            print('No Z-matrix available\n Load zmatrix file')
            exit()
        return self._z_matrix

    def _set_z_matrix(self, z_matrix):
        self._z_matrix = z_matrix

    def _get_int_label(self):
        return self._int_label

    def _set_int_label(self, int_label):
        self._int_label = int_label

    def _get_int_dict(self):
        self._internal_dict = {}
        for i, label in enumerate(self._get_int_label()[:, 0]):
            self._internal_dict.update({label:self._get_internal()[i, 0]})
        return self._internal_dict

    def _get_int_weights(self):
        return self._int_weights

    def _set_int_weights(self, int_weights):
        self._int_weights = int_weights

    def get_symbols_with_dummy(self):
       # print([i for i in self._atomic_elements if i != "X"])
       return self._symbols


    @property
    def name(self):
        """
        returns the name
        :return: structure name
        """
        return self._name

    @property
    def file_name(self):
        return self._file_name

    @file_name.setter
    def file_name(self, file_name):
        self._file_name = file_name

    @property
    def charge(self):
        """
        returns the charge
        :return: the charge
        """
        return self._charge

    @charge.setter
    def charge(self, charge):
        self._charge = charge

    @property
    def multiplicity(self):
        """
        returns the multiplicity

        :return: the multiplicity
        """
        return self._multiplicity

    @multiplicity.setter
    def multiplicity(self, multiplicity):
        self._multiplicity = multiplicity

    @property
    def number_of_electrons(self):
        """
        returns the total number of electrons

        :return: number of total electrons
        """
        return int(np.sum(self.get_atomic_numbers()) + self.charge)

    @property
    def alpha_electrons(self):
        """
        returns the alpha electrons

        :return: number of alpha electrons
        """
        alpha_unpaired = self.multiplicity // 2 + 1 if (self.number_of_electrons % 2) else self.multiplicity // 2
        return self.number_of_electrons // 2 + alpha_unpaired

    @property
    def beta_electrons(self):
        """
        returns the number of beta electrons

        :return: number of beta electrons
        """
        return self.number_of_electrons - self.alpha_electrons

    def _get_atom_types(self):
        if self._atom_types is None:
            print('No atom types available')
            exit()
        return self._atom_types

    def _set_atom_types(self, atom_types):
        self._atom_types = atom_types

    def get_atomic_numbers(self):
        """
        get the atomic numbers of the atoms of the molecule

        :return: list with the atomic numbers
        """
        if self._atomic_numbers is None:
            self._atomic_numbers = [[data[1].upper() for data in atom_data].index(element.upper())
                                    for element in self.get_symbols()]
        return self._atomic_numbers

    def set_atomic_numbers(self, atomic_numbers):
        self._atomic_numbers = atomic_numbers

    def get_symbols(self):
        """
        get the  atomic element symbols of the atoms of the molecule

        :return: list of symbols
        """
        if self._symbols is None:
            self._symbols = np.array(atom_data)[self.get_atomic_numbers()].T[1]
        return np.array([i for i in self._symbols if i != "X"], dtype=str)

    def set_atomic_elements(self, atomic_elements):
        self._symbols = atomic_elements

    def _get_connectivity(self):
        if self._connectivity is None:
            print('No atom connectivity available')
            exit()

        return self._connectivity

    def _set_connectivity(self, connectivity):
        self._connectivity = connectivity

#   Real methods
    def get_number_of_atoms(self):
        """
        get the number of atoms

        :return: number of atoms
        """
        if self._number_of_atoms is None:
            self._number_of_atoms = np.array(self.get_coordinates()).shape[0]

        return self._number_of_atoms

    def _get_number_of_internal(self):
        if self._number_of_internal is None:
            self._number_of_internal = self._get_internal().shape[0]

        return self._number_of_internal

    def _get_energy(self, method=None):
        if method is None:
            if len(self._energy) == 1:
                return self._energy.values()[0]
            raise Exception('No method defined')
        elif '{}'.format(hash(method)) not in self._energy:
            self._energy['{}'.format(hash(method))] = method.single_point(self)
        return self._energy['{}'.format(hash(method))]

    def _get_modes(self, method=None):
        if self._modes is None:
            if method is None:
                raise Exception('No method defined')
            self._modes, energy = method.vibrations(self)
            self._energy['{}'.format(method.multiplicity)] = energy
        return self._modes

    def get_atomic_masses(self):
        """
        get the atomic masses of the atoms of the molecule

        :return: list of atomic masses
        """
        if self._atomic_masses is None:

            try:
                masses_string = np.array(atom_data)[:, 3:4][[np.where(np.array(atom_data)==element)[0][0]
                                                             for element in self.get_symbols()]]
                self._atomic_masses = np.array(masses_string, dtype=float).T[0]
            except TypeError:
                print('Error reading element labels')
                exit()
        return self._atomic_masses

    def get_valence_electrons(self):
        """
        gets number of valence electrons

        :return: number of valence electrons
        """
        valence_electrons = 0
        for number in self.get_atomic_numbers():
            if 2 >= number > 0:
                valence_electrons += np.mod(number, 2)
            if 18 >= number > 2:
                valence_electrons += np.mod(number-2, 8)
            if 54 >= number > 18:
                valence_electrons += np.mod(number-18, 18)
            if 118 >= number > 54:
                valence_electrons += np.mod(number-54, 32)
            if number > 118:
                raise Exception('Atomic number size not implemented')

        valence_electrons -= self.charge

        return valence_electrons

    def get_xyz(self, title=''):
        """
        generates a XYZ formatted file

        :param title: title of the molecule
        :return: string with the formatted XYZ file
        """
        txt = '{}\n{}\n'.format(self.get_number_of_atoms(), title)
        for s, c in zip(self.get_symbols(), self.get_coordinates()):
            txt += '{:2} '.format(s) + '{:15.10f} {:15.10f} {:15.10f}\n'.format(*c)

        return txt

    def get_connectivity(self, thresh=1.2):
        from scipy.spatial import distance_matrix
        from warnings import warn

        try:
            radius = [atom_data[sym][4] for sym in self.get_atomic_numbers()]
        except KeyError:
            warn('failed to generate connectivity, no connectivity will be used')
            return None

        distances_matrix = distance_matrix(self.get_coordinates(), self.get_coordinates())

        radii_matrix = np.array([radius] * len(radius))
        radii_matrix = radii_matrix + radii_matrix.T

        try:
            relative_differences = np.abs(radii_matrix - distances_matrix) / radii_matrix
        except ValueError:
            warn('failed to generate connectivity')
            return None

        if not (np.array(np.where(relative_differences < thresh - 1)).T + 1).tolist():
            return None
        else:
            return (np.array(np.where(relative_differences < thresh - 1)).T + 1).tolist()


atom_data = [
    # atomic number, symbols, names, masses, bohr radius
    [  0, "X", "X",            0.000000, 0.000],  # 0
    [  1, "H", "Hydrogen",     1.007940, 0.324],  # 1
    [  2, "He", "Helium",      4.002602, 0.000],  # 2
    [  3, "Li", "Lithium",     6.941000, 1.271],  # 3
    [  4, "Be", "Beryllium",   9.012182, 0.927],  # 4
    [  5, "B", "Boron",       10.811000, 0.874],  # 5
    [  6, "C", "Carbon",      12.010700, 0.759],  # 6
    [  7, "N", "Nitrogen",    14.006700, 0.706],  # 7
    [  8, "O", "Oxygen",      15.999400, 0.678],  # 8
    [  9, "F", "Fluorine",    18.998403, 0.568],  # 9
    [ 10, "Ne", "Neon",       20.179700, 0.000],  # 10
    [ 11, "Na", "Sodium",     22.989769, 1.672],  # 11
    [ 12, "Mg", "Magnesium",  24.305000, 1.358],  # 12
    [ 13, "Al", "Aluminium",  26.981539, 1.218],  # 13
    [ 14, "Si", "Silicon",    28.085500, 1.187],  # 14
    [ 15, "P", "Phosphorus",  30.973762, 1.105],  # 15
    [ 16, "S", "Sulfur",      32.065000, 1.045],  # 16
    [ 17, "Cl", "Chlorine",   35.453000, 1.006],  # 17
    [ 18, "Ar", "Argon",      39.948000, 0.000],  # 18
    [ 19, "K", "Potassium",   39.098300, 2.247],  # 19
    [ 20, "Ca", "Calcium",    40.078000, 1.748],  # 20
    [ 21, "Sc", "Scandium",   44.955912, 1.664],  # 21
    [ 22, "Ti", "Titanium",   47.867000, 1.620],  # 22
    [ 23, "V", "Vanadium",    50.941500, 1.543],  # 23
    [ 24, "Cr", "Chromium",   51.996100, 1.418],  # 24
    [ 25, "Mn", "Manganese",  54.938045, 1.569],  # 25
    [ 26, "Fe", "Iron",       55.845000, 1.514],  # 26
    [ 27, "Co", "Cobalt",     58.933195, 1.385],  # 27
    [ 28, "Ni", "Nickel",     58.693400, 1.390],  # 28
    [ 29, "Cu", "Copper",     63.546000, 1.382],  # 29
    [ 30, "Zn", "Zinc",       65.380000, 1.416],  # 30
    [ 31, "Ga", "Gallium",    69.723000, 1.235],  # 31
    [ 32, "Ge", "Germanium",  72.640000, 1.201],  # 32
    [ 33, "As", "Arsenic",    74.921600, 1.232],  # 33
    [ 34, "Se", "Selenium",   78.960000, 1.210],  # 34
    [ 35, "Br", "Bromine",    79.904000, 1.190],  # 35
    [ 36, "Kr", "Krypton",    83.798000, 0.000],  # 36
    [ 37, "Rb", "Rubidium",   85.467800, 2.284],  # 37
    [ 38, "Sr", "Strontium",  87.620000, 1.942],  # 38
    [ 39, "Y", "Yttrium",     88.905850, 1.993],  # 39
    [ 40, "Zr", "Zirconium",  91.224000, 1.758],  # 40
    [ 41, "Nb", "Niobium",    92.906380, 1.610],  # 41
    [ 42, "Mo", "Molybdenum", 95.960000, 1.639],  # 42
    [ 43, "Tc", "Technetium",  0.000000, 1.493],  # 43
    [ 44, "Ru", "Ruthenium",  101.07000, 1.467],  # 44
    [ 45, "Rh", "Rhodium",    102.90550, 1.437],  # 45
    [ 46, "Pd", "Palladium",  106.42000, 1.422],  # 46
    [ 47, "Ag", "Silver",     107.86820, 1.466],  # 47
    [ 48, "Cd", "Cadmium",    112.41100, 1.441],  # 48
    [ 49, "In", "Indium",     114.81800, 1.421],  # 49
    [ 50, "Sn", "Tin",        118.71000, 1.408],  # 50
    [ 51, "Sb", "Antimony",   121.76000, 1.397],  # 51
    [ 52, "Te", "Tellurium",  127.60000, 1.395],  # 52
    [ 53, "I", "Iodine",      126.90447, 1.396],  # 53
    [ 54, "Xe", "Xenon",      131.29300, 1.336],  # 54
    [ 55, "Cs", "Caesium",    132.90545, 2.470],  # 55
    [ 56, "Ba", "Barium",     137.32700, 2.219],  # 56
    [ 57, "La", "Lanthanum",  138.90547, 2.089],  # 57
    [ 58, "Ce", "Cerium",     140.11600, 2.054],  # 58
    [ 59, "Pr", "Praseodymium", 140.90765, 1.979],  # 59
    [ 60, "Nd", "Neodymium",  144.24200, 0.000],  # 60
    [ 61, "Pm", "Promethium",   0.00000, 0.000],  # 61
    [ 62, "Sm", "Samarium",   150.36000, 2.535],  # 62
    [ 63, "Eu", "Europium",   151.96400, 0.000],  # 63
    [ 64, "Gd", "Gadolinium", 157.25000, 0.000],  # 64
    [ 65, "Tb", "Terbium",    158.92535, 0.000],  # 65
    [ 66, "Dy", "Dysprosium", 162.50000, 0.000],  # 66
    [ 67, "Ho", "Holmium",    164.93032, 0.000],  # 67
    [ 68, "Er", "Erbium",     167.25900, 0.000],  # 68
    [ 69, "Tm", "Thulium",    168.93421, 0.000],  # 69
    [ 70, "Yb", "Ytterbium",  173.05400, 0.000],  # 70
    [ 71, "Lu", "Lutetium",   174.96680, 0.000],  # 71
    [ 72, "Hf", "Hafnium",    178.49000, 1.779],  # 72
    [ 73, "Ta", "Tantalum",   180.94788, 1.723],  # 73
    [ 74, "W", "Tungsten",    183.84000, 1.627],  # 74
    [ 75, "Re", "Rhenium",    186.20700, 1.536],  # 75
    [ 76, "Os", "Osmium",     190.23000, 1.521],  # 76
    [ 77, "Ir", "Iridium",    192.21700, 1.456],  # 77
    [ 78, "Pt", "Platinum",   195.08400, 1.390],  # 78
    [ 79, "Au", "Gold",       196.96657, 1.402],  # 79
    [ 80, "Hg", "Mercury",    200.59000, 1.371],  # 80
    [ 81, "Tl", "Thallium",   204.38330, 1.384],  # 81
    [ 82, "Pb", "Lead",       207.20000, 1.820],  # 82
    [ 83, "Bi", "Bismuth",    208.98040, 1.507],  # 83
    [ 84, "Po", "Polonium",     0.00000, 0.000],  # 84
    [ 85, "At", "Astatine",     0.00000, 0.000],  # 85
    [ 86, "Rn", "Radon",        0.00000, 0.000],  # 86
    [ 87, "Fr", "Francium",     0.00000, 0.000],  # 87
    [ 88, "Ra", "Radium",       0.00000, 0.000],  # 88
    [ 89, "Ac", "Actinium",     0.00000, 0.000],  # 89
    [ 90, "Th", "Thorium",    232.03806, 0.000],  # 90
    [ 91, "Pa", "Protactinium",231.03588, 0.000], # 91
    [ 92, "U", "Uranium",     238.02891, 0.000],  # 92
    [ 93, "Np", "Neptunium",    0.00000, 0.000],  # 93
    [ 94, "Pu", "Plutonium",    0.00000, 0.000],  # 94
    [ 95, "Am", "Americium",    0.00000, 0.000],  # 95
    [ 96, "Cm", "Curium",       0.00000, 0.000],  # 96
    [ 97, "Bk", "Berkelium",    0.00000, 0.000],  # 97
    [ 98, "Cf", "Californium",  0.00000, 0.000],  # 98
    [ 99, "Es", "Einsteinium",  0.00000, 0.000],  # 99
    [100, "Fm", "Fermium",      0.00000, 0.000],  # 100
    [101, "Md", "Mendelevium",  0.00000, 0.000],  # 101
    [102, "No", "Nobelium",     0.00000, 0.000],  # 102
    [103, "Lr", "Lawrencium",   0.00000, 0.000],  # 103
    [104, "Rf", "Rutherfordium",0.00000, 0.000],  # 104
    [105, "Db", "Dubnium",      0.00000, 0.000],  # 105
    [106, "Sg", "Seaborgium",   0.00000, 0.000],  # 106
    [107, "Bh", "Bohrium",      0.00000, 0.000],  # 107
    [108, "Hs", "Hassium",      0.00000, 0.000],  # 108
    [109, "Mt", "Meitnerium",   0.00000, 0.000],  # 109
    [110, "Ds", "Darmstadtium", 0.00000, 0.000],  # 110
    [111, "Rg", "Roentgenium",  0.00000, 0.000],  # 111
    [112, "Cn", "Copernicium",  0.00000, 0.000],  # 112
    [113, "Uut", "Ununtrium",   0.00000, 0.000],  # 113
    [114, "Uuq", "Ununquadium", 0.00000, 0.000],  # 114
    [115, "Uup", "Ununpentium", 0.00000, 0.000],  # 115
    [116, "Uuh", "Ununhexium",  0.00000, 0.000],  # 116
    [117, "Uus", "Ununseptium", 0.00000, 0.000],  # 117
    [118, "Uuo", "Ununoctium",  0.00000, 0.000],  # 118
    ]

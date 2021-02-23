from pyqchem.structure import atom_data, Structure

import numpy as np
angstrom_to_bohr = 1/0.529177249


def mat_to_vect(matrix):
    matrix = np.array(matrix)
    tril_matrix = matrix[np.tril_indices_from(matrix)]
    return tril_matrix.tolist()


def get_array_txt(label, type, array, row_size=5):

    formats = {'R': '15.8e',
               'I': '11'}

    array = np.array(array, dtype=float if type == 'R' else int).flatten()
    n_elements = len(array)
    rows = int(np.ceil(n_elements/row_size))

    txt_fchk = '{:40}   {}   N=       {:5}\n'.format(label, type, n_elements)

    # print(rows)
    for i in range(rows):
        if (i+1)*row_size > n_elements:
            txt_fchk += (' {:{fmt}}'* (n_elements - i*row_size) + '\n').format(*array[i * row_size:n_elements],
                                                                               fmt=formats[type])
        else:
            txt_fchk += (' {:{fmt}}'* row_size + '\n').format(*array[i * row_size: (i+1)*row_size],
                                                              fmt=formats[type])

    return txt_fchk


def write_to_fchk(parsed_data, filename):
    txt = build_fchk(parsed_data)
    with open(filename, 'w') as f:
        f.write(txt)


def build_fchk(parsed_data):

    structure = parsed_data['structure']
    basis = parsed_data['basis']
    alpha_mo_coeff = parsed_data['coefficients']['alpha']
    alpha_mo_energies = parsed_data['mo_energies']['alpha']

    #overlap = parsed_data['overlap']
    #coor_shell = parsed_data['coor_shell']
    #core_hamiltonian = parsed_data['core_hamiltonian']

    number_of_basis_functions = len(alpha_mo_coeff)
    number_of_electrons = np.sum(structure.get_atomic_numbers()) - structure.charge
    if structure.multiplicity > 1:
        # raise Exception('1> multiplicity not yet implemented')
        pass

    alpha_mo_coeff = np.array(alpha_mo_coeff).flatten().tolist()

    shell_type_list = {'s':  {'type':  0, 'angular_momentum': 0},
                       'p':  {'type':  1, 'angular_momentum': 1},
                       'd':  {'type':  2, 'angular_momentum': 2},
                       'f':  {'type':  3, 'angular_momentum': 3},
                       'sp': {'type': -1, 'angular_momentum': 1},  # hybrid
                       'd_': {'type': -2, 'angular_momentum': 2},  # pure
                       'f_': {'type': -3, 'angular_momentum': 3}}  # pure

    shell_type = []
    p_exponents = []
    c_coefficients = []
    p_c_coefficients = []
    n_primitives = []
    atom_map = []

    largest_degree_of_contraction = 0
    highest_angular_momentum = 0
    number_of_contracted_shells = 0

    for i, atoms in enumerate(basis['atoms']):
        for shell in atoms['shells']:
            number_of_contracted_shells += 1
            st = shell['shell_type']
            shell_type.append(shell_type_list[st]['type'])
            n_primitives.append(len(shell['p_exponents']))
            atom_map.append(i+1)
            if highest_angular_momentum < shell_type_list[st]['angular_momentum']:
                highest_angular_momentum = shell_type_list[st]['angular_momentum']

            if len(shell['con_coefficients']) > largest_degree_of_contraction:
                    largest_degree_of_contraction = len(shell['con_coefficients'])

            for p in shell['p_exponents']:
                p_exponents.append(p)
            for c in shell['con_coefficients']:
                c_coefficients.append(c)
            for pc in shell['p_con_coefficients']:
                p_c_coefficients.append(pc)

    coordinates_list = angstrom_to_bohr*np.array(structure.get_coordinates()).flatten()

    txt_fchk = '{}\n'.format('filename')
    txt_fchk += 'SP        R                             {}\n'.format(basis['name'] if 'name' in basis else 'no_name')
    txt_fchk += 'Number of atoms                            I               {}\n'.format(structure.get_number_of_atoms())
    txt_fchk += 'Charge                                     I               {}\n'.format(structure.charge)
    txt_fchk += 'Multiplicity                               I               {}\n'.format(structure.multiplicity)
    txt_fchk += 'Number of electrons                        I               {}\n'.format(number_of_electrons)
    txt_fchk += 'Number of alpha electrons                  I               {}\n'.format(structure.alpha_electrons)
    txt_fchk += 'Number of beta electrons                   I               {}\n'.format(structure.beta_electrons)

    txt_fchk += get_array_txt('Atomic numbers', 'I', structure.get_atomic_numbers(), row_size=6)
    txt_fchk += get_array_txt('Current cartesian coordinates', 'R', coordinates_list)
    txt_fchk += get_array_txt('Nuclear charges', 'R', structure.get_atomic_numbers())

    txt_fchk += 'Number of basis functions                  I               {}\n'.format(number_of_basis_functions)
    txt_fchk += 'Number of contracted shells                I               {}\n'.format(number_of_contracted_shells)
    txt_fchk += 'Number of primitive shells                 I               {}\n'.format(np.sum(n_primitives))
    txt_fchk += 'Highest angular momentum                   I               {}\n'.format(highest_angular_momentum)
    txt_fchk += 'Largest degree of contraction              I               {}\n'.format(largest_degree_of_contraction)

    txt_fchk += get_array_txt('Shell types', 'I', shell_type, row_size=6)
    txt_fchk += get_array_txt('Number of primitives per shell', 'I', n_primitives, row_size=6)
    txt_fchk += get_array_txt('Shell to atom map', 'I', atom_map, row_size=6)
    txt_fchk += get_array_txt('Primitive exponents', 'R', p_exponents)
    txt_fchk += get_array_txt('Contraction coefficients', 'R', c_coefficients)
    txt_fchk += get_array_txt('P(S=P) Contraction coefficients', 'R', p_c_coefficients)
    # txt_fchk += get_array_txt('Coordinates of each shell', 'R', coor_shell) #
    # txt_fchk += get_array_txt('Overlap Matrix', 'R', overlap)
    #txt_fchk += get_array_txt('Core Hamiltonian Matrix', 'R', core_hamiltonian)
    txt_fchk += get_array_txt('Alpha Orbital Energies', 'R', alpha_mo_energies)
    # txt_fchk += get_array_txt('Beta Orbital Energies', 'R', beta_mo_energies)

    if 'scf_density' in parsed_data:
        scf_density = mat_to_vect(parsed_data['scf_density'])
        txt_fchk += get_array_txt('Total SCF Density', 'R', scf_density)

    txt_fchk += get_array_txt('Alpha MO coefficients', 'R', alpha_mo_coeff)

    if 'beta' in parsed_data['coefficients']:
        beta_mo_coeff = parsed_data['coefficients']['beta']
        beta_mo_coeff = np.array(beta_mo_coeff).flatten().tolist()
        beta_mo_energies = parsed_data['mo_energies']['beta']

        txt_fchk += get_array_txt('Beta MO coefficients', 'R', beta_mo_coeff)
        txt_fchk += get_array_txt('Beta Orbital Energies', 'R', beta_mo_energies)

    return txt_fchk


def read_structure_from_xyz(filename, read_multiple=False):
    """
    Read a XYZ file a return a Structure

    :param file_name: file name
    :param read_multiple: read multiple molecules
    :return: Structure or list of Structure
    """

    with open(filename, mode='r') as lines:
        file_lines = lines.readlines()
        i_start = 0
        structures = []
        while True:
            n_atoms = int(file_lines[i_start].split()[0])
            name = file_lines[i_start+ 1]
            body_txt = file_lines[i_start + 2: i_start+ n_atoms+2]

            structures.append(Structure(coordinates=np.array([c.split()[1:4] for c in body_txt], dtype=float).tolist(),
                                        symbols=[c.split()[0] for c in body_txt],
                                        name=name)
                              )

            if not read_multiple:
                return structures[0]

            i_start += len(body_txt) + 2
            if len(file_lines) - i_start < 3:
                break

    return structures


def write_structure_to_xyz(structures, filename):
    """
    Store structures in a XYZ file

    :param structures: Structure or list of Structure
    :param filename: name of the file
    :return:
    """

    # Get xyz file format from Molecule or Geometry (or list)

    if not isinstance(structures, list):
        structures = [structures]

    txt = ''
    for structure in structures:
        txt += '{}\n'.format(structure.get_number_of_atoms())
        txt += '{}\n'.format(structure.name if structure.name is not None else '')
        for idp, position in enumerate(structure.get_coordinates()):
            txt += '{:2} {:11.6f} {:11.6f} {:11.6f}\n'.format(structure.get_symbols()[idp],
                                                                  position[0], position[1], position[2])
    with open(filename, 'w') as f:
        f.write(txt)


if __name__ == '__main__':
    from pyqchem.parsers.parser_fchk import parser_fchk
    txt_fchk = open('qchem_temp_32947.fchk', 'r').read()
    parsed_data = parser_fchk(txt_fchk)
    txt_fchk_new = build_fchk(parsed_data)
    with open('test.fchk', 'w') as f:
        f.write(txt_fchk_new)

    #structure, basis, alpha_coeff, beta_coeff = get_data_from_file_fchk('qchem_temp_32947.fchk')
    #print(basis)
from pyqchem.structure import atom_data, Structure

import numpy as np


def reformat_input(array):
    flat_list = []
    for sublist in array:
        for item in sublist:
            if len(item) > 2:
                flat_list.append(item)
            else:
                flat_list.append(item)
    return flat_list


def basis_format(basis_set_name,
                 atomic_numbers,
                 atomic_symbols,
                 shell_type,
                 n_primitives,
                 atom_map,
                 p_exponents,
                 c_coefficients,
                 p_c_coefficients):

    # print(n_primitives)

    typeList = {'0': ['s', 1],
                '1': ['p', 3],
                '2': ['d', 6],
                '3': ['f', 10],
                '-1': ['sp', 4],
                '-2': ['d', 5],
                '-3': ['f', 7]}

    atomic_numbers = np.array(atomic_numbers, dtype=int)
    atom_map = np.array(atom_map, dtype=int)
    # print(atom_map)
    basis_set = {'name': basis_set_name,
                 'primitive_type': 'gaussian'}

    shell_type_index = [0] + np.cumsum([typeList['{}'.format(s)][1]
                                        for s in shell_type]).tolist()
    prim_from_shell_index = [0] + np.cumsum(np.array(n_primitives, dtype=int)).tolist()

    # print(shell_type_index)
    # print(prim_from_shell_index)

    atoms_data = []
    for iatom, atomic_number in enumerate(atomic_numbers):
        symbol = atomic_symbols[iatom]

        shell_from_atom_counts = np.unique(atom_map, return_counts=True)[1]
        shell_from_atom_index = np.unique(atom_map, return_index=True)[1]
        #print(shell_from_atom_counts)
        #print('atom_indexes', shell_from_atom_index)
        #print('atom_number', iatom)
        #print('shells index', shell_from_atom_index[iatom])
        #print('number of shells', shell_from_atom_counts[iatom])

        shells_data = []
        for ishell in range(shell_from_atom_counts[iatom]):
            st = typeList['{}'.format(shell_type[shell_from_atom_index[iatom] + ishell])]
            #print(st, ishell)
            ini_prim = prim_from_shell_index[shell_from_atom_index[iatom] + ishell]
            fin_prim = prim_from_shell_index[shell_from_atom_index[iatom] + ishell+1]
            #print(ini_prim)
            #print(fin_prim)

            shells_data.append({
                'shell_type': st[0],
                'p_exponents': p_exponents[ini_prim: fin_prim],
                'con_coefficients': c_coefficients[ini_prim: fin_prim],
                'p_con_coefficients': p_c_coefficients[ini_prim: fin_prim],
            })

        atoms_data.append({'shells': shells_data,
                           'symbol': symbol,
                           'atomic_number': atomic_number})

    basis_set['atoms'] = atoms_data

    return basis_set


def get_data_from_file_fchk(file_name):
    key_list = ['Charge', 'Multiplicity', 'Atomic numbers', 'Current cartesian coordinates',
                'Shell type', 'Number of primitives per shell', 'Shell to atom map', 'Primitive exponents',
                'Contraction coefficients', 'P(S=P) Contraction coefficients',
                'Alpha MO coefficients', 'Beta MO coefficients']
    input_molecule = [[] for _ in range(len(key_list))]
    read = False
    with open(file_name, mode='r') as lines:
        name = lines.readline()
        line = lines.readline().split()
        basis_set = line[-1]
        if 'R' in line[1]:
            del key_list[-1]

        n = 1
        options = True
        for line in lines:
            if read:
                try:
                    float(line.split()[0])
                    input_molecule[n].append(line.split())
                except ValueError:
                    input_molecule[n] = reformat_input(input_molecule[n])
                    read = False

            for idn, key in enumerate(key_list):
                if key in line:
                    if n == len(key_list) - 1:
                        break
                    if options and idn != 2:
                        input_molecule[idn].append(int(line.split()[-1]))
                        n = idn
                        break
                    else:
                        options = False
                    if n == idn:
                        n += 1
                    else:
                        n = idn
                    read = True
                    break

        bohr_to_angstrom = 0.529177249
        coordinates = np.array(input_molecule[3], dtype=float).reshape(-1, 3) * bohr_to_angstrom
        structure = Structure(coordinates=coordinates,
                              atomic_numbers=[int(num) for num in input_molecule[2]])

        print(input_molecule[7])

        basis = basis_format(basis_set_name=basis_set,
                             atomic_numbers=structure.get_atomic_numbers(),
                             atomic_symbols=structure.get_atomic_elements(),
                             shell_type=[int(num) for num in input_molecule[4]],
                             n_primitives=[int(num) for num in input_molecule[5]],
                             atom_map=[int(num) for num in input_molecule[6]],
                             p_exponents=[float(num) for num in input_molecule[7]],
                             c_coefficients=[float(num) for num in input_molecule[8]],
                             p_c_coefficients=[float(num) for num in input_molecule[9]])

        alpha_coeff = [float(num) for num in input_molecule[10]]
        beta_coeff = [float(num) for num in input_molecule[11]]
        if len(beta_coeff) == 0:
            beta_coeff = alpha_coeff

        return structure, basis, alpha_coeff, beta_coeff


def get_array_txt(label, type, array, row_size=5):

    formats = {'R': '15.8e',
               'I': '11'}

    n_elements = len(array)
    rows = int(np.ceil(n_elements/row_size))

    txt_fchk = '{:40}   {}   N=       {:5}\n'.format(label, type, n_elements)

    # print(rows)
    for i in range(rows):
        if (i+1)*row_size > n_elements:
            txt_fchk += (' {:{fmt}}'* (n_elements - i*row_size) + '\n').format(*array[i * row_size:n_elements],
                                                                               fmt=formats[type])
        else:
            txt_fchk += (' {:{fmt}}'* row_size  + '\n').format(*array[i * row_size: (i+1)*row_size],
                                                               fmt=formats[type])

    return txt_fchk


def build_fchk(parsed_data):

    structure = parsed_data['structure']
    basis = parsed_data['basis']
    alpha_mo_coeff = parsed_data['coefficients']['alpha']
    alpha_mo_energies = parsed_data['mo_energies']['alpha']

    #overlap = parsed_data['overlap']
    #coor_shell = parsed_data['coor_shell']
    #core_hamiltonian = parsed_data['core_hamiltonian']
    #scf_density = parsed_data['scf_density']

    number_of_basis_functions = len(alpha_mo_coeff)
    number_of_electrons = np.sum(structure.get_atomic_numbers()) - structure.charge
    if structure.multiplicity > 1:
        raise Exception('1> multiplicity not yet implemented')
    alpha_electrons = number_of_electrons // 2
    beta_electrons = number_of_electrons // 2
    #print(alpha_electrons)
    #print(number_of_electrons)

    alpha_mo_coeff = np.array(alpha_mo_coeff).flatten().tolist()

    if 'beta' in parsed_data['coefficients']:
        beta_mo_coeff = parsed_data['coefficients']['beta']
        beta_mo_coeff = np.array(beta_mo_coeff).flatten().tolist()

        beta_mo_energies = parsed_data['mo_energies']['beta']

    shell_type_list = {'s':  {'type':  0, 'angular_momentum': 0},
                       'p':  {'type':  1, 'angular_momentum': 1},
                       'd':  {'type':  2, 'angular_momentum': 2},
                       'f':  {'type':  3, 'angular_momentum': 3},
                       'sp': {'type': -1, 'angular_momentum': 1},  # hybrid
                       'dc': {'type': -2, 'angular_momentum': 2},  # pure
                       'fc': {'type': -3, 'angular_momentum': 3}}  # pure

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

    angstrom_to_bohr = 1/0.529177249
    coordinates_list = angstrom_to_bohr*structure.get_coordinates().flatten()

    txt_fchk = '{}\n'.format('filename')
    txt_fchk += 'SP        R                             {}\n'.format('6-31G')
    txt_fchk += 'Number of atoms                            I               {}\n'.format(structure.get_number_of_atoms())
    txt_fchk += 'Charge                                     I               {}\n'.format(structure.charge)
    txt_fchk += 'Multiplicity                               I               {}\n'.format(structure.multiplicity)
    txt_fchk += 'Number of electrons                        I               {}\n'.format(number_of_electrons)
    txt_fchk += 'Number of alpha electrons                  I               {}\n'.format(alpha_electrons)
    txt_fchk += 'Number of beta electrons                   I               {}\n'.format(beta_electrons)

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
    # txt_fchk += get_array_txt('Total SCF Density', 'R', scf_density)
    txt_fchk += get_array_txt('Alpha MO coefficients', 'R', alpha_mo_coeff)
    # txt_fchk += get_array_txt('Beta MO coefficients', 'R', beta_mo_coeff)

    return txt_fchk

if __name__ == '__main__':
    from pyqchem.parsers.parser_fchk import parser_fchk
    txt_fchk = open('qchem_temp_32947.fchk', 'r').read()
    parsed_data = parser_fchk(txt_fchk)
    txt_fchk_new = build_fchk(parsed_data)
    with open('test.fchk', 'w') as f:
        f.write(txt_fchk_new)

    #structure, basis, alpha_coeff, beta_coeff = get_data_from_file_fchk('qchem_temp_32947.fchk')
    #print(basis)
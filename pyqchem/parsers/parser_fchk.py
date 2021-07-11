import numpy as np
from pyqchem.structure import Structure


def basis_format(basis_set_name,
                 atomic_numbers,
                 atomic_symbols,
                 shell_type,
                 n_primitives,
                 atom_map,
                 p_exponents,
                 c_coefficients,
                 p_c_coefficients):
    """
    Function to generate standard basis dictionary

    :param basis_set_name:
    :param atomic_numbers:
    :param atomic_symbols:
    :param shell_type:
    :param n_primitives:
    :param atom_map:
    :param p_exponents:
    :param c_coefficients:
    :param p_c_coefficients:
    :return:
    """

    # print(n_primitives)

    typeList = {'0': ['s', 1],
                '1': ['p', 3],
                '2': ['d', 6],
                '3': ['f', 10],
                '-1': ['sp', 4],
                '-2': ['d_', 5],
                '-3': ['f_', 7]}

    atomic_numbers = [int(an) for an in atomic_numbers]
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
        symbol = str(atomic_symbols[iatom])

        shell_from_atom_counts = np.unique(atom_map, return_counts=True)[1]
        shell_from_atom_index = np.unique(atom_map, return_index=True)[1]
        # print(shell_from_atom_counts)
        # print('atom_indexes', shell_from_atom_index)
        # print('atom_number', iatom)
        # print('shells index', shell_from_atom_index[iatom])
        # print('number of shells', shell_from_atom_counts[iatom])

        shells_data = []
        for ishell in range(shell_from_atom_counts[iatom]):
            st = typeList['{}'.format(shell_type[shell_from_atom_index[iatom] + ishell])]
            # print(st, ishell)
            ini_prim = prim_from_shell_index[shell_from_atom_index[iatom] + ishell]
            fin_prim = prim_from_shell_index[shell_from_atom_index[iatom] + ishell+1]
            # print(ini_prim)
            # print(fin_prim)

            shells_data.append({
                'shell_type': st[0],
                'functions': st[1],
                'p_exponents': p_exponents[ini_prim: fin_prim],
                'con_coefficients': c_coefficients[ini_prim: fin_prim],
                'p_con_coefficients': p_c_coefficients[ini_prim: fin_prim],
            })

        atoms_data.append({'shells': shells_data,
                           'symbol': symbol,
                           'atomic_number': atomic_number})

    basis_set['atoms'] = atoms_data

    return basis_set


def _reformat_input(array):
    flat_list = []
    for sublist in array:
        for item in sublist:
            if len(item) > 2:
                flat_list.append(item)
            else:
                flat_list.append(item)
    return flat_list


def vect_to_mat(vector):
    n = int(np.sqrt(0.25 + 2 * len(vector)) - 0.5)

    k = 0
    matrix = np.zeros([n, n])
    for i in range(n):
        for j in range(0, i+1):
            matrix[i, j] = vector[k + j]
            matrix[j, i] = vector[k + j]
        k += i+1

    return matrix


def _get_all_nato(output):
    import re

    nato_coefficients_list = []
    nato_occupancies_list = []

    for m in re.finditer('Alpha NATO coefficients', output):
        n_elements = int(output[m.end():m.end() + 100].replace('\n', ' ').split()[2])
        nbas = int(np.sqrt(n_elements))
        data = output[m.end(): m.end() + n_elements*50].split()[3:n_elements+3]
        nato_coefficients_list.append({'alpha': np.array(data, dtype=float).reshape(nbas, nbas).tolist()})

    for m in re.finditer('Alpha Natural Orbital occupancies', output):
        n_elements = int(output[m.end():m.end() + 100].replace('\n', ' ').split()[2])
        data = output[m.end(): m.end() + n_elements*50].split()[3:n_elements+3]
        nato_occupancies_list.append({'alpha': np.array(data, dtype=float).tolist()})

    for i, m in enumerate(re.finditer('Beta NATO coefficients', output)):
        n_elements = int(output[m.end():m.end() + 100].replace('\n', ' ').split()[2])
        nbas = int(np.sqrt(n_elements))
        data = output[m.end(): m.end() + n_elements*50].split()[3:n_elements+3]
        nato_coefficients_list[i]['beta'] = np.array(data, dtype=float).reshape(nbas, nbas).tolist()

    for i, m in enumerate(re.finditer('Beta Natural Orbital occupancies', output)):
        n_elements = int(output[m.end():m.end() + 100].replace('\n', ' ').split()[2])
        data = output[m.end(): m.end() + n_elements*50].split()[3:n_elements+3]
        nato_occupancies_list[i]['beta'] = np.array(data, dtype=float).tolist()

    return nato_coefficients_list, nato_occupancies_list


def _get_all_nto(output):
    import re

    nto_coefficients_list = []
    nto_occupancies_list = []

    for m in re.finditer('Natural Transition Orbital U coefficients', output):
        n_elements = int(output[m.end():m.end() + 100].replace('\n', ' ').split()[2])
        nbas = int(np.sqrt(n_elements))
        data = output[m.end(): m.end() + n_elements*50].split()[3:n_elements+3]
        nto_coefficients_list.append({'U': np.array(data, dtype=float).reshape(nbas, nbas).tolist()})

    for m in re.finditer('Natural Transition Orbital occupancies', output):
        n_elements = int(output[m.end():m.end() + 100].replace('\n', ' ').split()[2])
        data = output[m.end(): m.end() + n_elements*50].split()[3:n_elements+3]
        nto_occupancies_list.append(np.array(data, dtype=float).tolist())

    for i, m in enumerate(re.finditer('Natural Transition Orbital V coefficients', output)):
        n_elements = int(output[m.end():m.end() + 100].replace('\n', ' ').split()[2])
        nbas = int(np.sqrt(n_elements))
        data = output[m.end(): m.end() + n_elements*50].split()[3:n_elements+3]
        nto_coefficients_list[i]['V'] = np.array(data, dtype=float).reshape(nbas, nbas).tolist()

    return nto_coefficients_list, nto_occupancies_list


def parser_fchk(output):

    def convert_to_type(item_type, item):
        item_types = {'I': int,
                      'R': float}

        if type(item) is list:
            return [item_types[item_type](e) for e in item]
        else:
            return item_types[item_type](item)

    key_list = ['Charge', 'Multiplicity', 'Number of alpha electrons', 'Number of beta electrons',
                'Atomic numbers', 'Current cartesian coordinates', 'Number of basis functions', 'Shell types',
                'Number of primitives per shell', 'Shell to atom map', 'Primitive exponents',
                'Contraction coefficients', 'P(S=P) Contraction coefficients', 'Alpha MO coefficients',
                'Beta MO coefficients', 'Coordinates of each shell', 'Overlap Matrix',
                'Core Hamiltonian Matrix', 'Alpha Orbital Energies', 'Beta Orbital Energies',
                'Total SCF Density', 'Alpha NATO coefficients', 'Beta NATO coefficients',
                'Alpha Natural Orbital occupancies', 'Beta Natural Orbital occupancies',
                'Natural Transition Orbital occupancies', 'Natural Transition Orbital U coefficients',
                'Natural Transition Orbital V coefficients'
                ]

    basis_set = output.split('\n')[1].split()[-1]
    words_output = output.replace('\n', ' ').split()

    data = {}
    nw = len(words_output)
    for key in key_list:
        wc = len(key.split())
        for i in range(nw):
            word = ' '.join(words_output[i:i+wc])
            if word == key:
                item_type = words_output[i+wc]
                if words_output[i + wc + 1] == 'N=':
                    n_elements = int(words_output[i + wc + 2])
                    data[word] = convert_to_type(item_type, words_output[i + wc + 3: i + wc + n_elements + 3])
                else:
                    data[word] = convert_to_type(item_type, words_output[i + wc + 1])
                break

    bohr_to_angstrom = 0.529177249

    coordinates = np.array(data['Current cartesian coordinates']).reshape(-1, 3) * bohr_to_angstrom
    structure = Structure(coordinates=coordinates,
                          atomic_numbers=data['Atomic numbers'],
                          multiplicity=data['Multiplicity'],
                          charge=data['Charge'])

    if not 'P(S=P) Contraction coefficients' in data:
        data['P(S=P) Contraction coefficients'] = np.zeros_like(data['Contraction coefficients']).tolist()

    basis = basis_format(basis_set_name=basis_set,
                         atomic_numbers=structure.get_atomic_numbers(),
                         atomic_symbols=structure.get_symbols(),
                         shell_type=data['Shell types'],
                         n_primitives=data['Number of primitives per shell'],
                         atom_map=data['Shell to atom map'],
                         p_exponents=data['Primitive exponents'],
                         c_coefficients=data['Contraction coefficients'],
                         p_c_coefficients=data['P(S=P) Contraction coefficients'])

    nbas = data['Number of basis functions']

    final_dict = {'structure': structure,
                  'basis': basis,
                  'number_of_electrons': {'alpha': data['Number of alpha electrons'],
                                          'beta': data['Number of beta electrons']}
                  }

    if 'Alpha MO coefficients' in data:
        final_dict['coefficients'] = {'alpha': np.array(data['Alpha MO coefficients']).reshape(nbas, -1).tolist()}
        final_dict['mo_energies'] = {'alpha': data['Alpha Orbital Energies']}

    if 'Beta MO coefficients' in data:
        final_dict['coefficients']['beta'] = np.array(data['Beta MO coefficients']).reshape(nbas, -1).tolist()
        final_dict['mo_energies']['beta'] = data['Beta Orbital Energies']

    if 'Total SCF Density' in data:
        final_dict['scf_density'] = vect_to_mat(data['Total SCF Density'])

    if 'Core Hamiltonian Matrix' in data:
        final_dict['scf_density'] = vect_to_mat(data['Core Hamiltonian Matrix']).tolist()

    if 'Overlap Matrix' in data:
        final_dict['overlap'] = vect_to_mat(data['Overlap Matrix']).tolist()

    if 'Alpha NATO coefficients' in data:
        final_dict['nato_coefficients'] = {'alpha': np.array(data['Alpha NATO coefficients']).reshape(nbas, -1).tolist()}
        final_dict['nato_occupancies'] = {'alpha': data['Alpha Natural Orbital occupancies']}

    if 'Beta NATO coefficients' in data:
        final_dict['nato_coefficients'].update({
            'beta': np.array(data['Beta NATO coefficients']).reshape(nbas, -1).tolist()})
        final_dict['nato_occupancies'].update({'beta': data['Beta Natural Orbital occupancies']})

    # check multiple NATO (may be improved)
    if 'Alpha NATO coefficients' in data:
        nato_coefficients_list, nato_occupancies_list = _get_all_nato(output)
        if len(nato_occupancies_list) > 1:
            final_dict['nato_coefficients_multi'] = nato_coefficients_list
            final_dict['nato_occupancies_multi'] = nato_occupancies_list

    if 'Natural Transition Orbital occupancies' in data:
        nat_coefficients_list, nat_occupancies_list = _get_all_nto(output)
        if len(nat_occupancies_list) > 1:
            final_dict['nto_coefficients_multi'] = nat_coefficients_list
            final_dict['nto_occupancies_multi'] = nat_occupancies_list

    return final_dict

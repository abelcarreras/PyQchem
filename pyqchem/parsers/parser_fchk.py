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

    print(n_primitives)

    typeList = {'0': ['s', 1],
                '1': ['p', 3],
                '2': ['d', 6],
                '3': ['f', 10],
                '-1': ['sp', 4],
                '-2': ['d', 5],
                '-3': ['f', 7]}

    atomic_numbers = np.array(atomic_numbers, dtype=int)
    atom_map = np.array(atom_map, dtype=int)
    print(atom_map)
    basis_set = {'name': basis_set_name,
                 'primitive_type': 'gaussian'}

    shell_type_index = [0] + np.cumsum([typeList['{}'.format(s)][1]
                                        for s in shell_type]).tolist()
    prim_from_shell_index = [0] + np.cumsum(np.array(n_primitives, dtype=int)).tolist()

    print(shell_type_index)
    print(prim_from_shell_index)

    atoms_data = []
    for iatom, atomic_number in enumerate(atomic_numbers):
        symbol = atomic_symbols[iatom]

        shell_from_atom_counts = np.unique(atom_map, return_counts=True)[1]
        shell_from_atom_index = np.unique(atom_map, return_index=True)[1]
        #print(shell_from_atom_counts)
        print('atom_indexes', shell_from_atom_index)
        print('atom_number', iatom)
        print('shells index', shell_from_atom_index[iatom])
        print('number of shells', shell_from_atom_counts[iatom])

        shells_data = []
        for ishell in range(shell_from_atom_counts[iatom]):
            st = typeList['{}'.format(shell_type[shell_from_atom_index[iatom] + ishell])]
            print(st, ishell)
            ini_prim = prim_from_shell_index[shell_from_atom_index[iatom] + ishell]
            fin_prim = prim_from_shell_index[shell_from_atom_index[iatom] + ishell+1]
            print(ini_prim)
            print(fin_prim)

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


def parser_fchk(output):
    key_list = ['Charge', 'Multiplicity', 'Atomic numbers', 'Current cartesian coordinates',
                'Shell type', 'Number of primitives per shell', 'Shell to atom map', 'Primitive exponents',
                'Contraction coefficients', 'P(S=P) Contraction coefficients',
                'Alpha MO coefficients', 'Beta MO coefficients',
                'Coordinates of each shell', 'Overlap Matrix', 'Core Hamiltonian Matrix', 'Alpha Orbital Energies', 'Beta Orbital Energies',
                'Total SCF Density'
                ]
    input_molecule = [[] for _ in range(len(key_list))]
    read = False
    #with open(file_name, mode='r') as lines:
    for line in output.split('\n'):
        #name = lines.readline()
        #line = lines.readline().split()
        basis_set = line[-1]
        if 'R' in line[1]:
            del key_list[-1]

        n = 1
        options = True
        for line in output.split('\n'):
            if read:
                try:
                    float(line.split()[0])
                    input_molecule[n].append(line.split())
                except ValueError:
                    input_molecule[n] = _reformat_input(input_molecule[n])
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

        print('multiplicity:', int(input_molecule[1][0]))
        print('charge:', int(input_molecule[0][0]))

        bohr_to_angstrom = 0.529177249
        coordinates = np.array(input_molecule[3], dtype=float).reshape(-1, 3) * bohr_to_angstrom
        structure = Structure(coordinates=coordinates,
                              atomic_numbers=[int(num) for num in input_molecule[2]],
                              multiplicity=int(input_molecule[1][0]),
                              charge=int(input_molecule[0][0]))

        basis = basis_format(basis_set_name=basis_set,
                             atomic_numbers=structure.get_atomic_numbers(),
                             atomic_symbols=structure.get_atomic_elements(),
                             shell_type=[int(num) for num in input_molecule[4]],
                             n_primitives=[int(num) for num in input_molecule[5]],
                             atom_map=[int(num) for num in input_molecule[6]],
                             p_exponents=[float(num) for num in input_molecule[7]],
                             c_coefficients=[float(num) for num in input_molecule[8]],
                             p_c_coefficients=[float(num) for num in input_molecule[9]])

        #alpha_mo_coeff = [float(num) for num in input_molecule[10]]
        #beta_mo_coeff = [float(num) for num in input_molecule[11]]

        nbas = int(np.sqrt(len(input_molecule[10])))
        alpha_mo_coeff = np.array(input_molecule[10], dtype=float).reshape(nbas, nbas).tolist()
        alpha_mo_energy = [float(num) for num in input_molecule[15]]

        if len(input_molecule[11]) == 0:
            beta_mo_coeff = alpha_mo_coeff
            beta_mo_energy = alpha_mo_energy
        else:
            beta_mo_coeff = np.array(input_molecule[11], dtype=float).reshape(nbas, nbas)
            beta_mo_energy = [float(num) for num in input_molecule[16]]

        # coor_shell = [float(num) for num in input_molecule[12]]
        # overlap = [float(num) for num in input_molecule[13]]
        # core_hamiltonian = [float(num) for num in input_molecule[14]]
        # scf_density = [float(num) for num in input_molecule[17]]

        return {'structure': structure,
                'basis': basis,
                'coefficients': {'alpha': alpha_mo_coeff,
                                 'beta': beta_mo_coeff},
                'mo_energies': {'alpha': alpha_mo_energy,
                                'beta': beta_mo_energy},
                #'overlap': overlap,
                #'coor_shell': coor_shell,
                #'core_hamiltonian': core_hamiltonian,
                #'scf_density': scf_density
                }

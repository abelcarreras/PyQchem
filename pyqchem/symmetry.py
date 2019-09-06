from wfnsympy import WfnSympy
import numpy as np


def get_wf_symmetry(structure,
                    basis,
                    mo_coeff,
                    center=(0., 0., 0.),
                    vaxis=(0., 0., 1.)):

    # print('valence', structure.get_valence_electrons())
    # print(structure.get_xyz())

    type_list = {'s': 0,
                 'p': 1,
                 'd': 2,
                 'f': 3,
                 'sp': -1,
                 'dc': -2,  # pure
                 'fc': -3}  # pure

    shell_type = []
    p_exponents = []
    c_coefficients = []
    p_c_coefficients = []
    n_primitives = []
    atom_map = []
    for i, atoms in enumerate(basis['atoms']):
        for shell in atoms['shells']:
            st = shell['shell_type']
            shell_type.append(type_list[st])
            n_primitives.append(len(shell['p_exponents']))
            atom_map.append(i+1)
            for p in shell['p_exponents']:
                p_exponents.append(p)
            for c in shell['con_coefficients']:
                c_coefficients.append(c)
            for pc in shell['p_con_coefficients']:
                p_c_coefficients.append(pc)

    alpha_mo_coeff = np.array(mo_coeff['alpha']).flatten().tolist()
    if 'beta' in mo_coeff:
        # print(mo_coeff)
        exit()
        beta_mo_coeff = np.array(mo_coeff['beta']).flatten().tolist()
    else:
        beta_mo_coeff = None

    # for s, c in zip(structure.get_atomic_elements(), structure.get_coordinates()):
    #     print('{:2} '.format(s) + '{:10.5f} {:10.5f} {:10.5f}'.format(*c))

    molsym = WfnSympy(coordinates=structure.get_coordinates(),
                      symbols=structure.get_atomic_elements(),
                      basis=basis,
                      center=center,
                      VAxis=vaxis,
                      alpha_mo_coeff=alpha_mo_coeff,
                      beta_mo_coeff=beta_mo_coeff,
                      charge=structure.charge,
                      multiplicity=structure.multiplicity,
                      group='C2h')

    return molsym


def get_orbital_type(molsym):

    sh_index = molsym.SymLab.index('s_h')
    orbital_type = []
    for i, overlap in enumerate(molsym.mo_SOEVs_a[:, sh_index]):
        if overlap < 0:
            orbital_type.append(['pi', np.abs(overlap)])
        else:
            orbital_type.append(['sigma', np.abs(overlap)])
    return orbital_type


def set_zero_coefficients(basis, mo_coeff, range_atoms):
    """
    set 0.0 to coefficients of functions centered in atoms 'range_atoms'

    :param basis: full basis dictionary
    :param mo_coeff: full molecular orbitals coefficients to be modiffied
    :param range_atoms: list containing the atom numbers whose coefficients will be set to zero
    :return:
    """
    functions_to_atom = []

    nat = len(basis['atoms'])
    for i in range(0, nat):
        nf = 0
        for shell in basis['atoms'][i]['shells']:
            nf += shell['functions']
        functions_to_atom.append(nf)
    functions_to_atom = functions_to_atom

    # print(funtions_to_atom)
    mo_coeff = np.array(mo_coeff)
    for i in range_atoms:
        ini = np.sum(functions_to_atom[:i], dtype=int)
        fin = np.sum(functions_to_atom[:i+1], dtype=int)
        # print('ini', ini, 'fin', fin)
        mo_coeff[:, ini: fin] *= 0.0

    return mo_coeff.tolist()


if __name__ == '__main__':
    from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input
    from pyqchem.structure import Structure
    from pyqchem.parsers.parser_fchk import parser_fchk
    from pyqchem.file_io import build_fchk

    benzene_coordinates = [[ 0.00000,  1.40272, 0.00000],
                           [ 0.00000,  2.49029, 0.00000],
                           [-1.21479,  0.70136, 0.00000],
                           [-2.15666,  1.24515, 0.00000],
                           [-1.21479, -0.70136, 0.00000],
                           [-2.15666, -1.24515, 0.00000],
                           [ 0.00000, -1.40272, 0.00000],
                           [ 0.00000, -2.49029, 0.00000],
                           [ 1.21479, -0.70136, 0.00000],
                           [ 2.15666, -1.24515, 0.00000],
                           [ 1.21479,  0.70136, 0.00000],
                           [ 2.15666,  1.24515, 0.00000]]

    symbols = ['C','H','C','H','C','H','C','H','C','H','C','F']

    benzene_coordinates = np.array(benzene_coordinates)
    molecule = Structure(coordinates=benzene_coordinates,
                         atomic_elements=symbols,
                         charge=0)

    parameters = {'jobtype': 'sp',
                  'exchange': 'hf',
                  'basis': '6-31G',
                  'gui': 2}

    txt_input = create_qchem_input(molecule, **parameters)

    txt_fchk = get_output_from_qchem(txt_input, processors=4, force_recalculation=True,
                                     parser=None, read_fchk=True)

    #txt_fchk = open('qchem_temp_32947.fchk', 'r').read()
    parsed_data = parser_fchk(txt_fchk)
    open('test.fchk', 'w').write(txt_fchk)

    alpha_mo_coeff = set_zero_coefficients(parsed_data['basis'],
                                           parsed_data['coefficients']['alpha'],
                                           range(6, 12))
    parsed_data['coefficients']['alpha'] = alpha_mo_coeff

    structure = parsed_data['structure']
    print(structure.get_xyz())

    txt_fchk = build_fchk(parsed_data)
    open('test2.fchk', 'w').write(txt_fchk)

    z_dist = structure.get_coordinates()[0][2]

    molsym = get_wf_symmetry(parsed_data['structure'],
                             parsed_data['basis'],
                             parsed_data['coefficients'],
                             center=[0.0, 0.0, z_dist])

    orbital_type = get_orbital_type(molsym)
    print('  {:5} {:5} {:5}'.format('num', 'type', 'overlap'))
    for i, ot in enumerate(orbital_type):
        print('{:5}:  {:5} {:5.2f}'.format(i+1, ot[0], ot[1]))

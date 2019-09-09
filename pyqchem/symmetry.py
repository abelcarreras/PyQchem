from wfnsympy import WfnSympy
import numpy as np


def get_wf_symmetry(structure,
                    basis,
                    mo_coeff,
                    center=(0., 0., 0.),
                    orientation=(0., 0., 1.),
                    group='C2h'):

    """
    Simplified interface between pyqchem anf wfnsympy

    :param structure: molecular geometry (3xn) array
    :param basis: dictionary containing the basis
    :param mo_coeff: molecular orbitals coefficients
    :param center: center of the molecule
    :param orientation: unit vector perpendicular to the plane of the molecule
    :param group: point symmetry group
    :return molsym: wfnsympy object
    """


    alpha_mo_coeff = np.array(mo_coeff['alpha']).flatten().tolist()
    if 'beta' in mo_coeff:
        beta_mo_coeff = np.array(mo_coeff['beta']).flatten().tolist()
    else:
        beta_mo_coeff = None

    molsym = WfnSympy(coordinates=structure.get_coordinates(),
                      symbols=structure.get_atomic_elements(),
                      basis=basis,
                      center=center,
                      VAxis=orientation,
                      alpha_mo_coeff=alpha_mo_coeff,
                      beta_mo_coeff=beta_mo_coeff,
                      charge=structure.charge,
                      multiplicity=structure.multiplicity,
                      group=group)
    return molsym


def get_orbital_classification(fchk_data,
                               center=(0., 0., 0.),
                               orientation=(0., 0., 1.)):
    """
    Classify orbitals in sigma/pi in planar molecules

    :param fchk_data: dictionary containing fchk data
    :param orientation: unit vector perpendicular to the plane of the molecule
    :return: list of labels 'sigma'/'pi' according to the order of the orbitals
    """
    molsym = get_wf_symmetry(fchk_data['structure'],
                             fchk_data['basis'],
                             fchk_data['coefficients'],
                             center=center,
                             orientation=orientation)

    sh_index = molsym.SymLab.index('s_h')
    orbital_type_alpha = []
    for i, overlap in enumerate(molsym.mo_SOEVs_a[:, sh_index]):
        if overlap < 0:
            orbital_type_alpha.append(['pi', np.abs(overlap)])
        else:
            orbital_type_alpha.append(['sigma', np.abs(overlap)])

    if 'beta' in fchk_data['coefficients']:
        orbital_type_beta = []
        for i, overlap in enumerate(molsym.mo_SOEVs_b[:, sh_index]):
            if overlap < 0:
                orbital_type_beta.append(['pi', np.abs(overlap)])
            else:
                orbital_type_beta.append(['sigma', np.abs(overlap)])
        return orbital_type_alpha, orbital_type_beta
    else:
        return orbital_type_alpha


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
                  'basis': '6-31G'}

    qc_input = create_qchem_input(molecule, **parameters)

    _, _, fchk_data = get_output_from_qchem(qc_input, processors=4, force_recalculation=True,
                                            read_fchk=True)

    from pyqchem.file_io import build_fchk
    open('test.fchk', 'w').write(build_fchk(fchk_data))

    alpha_mo_coeff = set_zero_coefficients(fchk_data['basis'],
                                           fchk_data['coefficients']['alpha'],
                                           range(6, 12))
    fchk_data['coefficients']['alpha'] = alpha_mo_coeff

    structure = fchk_data['structure']
    print(structure.get_xyz())

    txt_fchk = build_fchk(fchk_data)
    open('test2.fchk', 'w').write(txt_fchk)

    z_dist = structure.get_coordinates()[0][2]

    orbital_type = get_orbital_classification(fchk_data, center=[0.0, 0.0, z_dist])

    print('  {:5} {:5} {:5}'.format('num', 'type', 'overlap'))
    for i, ot in enumerate(orbital_type):
        print('{:5}:  {:5} {:5.2f}'.format(i+1, ot[0], ot[1]))

from wfnsympy import WfnSympy
import numpy as np
from pyqchem.utils import get_occupied_electrons, reorder_coefficients


def get_wf_symmetry(structure,
                    basis,
                    mo_coeff,
                    center=None,
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

    alpha_mo_coeff = np.array(mo_coeff['alpha']).tolist()
    if 'beta' in mo_coeff:
        beta_mo_coeff = np.array(mo_coeff['beta']).tolist()
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
        overlap = overlap / molsym.mo_SOEVs_a[i, molsym.SymLab.index('E')]  # normalize
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




def get_state_symmetry(parsed_fchk,
                       rasci_states,
                       center=None,
                       orientation=(0, 0, 1),
                       group='C2h',
                       extra_print=False,
                       amplitude_cutoff=0.0
                       ):

    sym_states = {}
    for i, state in enumerate(rasci_states):

        structure = parsed_fchk['structure']
        total_orbitals = len(parsed_fchk['coefficients']['alpha'])

        if extra_print:
            print('HF/KS Ground state\n---------------------------')
            molsym = get_wf_symmetry(structure,
                                     parsed_fchk['basis'],
                                     parsed_fchk['coefficients'],
                                     center=center,
                                     orientation=orientation,
                                     group=group)

            molsym.print_alpha_mo_IRD()
            molsym.print_beta_mo_IRD()
            molsym.print_wf_mo_IRD()

            print('\nRASCI Exited state {}\n---------------------------'.format(i+1))

        # print(configurations)
        occupations_list = []
        for configuration in state['configurations']:
            if np.abs(configuration['amplitude']) < amplitude_cutoff:
                continue

            occupied_orbitals = get_occupied_electrons(configuration, structure)
            n_extra = total_orbitals - occupied_orbitals - len(configuration['alpha'])
            vector_alpha = [1] * occupied_orbitals + [int(c) for c in configuration['alpha']] + [0] * n_extra

            n_extra = total_orbitals - occupied_orbitals - len(configuration['beta'])
            vector_beta = [1] * occupied_orbitals + [int(c) for c in configuration['beta']] + [0] * n_extra

            occupations_list.append({'alpha': vector_alpha, 'beta': vector_beta})

            if extra_print:
                print(configuration['alpha'], configuration['beta'], configuration['amplitude'])



        state_symmetry_list = []
        for occupations in occupations_list:
            # print('occupations', occupations)

            reordered_coefficients = reorder_coefficients(occupations, parsed_fchk['coefficients'])

            molsym = get_wf_symmetry(structure,
                                     parsed_fchk['basis'],
                                     reordered_coefficients,
                                     center=center,
                                     orientation=orientation,
                                     group=group)

            if extra_print:
                molsym.print_alpha_mo_IRD()
                molsym.print_beta_mo_IRD()
                molsym.print_wf_mo_IRD()

            state_symmetry = {}
            for label, measure in zip(molsym.IRLab, molsym.wf_IRd):
                state_symmetry[label] = measure

            # return maximum symmetry and value
            if extra_print:
                print([molsym.IRLab[np.argmax(molsym.wf_IRd)], np.max(molsym.wf_IRd)])
            state_symmetry_list.append([molsym.IRLab[np.argmax(molsym.wf_IRd)], np.max(molsym.wf_IRd)])

        if extra_print:
            print(state_symmetry_list)

        # Make sure symmetry of all configurations is the same
        assert len(np.unique([a[0] for a in state_symmetry_list])) == 1

        symmetry_label = state_symmetry_list[0][0]
        average_measure = np.average([a[1] for a in state_symmetry_list])
        sym_states['state {}'.format(i+1)] = [symmetry_label, average_measure]
    return sym_states


if __name__ == '__main__':
    from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input
    from pyqchem.structure import Structure
    from pyqchem.file_io import build_fchk
    from pyqchem.utils import set_zero_coefficients

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

    mo_coeff = set_zero_coefficients(fchk_data['basis'],
                                           fchk_data['coefficients'],
                                           range(6, 12))
    fchk_data['coefficients'] = mo_coeff

    structure = fchk_data['structure']
    print(structure.get_xyz())

    txt_fchk = build_fchk(fchk_data)
    open('test2.fchk', 'w').write(txt_fchk)

    z_dist = structure.get_coordinates()[0][2]

    orbital_type = get_orbital_classification(fchk_data, center=[0.0, 0.0, z_dist])

    print('  {:5} {:5} {:5}'.format('num', 'type', 'overlap'))
    for i, ot in enumerate(orbital_type):
        print('{:5}:  {:5} {:5.2f}'.format(i+1, ot[0], ot[1]))

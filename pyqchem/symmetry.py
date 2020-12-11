from wfnsympy import WfnSympy
import numpy as np
from pyqchem.utils import get_occupied_electrons, reorder_coefficients
from pyqchem.utils import classify_diabatic_states_of_fragment, get_basis_functions_ranges_by_atoms
from pyqchem.utils import get_inertia, get_plane
from pyqchem import Structure


def get_wf_symmetry(structure,
                    basis,
                    mo_coeff,
                    center=None,
                    orientation=(0., 0., 1.),
                    orientation2=(0, 1, 0),
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

    alpha_occupancy = [0] * len(mo_coeff['alpha'])
    alpha_occupancy[:int(structure.alpha_electrons)] = [1] * structure.alpha_electrons

    if 'beta' in mo_coeff:
        beta_mo_coeff = np.array(mo_coeff['beta']).tolist()
        beta_occupancy = [0] * len(mo_coeff['beta'])
        beta_occupancy[:int(structure.alpha_electrons)] = [1] * structure.alpha_electrons
    else:
        beta_mo_coeff = None
        beta_occupancy = None


    molsym = WfnSympy(coordinates=structure.get_coordinates(),
                      symbols=structure.get_symbols(),
                      basis=basis,
                      center=center,
                      axis=orientation,
                      axis2=orientation2,
                      alpha_mo_coeff=alpha_mo_coeff,
                      beta_mo_coeff=beta_mo_coeff,
                      alpha_occupancy=alpha_occupancy,
                      beta_occupancy=beta_occupancy,
                      #charge=structure.charge,
                      #multiplicity=structure.multiplicity,
                      group=group)
    return molsym


def get_orbital_classification(fchk_data,
                               center=None,
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
                       orientation2=(0, 1, 0),
                       group='C2h',
                       extra_print=False,
                       amplitude_cutoff=0.0
                       ):
    """
    Determines electronic state symmetry (only for closed shell configurations)
    :param parsed_fchk: electronic structure
    :param rasci_states: parsed rasci output dictionary
    :param center: center of the molecule
    :param orientation:  main symmetry axis
    :param orientation2: secondary symmetry axis
    :param group: symmetry group label
    :param extra_print: if True prints wf & orbital symmetry analysis
    :param amplitude_cutoff: determined the configurations to take in account in the symmetry measure
    :return:
    """
    sym_states = {}
    for i, state in enumerate(rasci_states):

        # print('HF/KS Ground state',i ,'\n***********************')

        structure = parsed_fchk['structure']
        total_orbitals = len(parsed_fchk['coefficients']['alpha'])

        if extra_print:
            print('HF/KS Ground state\n---------------------------')
            molsym = get_wf_symmetry(structure,
                                     parsed_fchk['basis'],
                                     parsed_fchk['coefficients'],
                                     center=center,
                                     orientation=orientation,
                                     orientation2=orientation2,
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

            occupations_list.append(configuration['occupations'])

            # print('occupied', occupied_orbitals, total_orbitals)
            # print('sum', np.sum(occupations_list[-1]['alpha']), np.sum(occupations_list[-1]['beta']))
            # print('alpha', occupations_list[-1]['alpha'])
            # print('beta', occupations_list[-1]['beta'])

            # print(configuration['hole'],'|' , configuration['alpha'], configuration['beta'], '|', configuration['part'])
            if extra_print:
                print('occ:', occupations_list[-1])
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
                                     orientation2=orientation2,
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
        if len(np.unique([a[0] for a in state_symmetry_list])) == 1:
            symmetry_label = state_symmetry_list[0][0]
            average_measure = np.average([a[1] for a in state_symmetry_list])
            sym_states['state {}'.format(i+1)] = [symmetry_label, average_measure]
        else:
            sym_states['state {}'.format(i+1)] = ['Undef', 0]

    return sym_states

def _indices_from_ranges(ranges):
    indices = []
    for i, j in ranges:
        indices += list(range(i, j))

    return indices


def get_symmetry_le(electronic_structure, data_rasci, fragment_atoms=(0), tol=0.1, group='D2h'):
    # This only works for singlets on close shell calculations

    types = classify_diabatic_states_of_fragment(data_rasci['diabatization']['diabatic_states'], fragment_atoms, tol=0.1)
    functions_range = get_basis_functions_ranges_by_atoms(electronic_structure['basis'], atoms_range=fragment_atoms)
    coordinates_frag = np.array(electronic_structure['structure'].get_coordinates())[fragment_atoms]
    labels_frag = electronic_structure['structure'].get_symbols()[fragment_atoms]

    inertia_moments, inertia_axis = get_inertia(Structure(coordinates=coordinates_frag, symbols=labels_frag))
    center_frag, normal_frag = get_plane(coordinates_frag)

    # print(np.array(inertia_axis))
    # print(inertia_moments)
    # print(inertia_axis[0])
    # print(inertia_axis[1])
    # print(inertia_axis[2])

    if 'LE' in types:
        indices = [i for i, x in enumerate(types) if x == "LE"]
        symmetry_labels = []
        for index in indices:
            print('\nLE state found at diabat {}!'.format(index+1))

            coefficients = electronic_structure['nato_coefficients_multi'][index+1]
            occupation = electronic_structure['nato_occupancies_multi'][index+1]
            overlap_matrix = np.array(electronic_structure['overlap'])

            functions_indices = _indices_from_ranges(functions_range)
            overlap_matrix = np.array([ovm[functions_indices] for ovm in overlap_matrix[functions_indices]])

            c_alpha = np.array(coefficients['alpha'])
            oc_alpha = np.array(occupation['alpha'])

            orbitals = []
            alpha = 0
            beta = 0

            coefficients_new = {'alpha': [], 'beta': []}
            for i, (va, o) in enumerate(zip(c_alpha, oc_alpha)):
                va_frag = va[functions_indices]
                frag_dot = np.dot(va_frag, np.dot(overlap_matrix, va_frag))

                if np.abs(frag_dot - 0.5) > tol:
                    if frag_dot > 0.5:
                        orbitals.append((i, np.round(o)))

                        orbital = np.zeros(len(c_alpha))
                        for j, val in zip(functions_indices, va_frag):
                            orbital[j] = val/np.sqrt(frag_dot)

                        if np.abs(o - 1) < tol and o < 2 - tol:
                            if alpha == beta:
                                alpha += 1
                                coefficients_new['alpha'].append(list(orbital))
                            else:
                                beta += 1
                                coefficients_new['beta'].append(list(orbital))
                    else:
                        pass

            if alpha == beta:
                multiplicity = 1
            else:
                multiplicity = 3

            n_electrons = alpha + beta

            try:
                # This only works for singlets
                molsym = get_wf_symmetry(electronic_structure['structure'],
                                         electronic_structure['basis'],
                                         coefficients_new,
                                         center=center_frag,
                                         orientation=inertia_axis[0],
                                         orientation2=inertia_axis[1],
                                         group=group
                                         )

                #molsym.print_alpha_mo_IRD()
                #molsym.print_beta_mo_IRD()
                molsym.print_wf_mo_IRD()
                symmetry_labels.append(molsym.IRLab[np.argmax(molsym.wf_IRd)])
            except:
                pass

        return symmetry_labels

    return None

if __name__ == '__main__':
    from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input
    from pyqchem.structure import Structure
    from pyqchem.file_io import build_fchk
    from pyqchem.utils import _set_zero_to_coefficients

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
                         symbols=symbols,
                         charge=0)

    parameters = {'jobtype': 'sp',
                  'exchange': 'hf',
                  'basis': '6-31G'}

    qc_input = create_qchem_input(molecule, **parameters)

    _, _, fchk_data = get_output_from_qchem(qc_input, processors=4, force_recalculation=True,
                                            read_fchk=True)

    from pyqchem.file_io import build_fchk
    open('test.fchk', 'w').write(build_fchk(fchk_data))

    mo_coeff = _set_zero_to_coefficients(fchk_data['basis'],
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

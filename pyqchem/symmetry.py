from pyqchem.utils import get_occupied_electrons, reorder_coefficients
from pyqchem.utils import classify_diabatic_states_of_fragment, get_basis_functions_ranges_by_atoms
from pyqchem.utils import get_inertia, get_plane
from pyqchem import Structure
import numpy as np


def get_orbitals(structure,
                 basis,
                 mo_coeff,
                 orbital_numbers=None):

    from posym.tools import get_basis_set, build_orbital
    from posym import SymmetryMolecule

    if orbital_numbers is None:
        orbital_numbers = [i+1 for i in range(len(mo_coeff))]

    basis_set = get_basis_set(structure.get_coordinates(), basis)

    orbitals = []
    for i, orbital_coeff in enumerate(mo_coeff):
        if i+1 in orbital_numbers:
            orbitals.append(build_orbital(basis_set, orbital_coeff))

    return orbitals


def get_orbitals_symmetry(structure,
                          basis,
                          mo_coeff,
                          center=None,
                          orientation=None,
                          group='C2h',
                          orbital_numbers=None):

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

    from posym import SymmetryGaussianLinear, SymmetryMolecule


    if orbital_numbers is None:
        orbital_numbers = [i+1 for i in range(len(mo_coeff))]

    if orientation is None:
        sm = SymmetryMolecule(group=group,
                              coordinates=structure.get_coordinates(),
                              symbols=structure.get_symbols(),
                              )

        orientation = sm.orientation_angles

    orbitals = get_orbitals(structure, basis, mo_coeff, orbital_numbers=orbital_numbers)

    sym_orbitals = []
    for orbital in orbitals:
        sym_orbital = SymmetryGaussianLinear(group, orbital, center=center, orientation_angles=orientation)
        sym_orbitals.append(sym_orbital)

    return sym_orbitals


def get_wf_symmetry(fchk_data,
                    alpha_electrons,
                    beta_electrons,
                    center=None,
                    orientation=None,
                    active_orbitals=0,
                    group='C2h'):
    """
    get symmetry of single determinat wave function

    :param fchk_data: electronic structure
    :param alpha_electrons: number of alpha electrons
    :param beta_electrons: number of beta electrons
    :param center: center of the group
    :param orientation: orientation angles
    :param active_orbitals: number of double occupied orbitals to use for symmetry calculation
    :param group: point group
    :return: symmetry
    """

    from posym import SymmetrySingleDeterminant

    base_electrons = np.min([alpha_electrons, beta_electrons])

    if active_orbitals < 1 and alpha_electrons- beta_electrons == 0:
        active_orbitals = 1

    alpha_numbers = list(range(base_electrons - active_orbitals + 1, alpha_electrons + 1))
    beta_numbers = list(range(base_electrons - active_orbitals + 1, beta_electrons + 1))

    coef_alpha = fchk_data['coefficients']['alpha']
    coef_beta = fchk_data['coefficients']['beta'] if 'beta' in fchk_data['coefficients'] else coef_alpha

    orbitals_alpha = get_orbitals(fchk_data['structure'],
                                  fchk_data['basis'],
                                  coef_alpha,
                                  orbital_numbers=alpha_numbers)

    orbitals_beta = get_orbitals(fchk_data['structure'],
                                 fchk_data['basis'],
                                 coef_beta,
                                 orbital_numbers=beta_numbers)

    symmetry = SymmetrySingleDeterminant(group=group,
                                         alpha_orbitals=orbitals_alpha,
                                         beta_orbitals=orbitals_beta,
                                         center=center,
                                         orientation_angles=orientation)

    return symmetry


def get_orbital_classification(fchk_data,
                               center=None,
                               orientation=(0, 0, 0)):
    """
    Classify orbitals in sigma/pi in planar molecules

    :param fchk_data: dictionary containing fchk data
    :param orientation: unit vector perpendicular to the plane of the molecule
    :return: list of labels 'sigma'/'pi' according to the order of the orbitals
    """

    orbitals_sym = get_orbitals_symmetry(fchk_data['structure'],
                                         fchk_data['basis'],
                                         fchk_data['coefficients']['alpha'],
                                         center=center,
                                         orientation=orientation)

    orbital_type_alpha = []
    for orb_sym in orbitals_sym:
        # use inversion operation to classify
        trace = orb_sym.get_reduced_op_representation()['sh']
        if trace < 0:
            orbital_type_alpha.append([' pi', np.abs(trace)])
        else:
            orbital_type_alpha.append([' sigma', np.abs(trace)])

    if 'beta' in fchk_data['coefficients']:

        orbitals_sym = get_orbitals_symmetry(fchk_data['structure'],
                                             fchk_data['basis'],
                                             fchk_data['coefficients']['beta'],
                                             center=center,
                                             orientation=orientation)

        orbital_type_beta = []
        for orb_sym in orbitals_sym:
            # use inversion operation to classify
            trace = orb_sym.get_reduced_op_representation()['sh']
            if trace < 0:
                orbital_type_beta.append([' pi', np.abs(trace)])
            else:
                orbital_type_beta.append([' sigma', np.abs(trace)])

        return orbital_type_alpha, orbital_type_beta
    else:
        return orbital_type_alpha


def get_state_symmetry(parsed_fchk,
                       rasci_states,
                       center=None,
                       group='C2h',
                       states_numbers=None):
    """
    Determine the symmetry of electronic states of multi-configurational states

    :param parsed_fchk: electronic structure
    :param rasci_states: parsed rasci output dictionary
    :param center: center of the molecule
    :param group: symmetry group label
    :param states_numbers: list of the states to be analyzed
    :return: states symmetry
    """

    from posym import SymmetryMultiDeterminant

    structure = parsed_fchk['structure']
    # total_orbitals = len(parsed_fchk['coefficients']['alpha'])

    total_states = len(rasci_states)
    if states_numbers is None:
        states_numbers = [i+1 for i in range(total_states)]

    def get_relevant_occupation_range(rasci_states):

        occupations_mat = []
        for istate, state in enumerate(rasci_states):
            for configuration in state['configurations']:
                # print('amplitude: {} '.format(configuration['occupations']['alpha']))
                occupations_mat.append(configuration['occupations']['alpha'])
                occupations_mat.append(configuration['occupations']['beta'])

        n_conf = len(occupations_mat)
        occupations_sum = np.sum(occupations_mat, axis=0)

        inferior = 0
        superior = n_conf

        for i, c in enumerate(occupations_sum):
            if c == n_conf:
                inferior = i+2
            if np.sum(occupations_sum[i:]) == 0:
                superior = i+1
                break
        return list(range(inferior, superior))

    relevant_range = get_relevant_occupation_range(rasci_states)

    orbitals = get_orbitals(structure,
                            parsed_fchk['basis'],
                            parsed_fchk['coefficients']['alpha'],
                            orbital_numbers=relevant_range)

    states_symmetry = []
    for istate, state in enumerate(rasci_states):

        if istate+1 not in states_numbers:
            continue

        inf_index = relevant_range[0]-1
        sup_index = relevant_range[-1]

        conf_reduced = []
        for configuration in state['configurations']:
            conf_reduced.append({'amplitude': configuration['amplitude'],
                                 'occupations': {'alpha': configuration['occupations']['alpha'][inf_index:sup_index],
                                                 'beta': configuration['occupations']['beta'][inf_index:sup_index]}
                                 })

        state_symmetry = SymmetryMultiDeterminant(group=group,
                                                  orbitals=orbitals,
                                                  center=center,
                                                  configurations=conf_reduced)

        states_symmetry.append(state_symmetry)

    return states_symmetry



# Alternative function that make use of WFNSYM - to be deprecated in the future

def get_state_symmetry_wfnsym(parsed_fchk,
                              rasci_states,
                              center=None,
                              orientation=(0, 0, 1),
                              orientation2=(0, 1, 0),
                              group='C2h',
                              extra_print=False,
                              amplitude_cutoff=0.0,
                              check_consistency=True,
                              full_vector=False):
    """
    Determine electronic state symmetry (only for closed shell configurations)

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
            molsym = get_wf_symmetry_wfsym(structure,
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

        occupations_list = []
        for configuration in state['configurations']:
            if np.abs(configuration['amplitude']) < amplitude_cutoff:
                continue

            occupations_list.append(configuration['occupations'])

            if extra_print:
                print('occ:', occupations_list[-1])
                print(configuration['occupations']['alpha'], configuration['occupations']['beta'], configuration['amplitude'])

        state_symmetry_list = []
        states_symmetry_list = []
        for occupations in occupations_list:

            molsym = get_wf_symmetry_wfsym(structure,
                                     parsed_fchk['basis'],
                                     parsed_fchk['coefficients'],
                                     center=center,
                                     orientation=orientation,
                                     orientation2=orientation2,
                                     occupancy=occupations,
                                     group=group)

            if extra_print:
                molsym.print_alpha_mo_IRD()
                molsym.print_beta_mo_IRD()
                molsym.print_wf_mo_IRD()
                print('center: ', molsym.center)
                print('orientation axis', molsym.axis, molsym.axis2)

            state_symmetry = {}
            for label, measure in zip(molsym.IRLab, molsym.wf_IRd):
                state_symmetry[label] = measure

            # return maximum symmetry and value
            if extra_print:
                print([molsym.IRLab[np.argmax(molsym.wf_IRd)], np.max(molsym.wf_IRd)])
            state_symmetry_list.append([molsym.IRLab[np.argmax(molsym.wf_IRd)], np.max(molsym.wf_IRd)])
            states_symmetry_list.append([molsym.IRLab, molsym.wf_IRd, configuration['amplitude']])

        if extra_print:
            print(state_symmetry_list)

        # Make sure symmetry of all configurations is the same
        if check_consistency:
            if len(np.unique([a[0] for a in state_symmetry_list])) == 1:
                symmetry_label = state_symmetry_list[0][0]
                average_measure = np.average([a[1] for a in state_symmetry_list])
                sym_states['state {}'.format(i+1)] = [symmetry_label, average_measure]
            else:
                sym_states['state {}'.format(i+1)] = ['Undef', 0]
        else:
            if not full_vector:
                sym_states['state {}'.format(i+1)] = state_symmetry_list[np.argmax([a[1] for a in state_symmetry_list])]
            else:
                sym_states['state {}'.format(i+1)] = states_symmetry_list[np.argmax([a[2]**2 for a in states_symmetry_list])]

                def average(states_symmetry_list):
                    amplitudes2 = [a[2]**2 for a in states_symmetry_list]
                    measure_average = np.average([a[1] for a in states_symmetry_list], axis=0, weights=amplitudes2)

                    return [states_symmetry_list[0][0], measure_average, np.average(amplitudes2)]

                # sym_states['state {}'.format(i+1)] = average(states_symmetry_list)
                sym_states['state {}'.format(i+1)] = states_symmetry_list[np.argmax([a[2]**2 for a in states_symmetry_list])]

    return [sym_states['state {}'.format(i+1)] for i in range(len(sym_states))]


def get_wf_symmetry_wfsym(structure,
                          basis,
                          mo_coeff,
                          center=None,
                          orientation=(0., 0., 1.),
                          orientation2=(0, 1, 0),
                          group='C2h',
                          occupancy=None):

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
    from wfnsympy import WfnSympy

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

    if occupancy is not None:
        alpha_occupancy = occupancy['alpha']
        if 'beta' in occupancy:
            beta_occupancy = occupancy['beta']

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


def get_orbital_classification_wfsym(fchk_data,
                                     center=None,
                                     orientation=(0., 0., 1.)):
    """
    Classify orbitals in sigma/pi in planar molecules

    :param fchk_data: dictionary containing fchk data
    :param orientation: unit vector perpendicular to the plane of the molecule
    :return: list of labels 'sigma'/'pi' according to the order of the orbitals
    """
    molsym = get_wf_symmetry_wfsym(fchk_data['structure'],
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
                molsym = get_wf_symmetry_wfsym(electronic_structure['structure'],
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

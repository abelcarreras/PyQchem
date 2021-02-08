import numpy as np

from pyqchem.symmetry import get_wf_symmetry
from pyqchem.utils import get_plane, _set_zero_to_coefficients
from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.structure import Structure
from pyqchem.file_io import build_fchk
from copy import deepcopy
from pyqchem.parsers.parser_rasci import parser_rasci
from pyqchem.utils import reorder_coefficients
from pyqchem.symmetry import get_state_symmetry


def inverse_cummulative(vector):
    vector = np.array(vector)

    for i in range(1, len(vector)):
        vector[-(i+1)] -= np.sum(vector[-i:])

    return vector


sym_matrix = {'Ag': np.diag([1, 1]),
              'Au': np.diag([1, -1]),
              'Bg': np.diag([-1,  1]),
              'Bu': np.diag([-1, -1])}


def get_location(occupancies, coefficients, overlap_matrix, error=0.2):
    # For singlets only
    alpha_list = []
    beta_list = []
    mul = 1
    i = 0
    for orbital in coefficients['alpha']:
        occ = occupancies['alpha'][i]
        dot = np.dot(np.dot(orbital, overlap_matrix), orbital)
        print('dot->',i+1,  dot)
        if dot < error:
            alpha_list.append(0)
            beta_list.append(0)
        else:
            if occ < error:
                # print('0')
                alpha_list.append(0)
                beta_list.append(0)
            elif 1+error > occ > error:
                # print('1')
                if mul == 1:
                    alpha_list.append(1)
                    beta_list.append(0)

                    mul += 1
                else:
                    alpha_list.append(0)
                    beta_list.append(1)
                    mul -= 1
            elif occ > 2-error:
                # print('2')
                alpha_list.append(1)
                beta_list.append(1)
            i +=1

    print(alpha_list)
    print(beta_list)

    return {'alpha': alpha_list, 'beta': beta_list}


def orbital_localization2(coefficients, overlap_matrix):
    n = 1
    for a_orbital, b_orbital in zip(coefficients['alpha'], coefficients['beta']):
        a_dot = np.dot(np.dot(a_orbital, overlap_matrix), a_orbital)
        b_dot = np.dot(np.dot(b_orbital, overlap_matrix), b_orbital)

        print('dot', n, a_dot, b_dot)
        n += 1
    return

def get_state_classification(structure,
                             basis,
                             nato_coeff,
                             nato_occup,
                             center,
                             orientation,
                             orbitals,
                             overlap_matrix):

    occupations = get_location(nato_occup, nato_coeff, overlap_matrix)

    reordered_coefficients = reorder_coefficients(occupations, nato_coeff)

    orbital_localization2(reordered_coefficients,overlap_matrix)

    molsym = get_wf_symmetry(structure,
                             basis,
                             reordered_coefficients,
                             center=center,
                             orientation=orientation,
                             group='D2h')

    # molsym.print_alpha_mo_IRD()
    # molsym.print_beta_mo_IRD()
    molsym.print_wf_mo_IRD()
    print('  wf' + '  '.join(['{:7.3f}'.format(s) for s in molsym.wf_IRd]))

    print('\nAlpha MOs: Irred. Rep. Decomposition')
    print('     ' + '  '.join(['{:^7}'.format(s) for s in molsym.IRLab]))
    for i, line in zip(orbitals, molsym.mo_IRd_a[orbitals]):
        print('{:4d}'.format(i + 1) + '  '.join(['{:7.3f}'.format(s) for s in line]))

    print('\nBeta MOs: Irred. Rep. Decomposition')
    print('     ' + '  '.join(['{:^7}'.format(s) for s in molsym.IRLab]))
    for i, line in zip(orbitals, molsym.mo_IRd_b[orbitals]):
        print('{:4d}'.format(i + 1) + '  '.join(['{:7.3f}'.format(s) for s in line]))

    print('occupations')
    print('alpha', np.array(nato_occup['alpha'])[orbitals])
    print('beta ', np.array(nato_occup['beta'])[orbitals])








ethene = [[0.0,  0.0000,   0.65750],
          [0.0,  0.0000,  -0.65750],
          [0.0,  0.92281,  1.22792],
          [0.0, -0.92281,  1.22792],
          [0.0, -0.92281, -1.22792],
          [0.0,  0.92281, -1.22792]]

symbols = ['C', 'C', 'H', 'H', 'H', 'H']


# create molecule
molecule = Structure(coordinates=ethene,
                     symbols=symbols,
                     charge=0,
                     multiplicity=1)


print('Optimized monomer structure')
print(molecule)

# Build dimer from monomer
coor_monomer2 = np.array(ethene)

from kimonet.utils.rotation import rotate_vector
coor_monomer2 = np.array([rotate_vector(coor, [0.0, 0.0, 0.0]) for coor in coor_monomer2])


coor_monomer2[:, 0] += 4.0  # monomers separation
coor_monomer2[:, 1] += 0.0  # monomers separation
coor_monomer2[:, 2] += 0.0  # monomers separation

coordinates = ethene + coor_monomer2.tolist()
symbols_dimer = symbols * 2

dimer = Structure(coordinates=coordinates,
                  symbols=symbols_dimer,
                  charge=0,
                  multiplicity=1)

print(dimer.get_xyz())

diabatization_scheme = [
                        # {'method': 'ER', 'states': [1, 2, 3, 4]},  # in order with respect to selected states
                        {'method': 'DQ', 'states': [1, 2, 3, 4], 'parameters': 0.0}
                        ]

# create Q-Chem input
qc_input = QchemInput(dimer,
                      jobtype='sp',
                      exchange='hf',
                      correlation='rasci',
                      ras_act=6,
                      ras_elec=4,
                      basis='sto-3g',
                      ras_spin_mult=1,  # singlets only
                      ras_roots=8,
                      ras_do_hole=False,
                      ras_do_part=False,
                      ras_diabatization_states=[3, 4, 5, 6],  # adiabatic states selected for diabatization
                      ras_diabatization_scheme=diabatization_scheme,
                      )

print(qc_input.get_txt())


# get data from Q-Chem calculation
output, electronic_structure = get_output_from_qchem(qc_input,
                                                          processors=4,
                                                          force_recalculation=False,
                                                          read_fchk=True,
                                                          parser=parser_rasci,
                                                          store_full_output=True)


def is_equal(conf1, conf2):
    lv = []
    for prop in ['alpha', 'beta', 'hole', 'part']:
        lv.append(conf1[prop] == conf2[prop])
    return np.all(lv)


def unify_configurations(configurations_list):
    new_configurations= []

    for conf in configurations_list:
        skip = False
        for n_conf in new_configurations:
            if is_equal(n_conf, conf):
                n_conf['amplitude'] += conf['amplitude']
                skip = True
        if not skip:
            new_configurations.append(conf)

    return new_configurations


print(output)
print(output['excited_states'])


print('++++++++')

print(len(electronic_structure['nato_coefficients_multi']))




# Just for testing
# electronic_structure['coefficients'] = electronic_structure['nato_coefficients']
# electronic_structure['mo_energies'] = electronic_structure['nato_occupancies']

print(len(electronic_structure['nato_coefficients_multi']))
print(len(electronic_structure['nato_occupancies_multi']))
print(electronic_structure['nato_occupancies_multi'])

# store original fchk info in file
es = deepcopy(electronic_structure)
es['coefficients'] = electronic_structure['nato_coefficients_multi'][1]
es['mo_energies'] = electronic_structure['nato_occupancies_multi'][1]

open('test.fchk', 'w').write(build_fchk(es))


# print(electronic_structure['nato_coefficients'])

range_f1 = range(0, 6)
range_f2 = range(6, 12)



mo_coeff_f1 = _set_zero_to_coefficients(electronic_structure['basis'],
                                        electronic_structure['nato_coefficients_multi'][1],
                                        range_f2)



# get symmetry classification
electronic_structure['coefficients'] = mo_coeff_f1
electronic_structure['mo_energies'] = electronic_structure['nato_occupancies']

# save test fchk file with new coefficients
open('test_f1.fchk', 'w').write(build_fchk(electronic_structure))


# get plane from coordinates
coordinates_f1 = np.array(electronic_structure['structure'].get_coordinates())[range_f1]
center_f1, normal_f1 = get_plane(coordinates_f1)


mulliken = output['diabatization']['mulliken_analysis'][1]
a = np.sum(np.array(mulliken['attach'])[range_f1])
print(a)
print(electronic_structure['overlap'])

print(mulliken)

def orbital_localization(coefficients, overlap_matrix):
    for a_orbital, b_orbital in zip(coefficients['alpha'], coefficients['beta']):
        a_dot = np.dot(np.dot(a_orbital, overlap_matrix), a_orbital)
        b_dot = np.dot(np.dot(b_orbital, overlap_matrix), b_orbital)

        print('dot', a_dot, b_dot)
    return

# orbital_localization(mo_coeff_f1, electronic_structure['overlap'])


# get classified orbitals
get_state_classification(electronic_structure['structure'],
                         electronic_structure['basis'],
                         mo_coeff_f1,
                         electronic_structure['nato_occupancies_multi'][1],
                         center=center_f1,
                         orientation=normal_f1,
                         orbitals=[10, 11, 12, 13, 14, 15],
                         overlap_matrix=electronic_structure['overlap'])


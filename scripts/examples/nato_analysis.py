import numpy as np

from pyqchem.symmetry import get_wf_symmetry
from pyqchem.utils import get_plane, set_zero_coefficients
from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.structure import Structure
from pyqchem.file_io import build_fchk
from copy import deepcopy
from pyqchem.parsers.parser_rasci import rasci as rasci_parser
from pyqchem.utils import reorder_coefficients

def inverse_cummulative(vector):
    vector = np.array(vector)

    for i in range(1, len(vector)):
        vector[-(i+1)] -= np.sum(vector[-i:])

    return vector


sym_matrix = {'Ag': np.diag([1, 1]),
              'Au': np.diag([1, -1]),
              'Bg': np.diag([-1,  1]),
              'Bu': np.diag([-1, -1])}


# Define custom classification function
def get_state_classification(parsed_fchk,
                             center=None,
                             orientation=(0, 0, 1),
                             configurations=None,
                             occupied=0
                             ):

    structure = parsed_fchk['structure']

    occupations_list = []
    for conf in configurations:
        vector_alpha = [1 for _ in range(occupied)] + [int(c) for c in conf['alpha']]
        vector_beta = [1 for _ in range(occupied)] + [int(c) for c in conf['beta']]
        occupations_list.append({'alpha': vector_alpha, 'beta': vector_beta})
        print(conf['alpha'], conf['beta'])
    # exit()

    molsym = get_wf_symmetry(structure,
                             parsed_fchk['basis'],
                             parsed_fchk['nato_coefficients'],
                             center=center,
                             orientation=orientation,
                             group='C2h')

    molsym.print_alpha_mo_IRD()
    molsym.print_beta_mo_IRD()
    molsym.print_wf_mo_IRD()

    for occupations in occupations_list:
        print(occupations)
        cumm_symmat = np.identity(2)
        print(occupations['alpha'])

        for i, occup in enumerate(occupations['alpha']):
            if occup == 1:
                symmat = np.zeros_like(cumm_symmat)
                for symlab, contribution in zip(molsym.IRLab, molsym.mo_IRd_a[i]):
                    symmat += np.dot(contribution, sym_matrix[symlab])

                cumm_symmat = np.dot(cumm_symmat, symmat)
                # print('alpha: ', i+1, cumm_symmat)
        for i, occup in enumerate(occupations['beta']):
            if occup == 1:

                symmat = np.zeros_like(cumm_symmat)
                for symlab, contribution in zip(molsym.IRLab, molsym.mo_IRd_b[i]):
                    symmat += np.dot(contribution, sym_matrix[symlab])

                cumm_symmat = np.dot(cumm_symmat, symmat)
                # print('beta: ', i+1, cumm_symmat)

        print('fin', cumm_symmat)

    print('------')


def get_state_symmetry(parsed_fchk,
                       center=None,
                       orientation=(0, 0, 1),
                       configurations=None,
                       occupied=0
                       ):

    structure = parsed_fchk['structure']
    n_orbitals = len(parsed_fchk['coefficients']['alpha'])

    print(configurations)
    occupations_list = []
    for conf in configurations:
        n_extra = n_orbitals - occupied - len(conf['alpha'])
        vector_alpha = [1 for _ in range(occupied)] + [int(c) for c in conf['alpha']] + [0] * n_extra
        vector_beta = [1 for _ in range(occupied)] + [int(c) for c in conf['beta']] + [0] * n_extra
        occupations_list.append({'alpha': vector_alpha, 'beta': vector_beta})
        print(conf['alpha'], conf['beta'], )

    for occupations in occupations_list:
        print('occupations', occupations)

        reordered_coefficients = reorder_coefficients(occupations, parsed_fchk['coefficients'])

        molsym = get_wf_symmetry(structure,
                                 parsed_fchk['basis'],
                                 reordered_coefficients,
                                 center=center,
                                 orientation=orientation,
                                 group='C2h')

        molsym.print_alpha_mo_IRD()
        molsym.print_beta_mo_IRD()
        molsym.print_wf_mo_IRD()

        state_symmetry = {}
        for label, measure in zip(molsym.IRLab, molsym.wf_IRd):
            state_symmetry[label] = measure

        # return maximum symmetry and value
        return [molsym.IRLab[np.argmax(molsym.wf_IRd)], np.max(molsym.wf_IRd)]


ethene = [[0.0,  0.0000,   0.65750],
          [0.0,  0.0000,  -0.65750],
          [0.0,  0.92281,  1.22792],
          [0.0, -0.92281,  1.22792],
          [0.0, -0.92281, -1.22792],
          [0.0,  0.92281, -1.22792]]

symbols = ['C', 'C', 'H', 'H', 'H', 'H']

range_f1 = range(0, 6)

n_state = 1  # First excited state

# create molecule
molecule = Structure(coordinates=ethene,
                     atomic_elements=symbols,
                     charge=0,
                     multiplicity=1)

# create Q-Chem input
qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      correlation='rasci',
                      ras_act=6,
                      ras_elec=6,
                      basis='sto-3g',
                      ras_roots=10,
                      ras_do_hole=False,
                      ras_do_part=False,
                      ras_natorb_state=n_state)

print(qc_input.get_txt())


# get data from Q-Chem calculation
output, err, electronic_structure = get_output_from_qchem(qc_input,
                                                          processors=4,
                                                          force_recalculation=False,
                                                          read_fchk=True,
                                                          parser=rasci_parser,
                                                          store_full_output=True)

print(output)
print(output['excited states rasci'])

# Just for testing
print(electronic_structure['nato_coefficients'])
# electronic_structure['coefficients'] = electronic_structure['nato_coefficients']
# electronic_structure['mo_energies'] = electronic_structure['nato_occupancies']

# store original fchk info in file
open('test.fchk', 'w').write(build_fchk(electronic_structure))


# print(electronic_structure['nato_coefficients'])


# get plane from coordinates
coordinates_f1 = np.array(electronic_structure['structure'].get_coordinates())[range_f1]
center_f1, normal_f1 = get_plane(coordinates_f1)


# get classified orbitals
# get_state_classification(electronic_structure,
#                          center=center_f1,
#                          orientation=normal_f1,
#                          configurations=output['excited states rasci'][n_state]['amplitudes'],
#                          occupied=qc_input._ras_occ)

wfsym = get_state_symmetry(electronic_structure,
                   center=center_f1,
                   orientation=normal_f1,
                   configurations=output['excited states rasci'][n_state]['amplitudes'],
                   occupied=qc_input._ras_occ)

print(wfsym)

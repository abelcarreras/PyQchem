import numpy as np

from pyqchem.symmetry import set_zero_coefficients, get_wf_symmetry, get_orbital_classification
from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input
from pyqchem.structure import Structure
from pyqchem.file_io import build_fchk


dimer_ethene = [[0.0,  0.0000,   0.65750],
                [0.0,  0.0000,  -0.65750],
                [0.0,  0.92281,  1.22792],
                [0.0, -0.92281,  1.22792],
                [0.0, -0.92281, -1.22792],
                [0.0,  0.92281, -1.22792],
                [3.7,  0.00000,  0.65750],
                [3.7,  0.00000, -0.65750],
                [3.7,  0.92281,  1.22792],
                [3.7, -0.92281,  1.22792],
                [3.7, -0.92281, -1.22792],
                [3.7,  0.92281, -1.22792]]

symbols = ['C', 'C', 'H', 'H', 'H', 'H', 'C', 'C', 'H', 'H', 'H', 'H']

# create molecule
molecule = Structure(coordinates=dimer_ethene,
                     atomic_elements=symbols,
                     charge=0,
                     multiplicity=1)

# create Q-Chem input
qc_input = create_qchem_input(molecule,
                              jobtype='sp',
                              exchange='hf',
                              basis='sto-3g')


# get data from Q-Chem calculation
output, err, parsed_fchk = get_output_from_qchem(qc_input,
                                                 processors=4,
                                                 force_recalculation=False,
                                                 read_fchk=True,
                                                 fchk_only=True)

# get partial wf localized in fragment 1
alpha_mo_coeff_f1 = set_zero_coefficients(parsed_fchk['basis'],
                                          parsed_fchk['coefficients']['alpha'],
                                          range(0, 6))

# get partial wf localized in fragment 2
alpha_mo_coeff_f2 = set_zero_coefficients(parsed_fchk['basis'],
                                          parsed_fchk['coefficients']['alpha'],
                                          range(6, 12))


# Define custom classification criteria
def get_custom_orbital_classification(parsed_fchk,
                                      center=(0.0, 0.0, -1.85),
                                      orientation=(0, 0, 1)
                                      ):

    molsym = get_wf_symmetry(parsed_fchk['structure'],
                             parsed_fchk['basis'],
                             parsed_fchk['coefficients'],
                             center=center,
                             orientation=orientation)

    sh_index = molsym.SymLab.index('i')  # operation used to separate orbitals
    orbital_type = []
    for i, overlap in enumerate(molsym.mo_SOEVs_a[:, sh_index]):
        if overlap < 0:
            orbital_type.append([' YES', np.abs(overlap)])
        else:
            orbital_type.append([' NOO', np.abs(overlap)])
    return orbital_type


parsed_fchk['coefficients']['alpha'] = alpha_mo_coeff_f1
open('test_f1.fchk', 'w').write(build_fchk(parsed_fchk))
orbital_type_f1 = get_custom_orbital_classification(parsed_fchk,
                                                    center=[0.0, 0.0, -1.85],
                                                    orientation=[0, 0, 1])

parsed_fchk['coefficients']['alpha'] = alpha_mo_coeff_f2
open('test_f2.fchk', 'w').write(build_fchk(parsed_fchk))
orbital_type_f2 = get_custom_orbital_classification(parsed_fchk,
                                                    center=[0.0, 0.0, 1.85],
                                                    orientation=[0, 0, 1])

frontier_orbitals = [11, 12, 13,  14, 15, 16, 17, 18, 19]

print('Inversion center?')
print('index  frag1  frag2')
for i in frontier_orbitals:
    print(' {}    {}   {}'.format(i+1, orbital_type_f1[i][0], orbital_type_f2[i][0]))


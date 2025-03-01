# Classify (quasi) molecular orbitals of the individual monomers from a dimer calculation
from pyqchem.symmetry import get_orbitals_symmetry
from pyqchem.utils import get_plane, crop_electronic_structure
from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input
from pyqchem.structure import Structure
from pyqchem.file_io import build_fchk
import numpy as np


# Define custom classification function
def get_custom_orbital_classification(parsed_fchk,
                                      center=None,
                                      orientation=(0, 0, 1),
                                      orbital_numbers=None
                                      ):

    n_orbitals = parsed_fchk['coefficients']['alpha']
    if orbital_numbers is None:
        orbital_numbers = list(np.range(len(n_orbitals)))

    orbitals_sym = get_orbitals_symmetry(parsed_fchk['structure'],
                                         parsed_fchk['basis'],
                                         parsed_fchk['coefficients']['alpha'],
                                         center=center,
                                         orientation=orientation,
                                         orbital_numbers=orbital_numbers)

    orbital_type = []
    for orb_sym in orbitals_sym:
        # use inversion operation to classify
        trace = orb_sym.get_reduced_op_representation()['i']
        if trace < 0:
            orbital_type.append([' NOO', np.abs(trace)])
        else:
            orbital_type.append([' YES', np.abs(trace)])
    return orbital_type


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

range_f1 = range(0, 6)
range_f2 = range(6, 12)


# create molecule
molecule = Structure(coordinates=dimer_ethene,
                     symbols=symbols,
                     charge=0,
                     multiplicity=1)

# create Q-Chem input
qc_input = create_qchem_input(molecule,
                              jobtype='sp',
                              exchange='hf',
                              basis='6-31G')

# get data from Q-Chem calculation
_, electronic_structure = get_output_from_qchem(qc_input,
                                                processors=4,
                                                force_recalculation=False,
                                                return_electronic_structure=True)

# store original fchk info in file
open('test.fchk', 'w').write(build_fchk(electronic_structure))

# get symmetry classification
electronic_structure_f1 = crop_electronic_structure(electronic_structure, range_f1)

# save test fchk file with new coefficients
open('test_f1.fchk', 'w').write(build_fchk(electronic_structure_f1))

# range of orbitals to show
frontier_orbitals = [12, 13,  14, 15, 16, 17, 18, 19, 20]

# get plane from coordinates
coordinates_f1 = electronic_structure['structure'].get_coordinates(fragment=range_f1)
center_f1, normal_f1 = get_plane(coordinates_f1)

# get classified orbitals
orbital_type_f1 = get_custom_orbital_classification(electronic_structure_f1,
                                                    center=center_f1,
                                                    orientation=normal_f1,
                                                    orbital_numbers=frontier_orbitals)

# get plane from coordinates
coordinates_f2 = electronic_structure['structure'].get_coordinates(fragment=range_f2)
center_f2, normal_f2 = get_plane(coordinates_f2)

electronic_structure_f2 = crop_electronic_structure(electronic_structure, range_f2)

# save test fchk file with new coefficients
open('test_f2.fchk', 'w').write(build_fchk(electronic_structure_f2))

# get classified orbitals
orbital_type_f2 = get_custom_orbital_classification(electronic_structure_f2,
                                                    center=center_f2,
                                                    orientation=normal_f2,
                                                    orbital_numbers=frontier_orbitals)

# Print results in table
print('Inversion center?')
print('index   fragment 1   fragment 2')
for i, (o_f1, o_f2) in enumerate(zip(orbital_type_f1, orbital_type_f2)):
    print(' {:4}  {:4}  {:4.3f}  {:4}  {:4.3f}'.format(frontier_orbitals[i],
                                                       o_f1[0], o_f1[1], o_f2[0], o_f2[1]))


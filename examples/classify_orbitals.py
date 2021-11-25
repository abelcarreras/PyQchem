# Classify molecular orbitals between sigma and pi
from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input
from pyqchem.structure import Structure
from pyqchem.symmetry import get_orbital_classification
from pyqchem.file_io import build_fchk


# define coordinates
benzene_coordinates = [[ 0.00000,  1.40272,  0.00000],
                       [ 0.00000,  2.49029,  0.00000],
                       [-1.21479,  0.70136,  0.00000],
                       [-2.15666,  1.24515,  0.00000],
                       [-1.21479, -0.70136,  0.00000],
                       [-2.15666, -1.24515,  0.00000],
                       [ 0.00000, -1.40272,  0.00000],
                       [ 0.00000, -2.49029,  0.00000],
                       [ 1.21479, -0.70136,  0.00000],
                       [ 2.15666, -1.24515,  0.00000],
                       [ 1.21479,  0.70136,  0.00000],
                       [ 2.15666,  1.24515,  0.00000]]

symbols = ['C', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'C', 'H']

# create molecule
molecule = Structure(coordinates=benzene_coordinates,
                     symbols=symbols,
                     charge=0,
                     multiplicity=1)

# define Q-Chem parameters
parameters = {'jobtype': 'sp',
              'exchange': 'hf',
              'basis': '6-31G'}

# create Q-Chem input
qc_input = create_qchem_input(molecule, **parameters)

# get data from Q-Chem calculation
output, parsed_fchk = get_output_from_qchem(qc_input,
                                            processors=4,
                                            force_recalculation=False,
                                            read_fchk=True,
                                            fchk_only=True)

# write .fchk file copy on disk (for checking purpose)
txt_fchk = build_fchk(parsed_fchk)
open('test_benzene.fchk', 'w').write(txt_fchk)

# get orbital type
orbital_types = get_orbital_classification(parsed_fchk,
                                           center=[0.0, 0.0, 0.0],
                                           orientation=[0.0, 0.0, 1.0])

# print results
print('  {:5} {:5} {:5}'.format('num', 'type', 'overlap'))
for i, ot in enumerate(orbital_types):
    print('{:5}:  {:5} {:5.3f}'.format(i + 1, ot[0], ot[1]))
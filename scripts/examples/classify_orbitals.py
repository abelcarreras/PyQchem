from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input
from pyqchem.structure import Structure
from pyqchem.parsers.parser_fchk import parser_fchk
from pyqchem.symmetry import get_wf_symmetry, get_orbital_type
from pyqchem.file_io import build_fchk

#define coordinates
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

symbols = ['C', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'C', 'F']

# create molecule
molecule = Structure(coordinates=benzene_coordinates,
                     atomic_elements=symbols,
                     charge=0,
                     multiplicity=1)

# define qchem parameters
parameters = {'jobtype': 'sp',
              'exchange': 'hf',
              'basis': '6-31G',
              'gui': 2}  # necessary to parse fchk!!

# create qchem input
txt_input = create_qchem_input(molecule, **parameters)

# get data from qchem calculation
parsed_data = get_output_from_qchem(txt_input, processors=4, force_recalculation=True,
                                    parser=parser_fchk, read_fchk=True)

# write .fchk file on disk (for checking purpose)
txt_fchk = build_fchk(parsed_data)
open('test_benzene.fchk', 'w').write(txt_fchk)

# get wavefuntion symmetry data
sym_data = get_wf_symmetry(parsed_data['structure'],
                           parsed_data['basis'],
                           parsed_data['coefficients'],
                           center=[0.0, 0.0, 0.0])
# get orbital type
orbital_type = get_orbital_type(sym_data)

# print results
print('  {:5} {:5} {:5}'.format('num', 'type', 'overlap'))
for i, ot in enumerate(orbital_type):
    print('{:5}:  {:5} {:5.2f}'.format(i + 1, ot[0], ot[1]))
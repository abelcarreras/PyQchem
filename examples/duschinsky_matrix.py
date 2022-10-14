# Example of the calculation of the Duschinsky matrix for methyl peroxy radical
from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.parsers.parser_frequencies import basic_frequencies
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.tools import get_geometry_from_pubchem
from pyqchem.tools.duschinsky import get_duschinsky
from pyqchem import Structure
import numpy as np

molecule = Structure(coordinates=[[ 1.004123,  -0.180454,   0.000000],
                                  [-0.246002,   0.596152,   0.000000],
                                  [-1.312366,  -0.230256,   0.000000],
                                  [ 1.810765,   0.567203,   0.000000],
                                  [ 1.036648,  -0.805445,  -0.904798],
                                  [ 1.036648,  -0.805445,   0.904798]],
                     symbols=['C', 'O', 'O', 'H', 'H', 'H'],
                     multiplicity=2)

basis_set = '6-31+G*'
n_state = 1  # excited state number

# Optimization of ground state geometry
qc_input = QchemInput(molecule,
                      jobtype='opt',
                      exchange='b3lyp',
                      basis=basis_set,
                      )

output = get_output_from_qchem(qc_input,
                               processors=4,
                               force_recalculation=False,
                               parser=basic_optimization,
                               )

# Ground state frequencies calculation
qc_input = QchemInput(output['optimized_molecule'],
                      jobtype='freq',
                      exchange='b3lyp',
                      basis=basis_set,
                      )

gs_output = get_output_from_qchem(qc_input,
                                  processors=4,
                                  force_recalculation=False,
                                  parser=basic_frequencies,
                                  store_full_output=True,
                                  )

# Optimization of the first excited state state geometry
qc_input = QchemInput(molecule,
                      jobtype='opt',
                      exchange='b3lyp',
                      basis=basis_set,
                      cis_n_roots=10,
                      cis_singlets=True,
                      cis_triplets=False,
                      #cis_convergence=6,
                      cis_state_deriv=n_state,
                      extra_rem_keywords={'XC_GRID': '000075000302'}
                      )

output = get_output_from_qchem(qc_input,
                               processors=4,
                               force_recalculation=False,
                               parser=basic_optimization,
                               )


# Excited state frequencies calculation
qc_input = QchemInput(output['optimized_molecule'],
                      jobtype='freq',
                      exchange='b3lyp',
                      basis=basis_set,
                      cis_n_roots=10,
                      cis_singlets=True,
                      cis_triplets=False,
                      #cis_convergence=6,
                      cis_state_deriv=n_state,
                      mem_static=200,
                      extra_rem_keywords={'XC_GRID': '000075000302'}
                      )

es_output = get_output_from_qchem(qc_input,
                                  processors=4,
                                  force_recalculation=False,
                                  parser=basic_frequencies,
                                  store_full_output=True,
                                  )

# Calculation of Duschinki rotation matrix according to the model where:
#
# q_f = S * q_i + d
#
# where q_i and q_f are the normal mode coordinates of the initial and final electronic states
# S the rotation matrix and d the vector of displacements

print('frequency         displacements')
for mode in gs_output['modes']:
    print('{:10}'.format(mode['frequency']), mode['displacement'])

# build duschinsky object
duschinsky = get_duschinsky(gs_output, es_output)

# align structures along principal axis of inertia
duschinsky.align_coordinates_pmi()

# compute and print data
with np.printoptions(precision=3, suppress=True, linewidth=100):
    s = duschinsky.get_s_matrix()
    print('Rotation matrix')
    print(s)

    d = duschinsky.get_d_vector()
    print('vector of displacements')
    print(d[None].T)

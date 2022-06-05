# Example of the calculation of the Duschinsky matrix for ethene
from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.parsers.parser_frequencies import basic_frequencies
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.tools import get_geometry_from_pubchem
from pyqchem.tools.duschinsky import get_duschinsky

import numpy as np


# define molecule
ethene = get_geometry_from_pubchem('ethene')

basis_set = 'cc-pVDZ'
n_state = 1  # excited state number

# Optimization of ground state geometry
qc_input = QchemInput(ethene,
                      jobtype='opt',
                      exchange='b3lyp',
                      basis=basis_set,
                      )

output = get_output_from_qchem(qc_input,
                               processors=4,
                               force_recalculation=False,
                               parser=basic_optimization,
                               )

gs_structure = output['optimized_molecule']

# Optimization of the first excited state state geometry
qc_input = QchemInput(ethene,
                      jobtype='opt',
                      exchange='b3lyp',
                      basis=basis_set,
                      cis_n_roots=8,
                      cis_singlets=True,
                      cis_triplets=False,
                      cis_convergence=8,
                      cis_state_deriv=n_state,
                      )

output = get_output_from_qchem(qc_input,
                               processors=4,
                               force_recalculation=False,
                               parser=basic_optimization,
                               )

es_structure = output['optimized_molecule']


# Ground state frequencies calculation
qc_input = QchemInput(gs_structure,
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

# Excited state frequencies calculation
qc_input = QchemInput(es_structure,
                      jobtype='freq',
                      exchange='b3lyp',
                      basis=basis_set,
                      cis_n_roots=8,
                      cis_singlets=True,
                      cis_triplets=False,
                      cis_convergence=8,
                      cis_state_deriv=n_state,
                      mem_static=200,
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

duschinsky = get_duschinsky(gs_output, es_output)

# align structures along principal axis of inertia
duschinsky.align_coordinates()

with np.printoptions(precision=3, suppress=True, linewidth=100):
    s = duschinsky.get_s_matrix()
    print('Rotation matrix')
    print(s)

    d = duschinsky.get_d_vector()
    print('vector of displacements')
    print(d[None].T)

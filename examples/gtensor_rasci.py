# example of the calculation of the g-tensor
from pyqchem.tools.gtensor import GTensor
from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.parsers.parser_rasci import parser_rasci
from pyqchem.structure import Structure
import numpy as np

# generate molecule
molecule = Structure(coordinates=[[ 0.0000000000,  0.0000000000,  0.1155512138],
                                  [-0.8125311840, -0.0000000000, -0.4622048550],
                                  [ 0.8125311840, -0.0000000000, -0.4622048550]],
                     symbols=['O', 'H', 'H'],
                     charge=1,
                     multiplicity=2)

# create qchem input
qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      unrestricted=False,
                      correlation='rasci',
                      basis='sto-3g',
                      scf_algorithm='gdm',
                      ras_act=7,
                      ras_elec_alpha=6,
                      ras_elec_beta=4,
                      set_iter=60,
                      max_scf_cycles=60,
                      ras_roots=5,
                      ras_occ=0,
                      calc_soc=True,
                      extra_rem_keywords={'RAS_SOC_SYM_DENS': True})

# calculate and parse qchem output
parsed_data = get_output_from_qchem(qc_input,
                                    processors=4,
                                    parser=parser_rasci,
                                    )


gtensor = GTensor(parsed_data,
                  # selected_states=[1, 2, 3],  # States to be included
                  )

g_tensor = gtensor.get_g_tensor()
g_shifts = gtensor.get_g_shift(use_ppm=False)
total_states = gtensor.get_number_of_states()
selected_states = gtensor.get_selected_states()


print('---------------------')
print('      G-TENSOR       ')
print('---------------------')
print('- Number of states in output: ', total_states)
print('- Selected states: ', selected_states)
print('g-tensor: ', np.real(g_tensor))
print('g-factor (x, y, z) [ppt]: ', np.real(g_shifts))
print()

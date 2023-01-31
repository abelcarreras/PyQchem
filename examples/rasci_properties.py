# Compute and extract plot properties of RASCI wave function
from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.parsers.parser_rasci import parser_rasci
from pyqchem.structure import Structure
import numpy as np

# generate molecule
molecule = Structure(coordinates=[[0.0, 0.0, 0.0],
                                  [0.0, 0.0, 0.9]],
                     symbols=['H', 'H'],
                     charge=0,
                     multiplicity=1)

# create qchem input
qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      correlation='rasci',
                      basis='sto-3g',
                      ras_act=2,
                      ras_elec=2,
                      ras_spin_mult=0,
                      ras_roots=2,
                      ras_occ=0,
                      state_analysis=True, # request WF analysis
                      ras_do_part=False,
                      ras_do_hole=False)

# calculate and parse qchem output
data, ee = get_output_from_qchem(qc_input,
                                 processors=4,
                                 parser=parser_rasci,
                                 return_electronic_structure=True
                                 )

# print properties
print('RASCI WF properties\n---------------')
print('Natural Orbitals Occupations: ', ee['nato_occupancies'])
print('Natural Orbitals Coefficients: ', ee['nato_coefficients'])
print('')

print('Natural Transition Orbitals')
for nto in ee['ntos_ras'].items():
    print('Transition ', nto[0])
    print('NTO occupancies: ', nto[1]['occupancies'])
    print('Origin orbital')
    print(np.array(nto[1]['V']))
    print('Target orbital')
    print(np.array(nto[1]['U']))

print('')
print('Fraction occupation density:\n', np.array(ee['fractional_occupation_density']))

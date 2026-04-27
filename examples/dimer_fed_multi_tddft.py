# Compute the excited states using TD-DFT
HARTREE_TO_EV = 27.211386245988
BOHR_TO_ANG = 0.529177249

from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.parsers.parser_cis import basic_cis
from pyqchem.structure import Structure
from pyqchem.tools import get_geometry_from_pubchem
from pyqchem.qc_input import CustomSection
import numpy as np


mol_naph = get_geometry_from_pubchem('naphthalene')

qc_input = QchemInput(mol_naph,
                      jobtype='sp',
                      exchange='b3lyp',
                      unrestricted=False,
                      basis='sto-3g',
                      cis_n_roots=6,
                      cis_convergence=8,
                      cis_singlets=True,
                      cis_triplets=False,
                      max_cis_cycles=100,
                      )

# calculate and parse qchem output
data = get_output_from_qchem(qc_input,
                             processors=4,
                             parser=basic_cis,
                             store_full_output=True,
                             )

print('\nMONOMER STATES')
for state in data['excited_states']:
    print(state['excitation_energy'], state['transition_moment'])

state_num = 2
print(state_num)
print('selected: ', data['excited_states'][state_num-1])

tdm = data['excited_states'][state_num-1]['transition_moment']


# define distances
distances = np.linspace(4, 14, 20)

couplings = []
couplings_dipole = []

print('\nDIMER STATES')
dist = 2

# generate molecule
monomer_1_c = mol_naph.get_coordinates()
monomer_2_c = [[c[0], c[1], c[2] + dist] for c in monomer_1_c]
monomer_s = mol_naph.get_symbols()

dimer = Structure(coordinates=monomer_1_c + monomer_2_c,
                  symbols=monomer_s + monomer_s,
                  charge=0,
                  multiplicity=1)

# dipole in X
# create qchem input
qc_input = QchemInput(dimer,
                      jobtype='sp',
                      exchange='b3lyp',
                      unrestricted=False,
                      basis='sto-3g',
                      cis_n_roots=6,
                      cis_convergence=14,
                      cis_singlets=True,
                      cis_triplets=False,
                      max_cis_cycles=100,
                      extra_rem_keywords={'STS_FED': True,
                                          'STS_MULTI_NROOTS': 4,
                                          'STS_MULTI_PRINT': 1,
                                          'STS_DONOR': '1-18',
                                          'STS_ACCEPTOR': '19-36'
                                          },
                      extra_sections=[CustomSection(title='localized_diabatization',
                                                    keywords={'adiabatic states\n': ' '.join(['1', '2', '3', '4'])}
                                                    )
                                      ],

                      )

# calculate and parse qchem output
data = get_output_from_qchem(qc_input,
                             processors=4,
                             # parser=basic_cis,
                             store_full_output=True,
                             )

print(data)

# Compute the excited states using TD-DFT
HARTREE_TO_EV = 27.211386245988
BOHR_TO_ANG = 0.529177249

from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.parsers.parser_cis import basic_cis
from pyqchem.structure import Structure
from pyqchem.tools import get_geometry_from_pubchem
import matplotlib.pyplot as plt
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
for dist in distances:

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
                                              'STS_DONOR': '1-18',
                                              'STS_ACCEPTOR': '19-36'}
                          )

    # calculate and parse qchem output
    data = get_output_from_qchem(qc_input,
                                 processors=4,
                                 parser=basic_cis,
                                 store_full_output=True,
                                 )

    state_pair = 3
    i, j = 2*(state_pair-1)+1, 2*(state_pair-1)+2
    print(i, j)
    couplings.append(abs(data['FED'][(i, j)]['coupling']))
    print('selected: ', data['excited_states'][i-1])
    print('selected: ', data['excited_states'][j-1])

    for state in data['excited_states']:
        print(state['excitation_energy'], state['transition_moment'])

    tdm = np.array(tdm, dtype=float)
    r_vec = np.array([0, 0, dist / BOHR_TO_ANG], dtype=float)

    r = np.linalg.norm(r_vec)
    r_hat = r_vec / r

    V_hartree = (np.dot(tdm, tdm) - 3 * np.dot(tdm, r_hat) * np.dot(tdm, r_hat)) / r**3
    couplings_dipole.append(V_hartree * HARTREE_TO_EV)

    print('-------------')

print(couplings)
plt.title('Couplings')
plt.xlabel('Distance [A]')
plt.ylabel('FED Coupling [Hartree]')
plt.yscale('log')
plt.plot(distances, couplings, label='FED')
plt.plot(distances, couplings_dipole, label='Dipole')
plt.legend()
plt.show()
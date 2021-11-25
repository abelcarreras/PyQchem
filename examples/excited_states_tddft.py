# Compute the excited states using TD-DFT
from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.parsers.parser_cis import basic_cis
from pyqchem.structure import Structure

import numpy as np
import matplotlib.pyplot as plt

# define distances
distances = np.arange(0.1, 3.5, 0.1)

energies = []
for dist in distances:

    # generate molecule
    molecule = Structure(coordinates=[[0.0, 0.0, 0.0],
                                      [0.0, 0.0, dist]],
                         symbols=['H', 'H'],
                         charge=0,
                         multiplicity=1)

    # create qchem input
    qc_input = QchemInput(molecule,
                          jobtype='sp',
                          exchange='b3lyp',
                          unrestricted=False,
                          basis='sto-3g',
                          cis_n_roots=4,
                          cis_convergence=14,
                          cis_singlets=True,
                          cis_triplets=True)

    # calculate and parse qchem output
    data = get_output_from_qchem(qc_input,
                                 processors=4,
                                 parser=basic_cis,
                                 store_full_output=True,
                                 )

    # store energies of excited states states in list
    energies.append([data['scf_energy']] + [state['total_energy'] for state in data['excited_states']])

multiplicity_list = ['Singlet'] + [state['multiplicity'] for state in data['excited_states']]

# transform energy list in array
energies = np.array(energies)

# print energies
print('\ndistance  state 1    state 2')
print('----------------------------')
for d, e in zip(distances, energies):
    print('{:4.1f} : {:10.5f} {:10.5f}'.format(d, e[0], e[1]))

# plot energies
for i, state_energy in enumerate(energies.T):
    plt.plot(distances, state_energy, label='state {} ({})'.format(i+1, multiplicity_list[i]))

plt.title('Dissociation of Hydrogen molecule')
plt.xlim([0, 3.5])
plt.ylim([-1.5, 1.0])
plt.xlabel('Distance [A]')
plt.ylabel('Energy [Hartree]')
plt.axhline(-1.0, color='grey', linestyle='--')

plt.legend()
plt.show()
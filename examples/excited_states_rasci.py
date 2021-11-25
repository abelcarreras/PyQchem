# Compute the excited states using RASCI method
from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.parsers.parser_rasci import parser_rasci
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
                          exchange='hf',
                          correlation='rasci',
                          basis='sto-3g',
                          ras_act=2,
                          ras_elec=2,
                          ras_spin_mult=0,
                          ras_roots=2,
                          ras_do_part=False,
                          ras_do_hole=False)

    # calculate and parse qchem output
    data = get_output_from_qchem(qc_input,
                                 processors=4,
                                 parser=parser_rasci,
                                 store_full_output=True,
                                 )

    # store energies of excited states states in list
    energies.append([state['total_energy'] for state in data['excited_states']])

multiplicity_list = [state['multiplicity'] for state in data['excited_states']]

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
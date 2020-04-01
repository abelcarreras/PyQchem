from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.parsers.basic import basic_parser_qchem
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
                          basis='6-31G(d,p)',
                          ras_act=2,
                          ras_elec=2,
                          ras_spin_mult=0,
                          ras_roots=2,
                          ras_do_hole=True,
                          ras_sts_tm=True,
                          # rasci sr-dft
                          ras_omega=400,
                          ras_srdft_cor='srpw92',
                          ras_srdft_exc='srpbe',
                          ras_natorb=False,
                          ras_print=0,
                          set_iter=30,
                          ras_srdft_damp=0.5)

    # calculate and parse qchem output

    data = get_output_from_qchem(qc_input,
                                 processors=4,
                                 parser=basic_parser_qchem)

    # store energies of excited states states in list
    energies.append([state['total energy'] for state in data['excited states']])

multiplicity_list = [state['multiplicity'] for state in data['excited states']]

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
plt.xlabel('Distance [Å]')
plt.ylabel('Energy [Hartree]')
plt.axhline(-1.0, color='grey', linestyle='--')

plt.legend()
plt.show()
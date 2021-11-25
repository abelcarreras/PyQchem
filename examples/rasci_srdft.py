# Example of the use of srDFT method in RASCI
from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.parsers.parser_rasci import parser_rasci
from pyqchem.structure import Structure
from pyqchem.tools import rotate_coordinates

import numpy as np
import matplotlib.pyplot as plt

# define distances
distances = np.arange(0.0, np.pi+np.pi/21, np.pi/21)

energies = []
for dist in distances:

    # generate molecule
    coordinates = [[ 0.6695,  0.0000,  0.0000],
                   [-0.6695,  0.0000,  0.0000],
                   [ 1.2321,  0.9289,  0.0000],
                   [ 1.2321, -0.9289,  0.0000],
                   [-1.2321,  0.9289,  0.0000],
                   [-1.2321, -0.9289,  0.0000]]

    coordinates = rotate_coordinates(coordinates, angle=dist, axis=[1, 0, 0], atoms_list=[4, 5])

    ethene = Structure(coordinates=coordinates,
                       symbols=['C', 'C', 'H', 'H', 'H', 'H'],
                       charge=0,
                       multiplicity=1)

    # create qchem input
    qc_input = QchemInput(ethene,
                          jobtype='sp',
                          exchange='hf',
                          correlation='rasci',
                          basis='6-31G',
                          ras_act=4,
                          ras_elec=2,
                          ras_spin_mult=0,
                          ras_roots=2,
                          ras_do_part=False,
                          ras_do_hole=False,
                          ras_srdft_spinpol=True,
                          ras_omega=400,
                          max_scf_cycles=100,
                          ras_srdft_cor='srPBE',
                          #ras_srdft_exc='srPBE'
                          )

    # calculate and parse qchem output
    data, ee = get_output_from_qchem(qc_input,
                                     processors=3,
                                     parser=parser_rasci,
                                     read_fchk=True,
                                     # store_full_output=True,
                                     )

    guess = ee['coefficients']

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
    plt.plot(distances, state_energy, label='state {} (S={})'.format(i+1, multiplicity_list[i]))

plt.title('Torsion of ethylene')
plt.xlabel('Torsion angle [Rad]')
plt.ylabel('Energy [Hartree]')

plt.legend()
plt.show()
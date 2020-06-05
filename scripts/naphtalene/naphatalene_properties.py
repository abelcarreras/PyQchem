from pyqchem import get_output_from_qchem, Structure, QchemInput
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.parsers.parser_rasci import parser_rasci
from pyqchem.file_io import build_fchk
from pyqchem.symmetry import get_orbital_classification
from pyqchem.qchem_core import redefine_calculation_data_filename
from pyqchem.utils import get_inertia
from pyqchem.units import DEBYE_TO_AU, EV_TO_NM
from pyqchem.tools import print_excited_states
from pyqchem.symmetry import get_state_symmetry

import numpy as np
import matplotlib.pyplot as plt


def draw_molecule(molecule):
    import mogli

    mol = mogli.Molecule(positions=molecule.get_coordinates(),
                         atomic_numbers=molecule.get_atomic_numbers())

    mogli.show(mol, bonds_param=1.5, title='Molecule')

redefine_calculation_data_filename('naphtalene.pkl')


# define molecule
coordinates = [[ 2.4610326539,  0.7054950347, -0.0070507104],
                 [ 1.2697800226,  1.4213478618,  0.0045894884],
                 [ 0.0071248839,  0.7134976955,  0.0071917580],
                 [-1.2465927908,  1.4207541565,  0.0039025332],
                 [ 2.4498274919, -0.7358510124,  0.0046346543],
                 [ 3.2528295760,  1.2280710625, -0.0312673955],
                 [ 1.3575083440,  2.3667492466,  0.0220260183],
                 [-1.2932627225,  2.3688000888, -0.0152164523],
                 [ 3.2670227933, -1.2176289251,  0.0251089819],
                 [-2.4610326539, -0.7054950347,  0.0070507104],
                 [-1.2697800226, -1.4213478618, -0.0045894884],
                 [-0.0071248839, -0.7134976955, -0.0071917580],
                 [ 1.2465927908, -1.4207541565, -0.0039025332],
                 [-2.4498274919,  0.7358510124, -0.0046346543],
                 [-3.2528295760, -1.2280710625,  0.0312673955],
                 [-1.3575083440, -2.3667492466, -0.0220260183],
                 [ 1.2932627225, -2.3688000888,  0.0152164523],
                 [-3.2670227933,  1.2176289251, -0.0251089819]]

symbols = ['C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H',
           'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H']

molecule = Structure(coordinates=coordinates,
                     symbols=symbols,
                     charge=0,
                     multiplicity=1)

#print(molecule.get_xyz())
# draw_molecule(molecule)

print('Initial structure')
print(molecule)


qc_input_hf = QchemInput(molecule,
                         unrestricted=True,
                         jobtype='sp',
                         exchange='hf',
                         basis='sto-3g'
                         )

#print(qc_input.get_txt())
parsed_data, electronic_structure = get_output_from_qchem(qc_input_hf,
                                                          processors=4,
                                                          force_recalculation=False,
                                                          # parser=rasci,
                                                          read_fchk=True,
                                                          store_full_output=True,
                                                          )

# create qchem input
qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      unrestricted=False,
                      correlation='rasci',
                      basis='sto-3g',
                      scf_guess=electronic_structure['coefficients'],
                      ras_act=2,
                      ras_elec=2,
                      # ras_occ=31,
                      ras_spin_mult=0,
                      ras_roots=15,
                      ras_do_hole=True,
                      ras_do_part=True,
                      # ras_srdft_cor='srpbe',
                      # ras_srdft_exc='srpbe',
                      # ras_omega=300,
                      # skip_scfman=True,
                      # mom_start=True,
                      set_iter=60)


#print(qc_input.get_txt())
parsed_data, electronic_structure = get_output_from_qchem(qc_input,
                                                          processors=4,
                                                          force_recalculation=False,
                                                          parser=parser_rasci,
                                                          read_fchk=True,
                                                          store_full_output=True,
                                                          )

# open('reference.fchk', 'w').write(build_fchk(electronic_structure))

# get orbital type
orbital_types = get_orbital_classification(electronic_structure,
                                           center=[0.0, 0.0, 0.0],
                                           orientation=[0.0, 0.0, 1.0])

# print results
print('  {:5} {:5} {:5}'.format('num', 'type', 'overlap'))
for i, ot in enumerate(orbital_types):
    print('{:5}:  {:5} {:5.3f}'.format(i + 1, ot[0], ot[1]))

# print(molecule.alpha_electrons, molecule.beta_electrons)

print('\nAdiabatic states\n--------------------')
print_excited_states(parsed_data['excited_states'], include_conf_rasci=True)

inertia_moments, inertia_axis = get_inertia(electronic_structure['structure'])

#print(im)
#print(ivec[0])
#print(ivec[1])
#print(ivec[2])

print('\nInertia\n----------')
print(inertia_moments)
print(np.array(inertia_axis))

print(electronic_structure['structure'])
symmetry_measures = get_state_symmetry(electronic_structure,
                                       parsed_data['excited_states'],
                                       center=[0, 0, 0],
                                       orientation=inertia_axis[0],
                                       orientation2=inertia_axis[1],
                                       group='D2h',
                                       extra_print=False
                                       )

energies = [state['excitation_energy'] for state in parsed_data['excited_states']]

print('\nSymmetry of RASCI excited states\n--------------------------------')
for energy, state in zip(energies, symmetry_measures.items()):
    print('{:8} '.format(state[0]) + ' {:5} {:5.3f}'.format(*state[1]) + '  {:5.3f} eV'.format(energy))


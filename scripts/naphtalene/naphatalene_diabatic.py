from pyqchem.parsers.parser_rasci import parser_rasci as rasci_parser
from pyqchem import get_output_from_qchem, Structure, QchemInput
from pyqchem.file_io import build_fchk
from pyqchem.symmetry import get_symmetry_le
from pyqchem.qchem_core import redefine_calculation_data_filename
from pyqchem.utils import get_plane
from pyqchem.utils import get_ratio_of_condition
from pyqchem.tools import print_excited_states, plot_rasci_state_configurations
from pyqchem.plots import plot_diabatization
import numpy as np


redefine_calculation_data_filename('diabatic_naphtalene.pkl')


# define dimer
coor_monomer1 = [[ 2.4610326539,  0.7054950347, -0.0070507104],
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

# set dimer geometry
coor_monomer2 = np.array(coor_monomer1).copy()
coor_monomer2[:, 2] += 4.0
# coor_monomer2[:, 1] -= 5.0
coor_monomer2 = list(coor_monomer2)

symbols_monomer = ['C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H',
                   'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H']

coordinates = coor_monomer1 + coor_monomer2
symbols_dimer = symbols_monomer * 2

dimer = Structure(coordinates=coordinates,
                  symbols=symbols_dimer,
                  charge=0,
                  multiplicity=1)

range_f1 = range(len(coor_monomer1))
range_f2 = range(len(coor_monomer1), len(coordinates))

print('Dimer structure')
print(dimer)

# RASCI qchem input
qc_input = QchemInput(dimer,
                      unrestricted=False,
                      jobtype='sp',
                      exchange='hf',
                      correlation='rasci',
                      basis='sto-3g',
                      ras_act=4,
                      ras_elec=4,
                      # ras_occ=66,
                      ras_spin_mult=1,  # singlets only
                      ras_roots=10,      # calculate 8 states
                      ras_do_hole=True,
                      ras_do_part=True,
                      # ras_srdft_cor='srpbe',
                      # ras_srdft_exc='srpbe',
                      # ras_omega=300,
                      set_iter=60)


parsed_data, electronic_structure = get_output_from_qchem(qc_input,
                                                          processors=4,
                                                          force_recalculation=False,
                                                          parser=rasci_parser,
                                                          read_fchk=True,
                                                          store_full_output=True
                                                          )

print('final structure')
print(electronic_structure['structure'])

# Analysis of diabatic states to use in diabatization
print('\nAdiabatic states to use in diabatization (1e, max_jump 4)')
list_diabatic = []
for i, state in enumerate(parsed_data['excited_states']):
    ratio = get_ratio_of_condition(state, n_electron=1, max_jump=4)

    if ratio > 0.1:
        mark = 'X'
        list_diabatic.append(i+1)
    else:
        mark = ''

    print('State {}: {:4.3f}  {}'.format(i+1, ratio, mark))


# sequential diabatization scheme (2 steps)
diabatization_scheme = [
                        # {'method': 'ER', 'states': list(range(1, len(list_diabatic)+1))},
                        {'method': 'Boys', 'states': list(range(1, len(list_diabatic)+1))},
                        # {'method': 'ER', 'states': [3, 4, 5, 6]},
                        # {'method': 'DQ', 'states': list(range(1, len(list_diabatic)+1)), 'paramters': [0.9]}
                        ]

# RASCI qchem input
qc_input = QchemInput(electronic_structure['structure'],
                      jobtype='sp',
                      exchange='hf',
                      correlation='rasci',
                      basis='sto-3g',
                      scf_guess=electronic_structure['coefficients'],
                      ras_act=4,
                      ras_elec=4,
                      # ras_occ=66,
                      ras_spin_mult=1,  # singlets only
                      ras_roots=10,      # calculate 8 states
                      ras_do_hole=True,
                      ras_do_part=True,
                      # ras_srdft_cor='srpbe',
                      # ras_srdft_exc='srpbe',
                      # ras_omega=300,
                      ras_diabatization_states=list_diabatic,  # adiabatic states selected for diabatization
                      ras_diabatization_scheme=diabatization_scheme,
                      set_iter=60)

# print(qc_input.get_txt())

parsed_data, electronic_structure = get_output_from_qchem(qc_input,
                                                          processors=14,
                                                          force_recalculation=False,
                                                          parser=rasci_parser,
                                                          read_fchk=True,
                                                          store_full_output=True
                                                          )

#print(parsed_data)
#exit()
np.set_printoptions(precision=3)

print('\nAdiabatic states\n--------------------')
print_excited_states(parsed_data['excited_states'])

#plot_rasci_state_configurations(parsed_data['excited_states'])

# diabatization analysis
diabatization = parsed_data['diabatization']

print('\nRotation Matrix')
print(np.array(diabatization['rot_matrix']))

print('\nAdiabatic Matrix')
print(np.array(diabatization['adiabatic_matrix']))

print('\nDiabatic Matrix')
print(np.array(diabatization['diabatic_matrix']))

print('\nDiabatic states\n--------------------')
print_excited_states(diabatization['diabatic_states'], include_mulliken_rasci=True)


plot_diabatization(diabatization['diabatic_states'], atoms_ranges=[dimer.get_number_of_atoms()/2,
                                                                   dimer.get_number_of_atoms()])

coordinates_f1 = np.array(electronic_structure['structure'].get_coordinates())[range_f1]
center_f1, normal_f1 = get_plane(coordinates_f1)

coordinates_f2 = np.array(electronic_structure['structure'].get_coordinates())[range_f2]
center_f2, normal_f2 = get_plane(coordinates_f2)

print('\nDiabatic states symmetry analysis (monomer 1)')
sym_data = get_symmetry_le(electronic_structure, parsed_data, fragment_atoms=range_f2, group='D2h')
print('Symmetry: ', sym_data)

print('\nDiabatic states symmetry analysis (monomer 2)')
sym_data = get_symmetry_le(electronic_structure, parsed_data, fragment_atoms=range_f1, group='D2h')
print('Symmetry: ', sym_data)

from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.parsers.parser_rasci import parser_rasci
from pyqchem.structure import Structure
from pyqchem.plots import plot_diabatization, plot_state
from pyqchem.utils import get_ratio_of_condition

import numpy as np
import matplotlib.pyplot as plt


# define monomer
coor_monomer = [[ 0.6695,  0.0000,  0.0000],
                [-0.6695,  0.0000,  0.0000],
                [ 1.2321,  0.9289,  0.0000],
                [ 1.2321, -0.9289,  0.0000],
                [-1.2321,  0.9289,  0.0000],
                [-1.2321, -0.9289,  0.0000]]

symbols_monomer = ['C', 'C', 'H', 'H', 'H', 'H']

monomer = Structure(coordinates=coor_monomer,
                    symbols=symbols_monomer,
                    charge=0,
                    multiplicity=1)

# optimization qchem input
qc_input = QchemInput(monomer,
                      jobtype='opt',
                      exchange='hf',
                      basis='sto-3g',
                      geom_opt_tol_gradient=300,
                      geom_opt_tol_energy=100,
                      geom_opt_coords=-1,
                      geom_opt_tol_displacement=1200)

parsed_data = get_output_from_qchem(qc_input,
                                    processors=4,
                                    parser=basic_optimization)

opt_monomer = parsed_data['optimized_molecule']

print('Optimized monomer structure')
print(opt_monomer)

# Build dimer from monomer
coor_monomer2 = np.array(opt_monomer.get_coordinates())
coor_monomer2[:, 2] += 4.0  # monomers separation

coordinates = opt_monomer.get_coordinates() + coor_monomer2.tolist()
symbols_dimer = symbols_monomer * 2

dimer = Structure(coordinates=coordinates,
                  symbols=symbols_dimer,
                  charge=0,
                  multiplicity=1)

print('Dimer structure')
print(dimer)

# sequential diabatization scheme (2 steps)
diabatization_scheme = [{'method': 'ER', # first step
                         'states': [1, 2, 3, 4]},  # in order with respect to selected states
                        {'method': 'DQ', # Note: This second step is not really needed in this case, just for demostration
                         'states': [3, 4],  # in energy order (low->high) with respect to diabatic states in previous step
                         'parameters': 0.0}]

# RASCI qchem input
qc_input = QchemInput(dimer,
                      jobtype='sp',
                      exchange='hf',
                      correlation='rasci',
                      basis='sto-3g',
                      ras_act=6,
                      ras_elec=4,
                      ras_spin_mult=1,  # singlets only
                      ras_roots=8,      # calculate 8 states
                      ras_do_hole=False,
                      ras_do_part=False,
                      ras_diabatization_states=[3, 4, 5, 6],  # adiabatic states selected for diabatization
                      ras_diabatization_scheme=diabatization_scheme,
                      set_iter=30)

# print(qc_input.get_txt())

parsed_data = get_output_from_qchem(qc_input,
                                    processors=4,
                                    force_recalculation=False,
                                    parser=parser_rasci
                                    )

# parsed_data = rasci_parser(parsed_data)
# print(parsed_data)

# print adiabatic states
print('Adiabatic states dimer\n--------------------')
for i, state in enumerate(parsed_data['excited_states']):
    print('\nState {}'.format(i+1))
    print('Transition DM: ', state['transition_moment'])
    print('Energy: ', state['excitation_energy'])
    print(' Alpha  Beta   Amplitude')
    for j, conf in enumerate(state['configurations']):
        print('  {}  {} {:8.3f}'.format(conf['alpha'], conf['beta'], conf['amplitude']))

# plot adiabatic states
for i, state in enumerate(parsed_data['excited_states']):
    plot_state(state, with_amplitude=True, orbital_range=[qc_input._ras_occ,
                                                          qc_input._ras_occ + qc_input._ras_act])
    plt.title('Adiabatic State {}'.format(i+1))

plt.show()


# Analysis of diabatic states to use in diabatization
print('\nAdiabatic states to use in diabatization (1e, max_jump 4)')
for i, state in enumerate(parsed_data['excited_states']):
    ratio = get_ratio_of_condition(state, n_electron=1, max_jump=4)
    mark = 'X' if ratio > 0.5 else ''
    print('State {}: {:4.3f}  {}'.format(i+1, ratio, mark))


# diabatization analysis
diabatization = parsed_data['diabatization']

print('\nRotation Matrix')
print(np.array(diabatization['rot_matrix']))

print('\nAdiabatic Matrix')
print(np.array(diabatization['adiabatic_matrix']))

print('\nDiabatic Matrix')
print(np.array(diabatization['diabatic_matrix']))


print('\nDiabatic states dimer\n--------------------')

for i, state in enumerate(diabatization['diabatic_states']):
    print('\nState {}'.format(i+1))
    print('Transition DM: ', state['transition_moment'])
    print('Energy: ', state['excitation_energy'])

plot_diabatization(diabatization['diabatic_states'], atoms_ranges=[dimer.get_number_of_atoms()/2,
                                                                   dimer.get_number_of_atoms()])


# Monomer adiabatic states (extra test)
qc_input = QchemInput(opt_monomer,
                      jobtype='sp',
                      exchange='hf',
                      correlation='rasci',
                      basis='sto-3g',
                      ras_act=6,
                      ras_elec=4,
                      ras_spin_mult=1,  # singlets only
                      ras_roots=4,      # calculate 8 states
                      ras_do_hole=False,
                      ras_do_part=False,
                      set_iter=30)

parsed_data = get_output_from_qchem(qc_input,
                                    processors=14,
                                    force_recalculation=False,
                                    parser=parser_rasci
                                    )

print('\nAdiabatic states monomer\n--------------------')
for i, state in enumerate(parsed_data['excited_states']):
    print('\nState {}'.format(i+1))
    print('Transition DM: ', state['transition_moment'])
    print('Energy: ', state['excitation_energy'])
    print(' Alpha  Beta   Amplitude')
    for j, conf in enumerate(state['configurations']):
        print('  {}  {} {:8.3f}'.format(conf['alpha'], conf['beta'], conf['amplitude']))

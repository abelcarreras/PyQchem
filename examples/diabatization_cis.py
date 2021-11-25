# Application of diabatization method implemented in CIS for ethylene dimer
from pyqchem import Structure, QchemInput, get_output_from_qchem
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.parsers.parser_cis import basic_cis
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
coor_monomer2[:, 2] += 6.0  # monomers separation

coordinates = opt_monomer.get_coordinates() + coor_monomer2.tolist()
symbols_dimer = symbols_monomer * 2

dimer = Structure(coordinates=coordinates,
                  symbols=symbols_dimer,
                  charge=0,
                  multiplicity=1)

print('Dimer structure')
print(dimer)


# CIS qchem input
qc_input = QchemInput(dimer,
                      exchange='hf',
                      basis='sto-3g',
                      #unrestricted=True,
                      cis_n_roots=8,
                      cis_convergence=8,
                      cis_singlets=True,
                      cis_triplets=False,
                      cis_ampl_anal=True,
                      RPA=False,
                      mem_static=900,
                      set_iter=30
                      )

print(qc_input.get_txt())

parsed_data = get_output_from_qchem(qc_input,
                                    processors=4,
                                    force_recalculation=False,
                                    parser=basic_cis
                                    )


# Analysis of diabatic states to use in diabatization
print('\nStates to use in diabatization (1e, max_jump 3)')
list_diabatic = []
for i, state in enumerate(parsed_data['excited_states']):
    ratio = get_ratio_of_condition(state, n_electron=1, max_jump=3)
    if ratio > 0.9:
        list_diabatic.append(i+1)
        mark = 'X'
    else:
        mark = ''
    print('State {}: {:4.3f}  {}'.format(i+1, ratio, mark))


# RASCI qchem input
qc_input = QchemInput(dimer,
                      exchange='hf',
                      basis='sto-3g',
                      #unrestricted=True,
                      cis_n_roots=8,
                      cis_convergence=8,
                      cis_singlets=True,
                      cis_triplets=False,
                      cis_ampl_anal=True,
                      loc_cis_ov_separate=False,
                      namd_nsurfaces=0,
                      er_cis_numstate=len(list_diabatic),
                      cis_diabath_decompose=True,
                      localized_diabatization=list_diabatic,
                      RPA=False,
                      mem_static=900,
                      set_iter=30
                      )

# print(qc_input.get_txt())

parsed_data = get_output_from_qchem(qc_input,
                                    processors=4,
                                    force_recalculation=False,
                                    parser=basic_cis,
                                    store_full_output=True,
                                    )

# parsed_data = rasci_parser(parsed_data)
# print(parsed_data)

# print adiabatic states
print('\nAdiabatic states dimer\n--------------------')
for i, state in enumerate(parsed_data['excited_states']):
    print('\nState {}'.format(i+1))
    print('Transition DM: ', state['transition_moment'])
    print('Energy: ', state['excitation_energy'])
    print(' Origin  Target   Amplitude')
    for j, conf in enumerate(state['configurations']):
        print('{:^9} {:^7} {:8.3f}'.format(conf['origin'], conf['target'], conf['amplitude']))

# plot adiabatic states
for i, state in enumerate(parsed_data['excited_states']):
    plot_state(state, with_amplitude=True, orbital_range=[dimer.alpha_electrons-4, dimer.alpha_electrons+4])
    plt.title('Adiabatic State {}'.format(i+1))
plt.show()

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
    mat = [state['transition_moment'] for i, state in enumerate(parsed_data['excited_states']) if i+1 in list_diabatic]
    tdm = np.dot(np.array(diabatization['rot_matrix']).T, mat)[i]
    print('Transition DM: ', tdm)
    print('Energy: ', state['excitation_energy'])

plot_diabatization(diabatization['diabatic_states'], atoms_ranges=[dimer.get_number_of_atoms()/2,
                                                                   dimer.get_number_of_atoms()])


# Monomer adiabatic states (extra test)
qc_input = QchemInput(dimer,
                      #method='cis',
                      exchange='hf',
                      basis='sto-3g',
                      #unrestricted=True,
                      cis_n_roots=10,
                      cis_convergence=8,
                      cis_singlets=True,
                      cis_triplets=False,
                      cis_ampl_anal=True,
                      RPA=False,
                      mem_static=900,
                      set_iter=30
                      )

parsed_data = get_output_from_qchem(qc_input,
                                    processors=14,
                                    force_recalculation=False,
                                    parser=basic_cis
                                    )

print('\nAdiabatic states monomer\n--------------------')
for i, state in enumerate(parsed_data['excited_states']):
    print('\nState {}'.format(i+1))
    print('Transition DM: ', state['transition_moment'])
    print('Energy: ', state['excitation_energy'])
    print(' Origin  Target   Amplitude')
    for j, conf in enumerate(state['configurations']):
        print('{:^9} {:^7} {:8.3f}'.format(conf['origin'], conf['target'], conf['amplitude']))

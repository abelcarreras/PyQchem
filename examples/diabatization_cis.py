# Application of diabatization method implemented in CIS for ethylene dimer
from pyqchem import Structure, QchemInput, get_output_from_qchem
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.parsers.parser_cis import basic_cis
from pyqchem.plots import plot_diabatization, plot_state
from pyqchem.utils import get_ratio_of_condition
from pyqchem.symmetry import get_state_symmetry

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
                      basis='6-311g',
                      cis_n_roots=8,
                      cis_convergence=8,
                      cis_singlets=True,
                      cis_triplets=False,
                      cis_ampl_anal=True,
                      RPA=False,
                      mem_static=900,
                      set_iter=30
                      )

# print(qc_input.get_txt())

parsed_data, ee = get_output_from_qchem(qc_input,
                                        processors=4,
                                        force_recalculation=False,
                                        parser=basic_cis,
                                        return_electronic_structure=True
                                        )


symmetry_measures = get_state_symmetry(ee,
                                       parsed_data['excited_states'],
                                       group='D2h',
                                       orientation=(1, 0, 0),
                                       orientation2=(0, 1, 0),
                                       )

# Analysis of adiabatic states
print('\nSymmetry analysis of adiabatic states')

for i, state in enumerate(parsed_data['excited_states']):
    sym_lab = symmetry_measures[i][0]
    print('State {}: {:3}  {:6.2f} {}'.format(i+1, sym_lab, state['excitation_energy'], state['excitation_energy_units']))


# select states to diabatize
# CT
list_diabatic = [3, 4, 5, 6]

# LE
# list_diabatic = [1, 2, 7, 8]

# all
# list_diabatic = [1, 2, 3, 4, 5, 6, 7, 8]

# CIS qchem input
qc_input = QchemInput(dimer,
                      exchange='hf',
                      basis='6-311g',
                      cis_n_roots=8,
                      cis_convergence=8,
                      cis_singlets=True,
                      cis_triplets=False,
                      cis_ampl_anal=True,
                      loc_cis_ov_separate=False,
                      namd_nsurfaces=0,
                      boys_cis_numstate=len(list_diabatic),  # use boys method
                      cis_diabath_decompose=True,
                      localized_diabatization=list_diabatic,
                      extra_rem_keywords={'dft_d': 'd3'},
                      RPA=False,
                      mem_static=900,
                      set_iter=30
                      )

parsed_data = get_output_from_qchem(qc_input,
                                    processors=4,
                                    force_recalculation=True,
                                    parser=basic_cis,
                                    store_full_output=True,
                                    )

# parsed_data = basic_cis(parsed_data)
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


print('\nDiabatic states dimer\n---------------------')

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
                      exchange='hf',
                      basis='6-311g',
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

parsed_data, ee = get_output_from_qchem(qc_input,
                                        processors=14,
                                        force_recalculation=False,
                                        parser=basic_cis,
                                        return_electronic_structure=True
                                        )

symmetry_measures = get_state_symmetry(ee,
                                       parsed_data['excited_states'],
                                       group='D2h',
                                       orientation=(0, 0, 1),
                                       orientation2=(1, 0, 0),
                                       )


print('\nAdiabatic states monomer\n------------------------')
for i, state in enumerate(parsed_data['excited_states']):
    print('\nState {}'.format(i+1))
    print('Transition DM: ', state['transition_moment'])
    print('Energy: ', state['excitation_energy'])
    print('Symmetry: ', symmetry_measures[i][0])
    print(' Origin  Target   Amplitude')
    for j, conf in enumerate(state['configurations']):
        print('{:^9} {:^7} {:8.3f}'.format(conf['origin'], conf['target'], conf['amplitude']))

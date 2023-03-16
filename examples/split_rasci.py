# splitting the calculation in Optimization + Reference + 2 x Post HF (RASCI)
from pyqchem import get_output_from_qchem, QchemInput
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.parsers.parser_rasci import parser_rasci
from pyqchem.structure import Structure


def print_formatted_data(parsed_data):
    for i, state in enumerate(parsed_data['excited_states']):
        print('\n - State {} -'.format(i + 1))
        if i > 0:
            print('Transition DM: ', state['transition_moment'])
        print('Excitation energy: {:5.2f} eV'.format(state['excitation_energy']))
        print('   Alpha  Beta   Amplitude')
        for j, conf in enumerate(state['configurations']):
            print('   {:^5}  {:^4} {:8.3f}'.format(conf['alpha'], conf['beta'], conf['amplitude']))


# define molecule (H2O)
molecule = Structure(coordinates=[[0.0,  1.0,  0.0],
                                  [0.0,  0.0,  1.0],
                                  [0.0,  0.0, -1.0]],
                     symbols=['O', 'H', 'H'],
                     charge=0,
                     multiplicity=1)

print('Initial structure')
print(molecule)

# optimize geometry
qc_input = QchemInput(molecule,
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

opt_molecule = parsed_data['optimized_molecule']

print('Optimized structure')
print(opt_molecule)

# Compute the HF reference
qc_input = QchemInput(opt_molecule,
                      jobtype='sp',
                      exchange='hf',
                      basis='sto-3g',
                      )

data_ref, ee_reference = get_output_from_qchem(qc_input,
                                               processors=4,
                                               return_electronic_structure=True)

# Compute post HF (RASCI) skipping SCF
qc_input = QchemInput(opt_molecule,
                      jobtype='sp',
                      exchange='hf',
                      correlation='rasci',
                      max_scf_cycles=0,
                      basis=ee_reference['basis'],
                      scf_guess=ee_reference['coefficients'],
                      ras_roots=3,
                      ras_elec=2,
                      ras_act=2)

parsed_data = get_output_from_qchem(qc_input,
                                    parser=parser_rasci,
                                    processors=4)

print('\nRASCI states (Active space 1)\n' + '-'*31)
print_formatted_data(parsed_data)

qc_input = QchemInput(opt_molecule,
                      jobtype='sp',
                      exchange='hf',
                      correlation='rasci',
                      max_scf_cycles=0,
                      basis=ee_reference['basis'],
                      scf_guess=ee_reference['coefficients'],
                      ras_roots=3,
                      ras_elec=2,
                      ras_act=3)

parsed_data = get_output_from_qchem(qc_input,
                                    parser=parser_rasci,
                                    processors=4)

print('\nRASCI states (Active space 2)\n' + '-'*31)
print_formatted_data(parsed_data)

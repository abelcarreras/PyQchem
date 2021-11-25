# Use of basis2 to use a guess from a previous calculation with a different basis set
from pyqchem import get_output_from_qchem, Structure, QchemInput
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.parsers.basic import basic_parser_qchem


# define molecule
coordinates = [[0.0,  1.0,  0.0],
               [0.0,  0.0,  1.0],
               [0.0,  0.0, -1.0]]

symbols = ['O', 'H', 'H']

molecule = Structure(coordinates=coordinates,
                     symbols=symbols,
                     charge=0,
                     multiplicity=1)

print('Initial structure')
print(molecule)

# Initial calculation with STO-3G basis
qc_input = QchemInput(molecule,
                      jobtype='opt',
                      exchange='hf',
                      basis='sto-3g',
                      )

parsed_data, ee = get_output_from_qchem(qc_input,
                                        processors=4,
                                        parser=basic_optimization,
                                        force_recalculation=False,
                                        read_fchk=True)

opt_molecule = parsed_data['optimized_molecule']

print('Optimized structure')
print(opt_molecule)
print('Final energy:', parsed_data['energy'])

# Precise calculation with larger 6-31G basis using previous MO as guess
qc_input = QchemInput(opt_molecule,
                      jobtype='sp',
                      exchange='hf',
                      basis='6-31g',
                      basis2=ee['basis'],
                      scf_guess=ee['coefficients']
                      )


parsed_data = get_output_from_qchem(qc_input,
                                    processors=4,
                                    parser=basic_parser_qchem,
                                    force_recalculation=True)

print('results: ', parsed_data)
from pyqchem import Structure, QchemInput, get_output_from_qchem
from pyqchem.qchem_core import redefine_calculation_data_filename
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.parsers.parser_cis_basic import basic_cis

redefine_calculation_data_filename('soc_benzo.pkl')

coordinates = [[0.0000000,	0.0000000,  0.6946470],
               [0.0000000,	0.0000000,  1.9169830],
               [0.0000000,	1.4037910,  0.0947360],
               [0.0000000, -1.4037910,	0.0947360],
               [0.0000000,	2.3981510,  1.0892960],
               [0.0000000, -2.3981510,	1.0892960],
               [0.0000000,	1.8577290, -1.2325690],
               [0.0000000, -1.8577290, -1.2325690],
               [0.0000000,	3.7488810,  0.7875230],
               [0.0000000, -3.7488810,	0.7875230],
               [0.0000000,	3.2140210, -1.5415030],
               [0.0000000, -3.2140210, -1.5415030],
               [0.0000000,	4.1665510, -0.5366730],
               [0.0000000, -4.1665510, -0.5366730],
               [0.0000000,	2.0768750,  2.1198300],
               [0.0000000, -2.0768750,	2.1198300],
               [0.0000000,	1.1928770, -2.0639740],
               [0.0000000, -1.1928770, -2.0639740],
               [0.0000000,	4.4784890,  1.5866060],
               [0.0000000, -4.4784890,	1.5866060],
               [0.0000000,	3.5236900, -2.5782740],
               [0.0000000, -3.5236900, -2.5782740],
               [0.0000000,	5.2210960, -0.7809190],
               [0.0000000, -5.2210960, -0.7809190]]

symbols = ['C', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
           'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']

molecule = Structure(coordinates=coordinates,
                     atomic_elements=symbols,
                     charge=0,
                     multiplicity=1)

qc_input = QchemInput(molecule,
                      jobtype='opt',
                      exchange='b3lyp',
                      basis='sto-3g',
                      geom_opt_tol_gradient=300,
                      geom_opt_tol_energy=100,
                      geom_opt_tol_displacement=1200,
                      geom_opt_max_cycles=50,
                      )

parsed_data = get_output_from_qchem(qc_input,
                                    processors=4,
                                    parser=basic_optimization,
                                    )

opt_molecule = parsed_data['optimized_molecule']

qc_input = QchemInput(opt_molecule,
                      jobtype='sp',
                      exchange='b3lyp',
                      basis='sto-3g',
                      cis_n_roots=5,
                      cis_singlets=True,
                      cis_triplets=True,
                      )

parsed_data, electronic_structure = get_output_from_qchem(qc_input,
                                                          processors=4,
                                                          parser=basic_cis,
                                                          store_full_output=True
                                                          )

print(parsed_data)

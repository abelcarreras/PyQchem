from pyqchem import Structure, QchemInput, get_output_from_qchem
from pyqchem.qchem_core import redefine_calculation_data_filename
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.parsers.parser_cis import basic_cis
from pyqchem.parsers.parser_rasci import parser_rasci
from pyqchem.basis import get_basis_from_ccRepo
from pyqchem.file_io import build_fchk
from pyqchem.errors import OutputError

redefine_calculation_data_filename('soc_benzo_2_ras.pkl')

coordinates = [[-0.49800797,  0.41832669,  0.00000000],
               [ 0.72930903,  0.41832669,  0.00000000],
               [-1.31920411,  1.72110714,  0.00000000],
               [-2.71298836,  1.66723263,  0.00000000],
               [-0.66872148,  2.95532032,  0.00026620],
               [-3.45673426,  2.84793521, -0.00001025],
               [-1.41217483,  4.13550038, -0.00026053],
               [ 0.43007141,  2.99754789,  0.00046401],
               [-2.80629666,  4.08169051, -0.00073971],
               [-4.55559312,  2.80544226,  0.00024612],
               [-0.89955468,  5.10831249, -0.00038580],
               [-3.39213342,  5.01230158, -0.00037769],
               [-1.31920411, -0.88445377,  0.00005409],
               [-2.71298836, -0.83057926,  0.00005185],
               [-0.66872148, -2.11866695, -0.00016087],
               [-3.45673426, -2.01128184,  0.00011112],
               [-1.41217483, -3.29884699,  0.00041485],
               [ 0.43007141, -2.16089453, -0.00035693],
               [-2.80629666, -3.24503710,  0.00089180],
               [-4.55559312, -1.96878890, -0.00014701],
               [-0.89955468, -4.27165909,  0.00058051],
               [-3.39213342, -4.17564818,  0.00056841],
               [-3.43149206,  0.53153375, -0.00072728],
               [-4.06231571,  0.55638246, -0.86463934],
               [-4.06411994,  0.55641830,  0.86186343]]

symbols = ['C', 'O', 'C', 'C', 'C', 'C', 'C', 'H', 'C', 'H', 'H', 'H',
           'C', 'C', 'C', 'C', 'C', 'H', 'C', 'H', 'H', 'H', 'C', 'H', 'H']

molecule = Structure(coordinates=coordinates,
                     symbols=symbols,
                     charge=0,
                     multiplicity=1)

basis_custom_repo = get_basis_from_ccRepo(molecule, 'cc-pVDZ')

qc_input = QchemInput(molecule,
                      jobtype='opt',
                      exchange='b3lyp',
                      basis='sto-3g',
                      geom_opt_tol_gradient=300,
                      geom_opt_tol_energy=100,
                      geom_opt_tol_displacement=1200,
                      geom_opt_max_cycles=50,
                      # n_frozen_core=0,
                      )

try:
    parsed_data, electronic_structure = get_output_from_qchem(qc_input,
                                                              processors=4,
                                                              parser=basic_optimization,
                                                              read_fchk=True,
                                                              store_full_output=True
                                                              )
except OutputError as e:
    print(e.error_lines)
    exit()

opt_molecule = parsed_data['optimized_molecule']

print('Optimized molecule')
print(opt_molecule.get_xyz())

basis_custom_repo = get_basis_from_ccRepo(molecule, 'cc-pVTZ')

qc_input = QchemInput(opt_molecule,
                      jobtype='sp',
                      exchange='b3lyp',
                      basis='sto-3g',
                      cis_n_roots=5,
                      cis_singlets=True,
                      cis_triplets=True,
                      calc_soc=1,
                      n_frozen_core=0,
                      )

print(qc_input.get_txt())

try:
    output, electronic_structure = get_output_from_qchem(qc_input,
                                                         processors=4,
                                                         # parser=basic_cis,
                                                         read_fchk=True,
                                                         store_full_output=True,
                                                         )
except OutputError as e:
    print(e.error_lines)
    exit()

print(output)
open('benzo_2.fchk', 'w').write(build_fchk(electronic_structure))

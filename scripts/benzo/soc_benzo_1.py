from pyqchem import Structure, QchemInput, get_output_from_qchem
from pyqchem.qchem_core import redefine_calculation_data_filename
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.parsers.parser_cis import basic_cis
from pyqchem.basis import get_basis_from_ccRepo
from pyqchem.file_io import build_fchk
from pyqchem.errors import OutputError

redefine_calculation_data_filename('soc_benzo.pkl')

coordinates = [[ 0.0000117181, 1.2100956095, 0.0061696737],
               [-0.0001093713, 2.4729755203,-0.0021015137],
               [ 1.3279168762, 0.4232986337,-0.0171965889],
               [-1.3278944069, 0.4235500050, 0.0399196965],
               [ 2.4807806771, 1.0511169801, 0.5180739582],
               [-2.4803868609, 1.0433998845,-0.5053843149],
               [ 1.4563127611,-0.8547043712,-0.6143971176],
               [-1.4567343496,-0.8455806492, 0.6556710511],
               [ 3.7292748613, 0.4051619446, 0.4862164494],
               [-3.7288476489, 0.3979028369,-0.4649885951],
               [ 2.7089790351,-1.4969333402,-0.6545726840],
               [-2.7093881412,-1.4872657442, 0.7042621010],
               [ 3.8442224305,-0.8731993475,-0.0989159605],
               [-3.8442208947,-0.8718029198, 0.1386676937],
               [ 2.3775425507, 2.0542967591, 0.9561195989],
               [-2.3769578269, 2.0399945699,-0.9581358524],
               [ 0.5782061601,-1.3390618566,-1.0619687011],
               [-0.5789352758,-1.3233795615, 1.1109066649],
               [ 4.6141668645, 0.8967520518, 0.9128524252],
               [-4.6134620945, 0.8832195307,-0.8993416880],
               [ 2.8007546220,-2.4859078357,-1.1237708580],
               [-2.8014425746,-2.4692672288, 1.1878232884],
               [ 4.8189374990,-1.3790371039,-0.1285568929],
               [-4.8189475688,-1.3771762516, 0.1749331354]]

symbols = ['C', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
           'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']

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
                                                         store_full_output=True
                                                         )
except OutputError as e:
    print(e.error_lines)
    exit()

open('benzo_1.fchk', 'w').write(build_fchk(electronic_structure))
print(output)

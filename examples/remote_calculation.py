# example of remote calculation
from pyqchem import get_output_from_qchem, QchemInput
from pyqchem.structure import Structure
from pyqchem.parsers.basic import basic_parser_qchem


remote_data_pc = {'hostname': '111.222.333.44',  # IP remote computer
                  'port': 22,                    # ssh port
                  'username': 'user',            # user in remote computer
                  'password': '**************',  # password to remote computer
                  'allow_agent': False,
                  'look_for_keys': False,
                  'precommand': ['module load qchem/qchem_group'], # commands (if) needed to load qchem in remote computer
                  'remote_scratch': '/home/user/'}  # scratch directory in remote computer

# molecule
ethene = [[0.0,  0.0000,   0.65750],
          [0.0,  0.0000,  -0.65750],
          [0.0,  0.92281,  1.22792],
          [0.0, -0.92281,  1.22792],
          [0.0, -0.92281, -1.22792],
          [0.0,  0.92281, -1.22792]]

symbols = ['C', 'C', 'H', 'H', 'H', 'H']

# create molecule
molecule = Structure(coordinates=ethene,
                     symbols=symbols,
                     charge=0,
                     multiplicity=1)

# create Q-Chem input
qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      basis='sto-3g')

print(qc_input.get_txt())

# get data from Q-Chem calculation
output = get_output_from_qchem(qc_input,
                               processors=1,
                               force_recalculation=True,
                               remote=remote_data_pc,
                               scratch='/Users/abel/Scratch/export/',
                               parser=basic_parser_qchem
                               )

# Get data
print('OUTPUT')
print(output)

print('the energy is: ', output['scf_energy'])
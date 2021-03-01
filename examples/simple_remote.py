from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input
from pyqchem.parsers.parser_frequencies import basic_frequencies
from pyqchem.structure import Structure


remote_data = {'hostname': '111.222.333.44', # the remote machine host
               'port': 22,                   # port for SSH connection
               'username': 'user',           # your username
               'password': 'thepassword',    # your password
               'allow_agent': False,
               'look_for_keys': False,
               'precommand': ['module load qchem']}  # some prepend lines you may need to load qchem

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

# calculation
qc_input = create_qchem_input(molecule,
                              jobtype='freq',
                              exchange='hf',
                              basis='sto-3g')

out, electronic_structure = get_output_from_qchem(qc_input,
                                                  processors=4,
                                                  read_fchk=True,
                                                  force_recalculation=True,
                                                  remote=remote_data,  # Set remote data
                                                  parser=basic_frequencies
                                                  )

# Show output
print(out)

# Show Electronic structure data
print(electronic_structure)
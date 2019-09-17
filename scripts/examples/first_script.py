from pyqchem.qchem_core import create_qchem_input, get_output_from_qchem
from pyqchem.structure import Structure
from pyqchem.parsers.basic import basic_parser_qchem


ethene = [[0.0,  0.0000,   0.65750],
          [0.0,  0.0000,  -0.65750],
          [0.0,  0.92281,  1.22792],
          [0.0, -0.92281,  1.22792],
          [0.0, -0.92281, -1.22792],
          [0.0,  0.92281, -1.22792]]

symbols = ['C', 'C', 'H', 'H', 'H', 'H']

# create molecule
molecule = Structure(coordinates=ethene,
                     atomic_elements=symbols,
                     charge=0,
                     multiplicity=1)

# create Q-Chem input
qc_input = create_qchem_input(molecule,
                              jobtype='sp',
                              exchange='hf',
                              basis='sto-3g')

# get data from Q-Chem calculation
output, err = get_output_from_qchem(qc_input,
                                    processors=1,
                                    force_recalculation=True,
                                    # parser=basic_parser_qchem
                                    )

# Get data
print('OUTPUT')
print(output)
print('ERROR')
print(err)
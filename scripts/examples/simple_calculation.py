from pyqchem import get_output_from_qchem, QchemInput, Structure
from pyqchem.parsers.basic import basic_parser_qchem

molecule = Structure(coordinates=[[0.0, 0.0, 0.0],
                                  [0.0, 0.0, 0.9]],
                     symbols=['O', 'H'],
                     charge=-1,
                     multiplicity=1)

# create qchem input
qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      basis='6-31G',
                      unrestricted=True)

# calculate and parse qchem output

data = get_output_from_qchem(qc_input,
                             processors=4,
                             parser=basic_parser_qchem,
                             )

print(data)
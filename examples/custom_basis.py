# Use of custom basis set obtained from ccRepo
from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.structure import Structure
from pyqchem.parsers.basic import basic_parser_qchem
from pyqchem.basis import get_basis_from_ccRepo
from pyqchem.basis import get_basis_from_BSE


# create molecule
molecule = Structure(coordinates=[[0.0, 0.0, 0.0000],
                                  [0.0, 0.0, 1.5811]],
                     symbols=['Se', 'H'],
                     charge=-1,
                     multiplicity=1)

# Standard input using 6-31G basis
qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      basis='6-31G')

# print input
print(qc_input.get_txt())

output, electronic_structure = get_output_from_qchem(qc_input,
                                                     processors=4,
                                                     force_recalculation=False,
                                                     read_fchk=True,
                                                     parser=basic_parser_qchem,
                                                     store_full_output=False)


print('scf_energy 6-31G: ', output['scf_energy'])

basis_custom = electronic_structure['basis']

# Standard input using a custom basis obtained from previous calculation
qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      basis=basis_custom)

output = get_output_from_qchem(qc_input,
                               processors=4,
                               force_recalculation=False,
                               parser=basic_parser_qchem
                               )

print('scf_energy (custom basis: 6-31G): ', output['scf_energy'])

# Get custom basis from ccRepo online repository
basis_custom_repo = get_basis_from_ccRepo(molecule, 'cc-pVTZ')

qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      basis=basis_custom_repo)


output = get_output_from_qchem(qc_input,
                               processors=4,
                               force_recalculation=False,
                               parser=basic_parser_qchem
                               )

print('scf_energy (custom basis: cc-pVTZ): ', output['scf_energy'])

basis_custom_repo = get_basis_from_BSE(molecule, '4zapa-nr', if_missing='6-31g')

qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      basis=basis_custom_repo)


output = get_output_from_qchem(qc_input,
                               processors=4,
                               force_recalculation=False,
                               parser=basic_parser_qchem
                               )

print('scf_energy (custom basis: coemd): ', output['scf_energy'])

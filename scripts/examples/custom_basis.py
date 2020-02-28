from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.structure import Structure
from pyqchem.parsers.basic import basic_parser_qchem
from pyqchem.basis import get_basis_from_ccRepo


# create molecule
molecule = Structure(coordinates=[[0.0, 0.0, 0.0000],
                                  [0.0, 0.0, 1.5811]],
                     atomic_elements=['Se', 'H'],
                     charge=-1,
                     multiplicity=1)

# Standard input using 6-31G basis
qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      basis='6-31G')

# print input
print(qc_input.get_txt())

output, err, electronic_structure = get_output_from_qchem(qc_input,
                                                          processors=4,
                                                          force_recalculation=False,
                                                          read_fchk=True,
                                                          parser=basic_parser_qchem,
                                                          store_full_output=False)


print('scf_energy 6-31G: ', output['scf energy'])

basis_custom = electronic_structure['basis']

# Standard input using a custom basis obtained from previous calculation
qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      basis=basis_custom)

output, _ = get_output_from_qchem(qc_input,
                                    processors=4,
                                    force_recalculation=False,
                                    parser=basic_parser_qchem
                                    )

print('scf_energy (custom basis: 6-31G): ', output['scf energy'])

# Get custom basis from ccRepo online repository
basis_custom_repo = get_basis_from_ccRepo(molecule, 'cc-pVTZ')

qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      basis=basis_custom_repo)


output, _ = get_output_from_qchem(qc_input,
                                  processors=4,
                                  force_recalculation=False,
                                  parser=basic_parser_qchem
                                  )

print('scf_energy (custom basis: cc-pVTZ): ', output['scf energy'])
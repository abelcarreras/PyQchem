import numpy as np

from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.structure import Structure
from pyqchem.parsers.parser_rasci import rasci as rasci_parser
from pyqchem.basis import get_basis_from_ccRepo

# create molecule
molecule = Structure(coordinates=[[0.0, 0.0, 0.0000],
                                  [0.0, 0.0, 1.5811]],
                     atomic_elements=['Se', 'H'],
                     charge=-1,
                     multiplicity=1)


# create Q-Chem input
qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      correlation='rasci',
                      ras_act=3,
                      ras_elec=2,
                      ras_occ=16,
                      basis='sto-3g',
                      ras_roots=2,
                      ras_do_hole=False,
                      ras_do_part=False)

print(qc_input.get_txt())


# get data from Q-Chem calculation
output, err, electronic_structure = get_output_from_qchem(qc_input,
                                                          processors=4,
                                                          force_recalculation=False,
                                                          read_fchk=True,
                                                          parser=rasci_parser,
                                                          store_full_output=True)


energies = [state['excitation_energy'] for state in output['excited states rasci']]

print('scf_energy: ', output['scf energy'])
print(energies)
# print(electronic_structure['basis'])

basis_set = get_basis_from_ccRepo(molecule, 'cc-pVTZ')

# create Q-Chem input
qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      correlation='rasci',
                      ras_act=3,
                      ras_elec=2,
                      ras_occ=16,
                      # basis='6-31g',
                      basis=basis_set,
                      ras_roots=2,
                      ras_do_hole=False,
                      ras_do_part=False)

print('---------------')
#print(qc_input.get_txt())


output, err, electronic_structure = get_output_from_qchem(qc_input,
                                                          processors=4,
                                                          force_recalculation=False,
                                                          read_fchk=True,
                                                          parser=rasci_parser,
                                                          store_full_output=True
                                                          )
print(output)


energies = [state['excitation_energy'] for state in output['excited states rasci']]

print('scf_energy: ', output['scf energy'])
print(energies)

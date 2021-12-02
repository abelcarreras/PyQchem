# COmpute the wave function symmetry of excited states computed with CIS and RASCI
from pyqchem.utils import get_plane
from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.structure import Structure
from pyqchem.parsers.parser_rasci import parser_rasci
from pyqchem.parsers.parser_cis import basic_cis
from pyqchem.symmetry import get_state_symmetry
import numpy as np


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

point_group = molecule.get_point_symmetry()
print('Point group: {}\n'.format(point_group))

# create Q-Chem input
qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      correlation='rasci',
                      ras_act=6,
                      ras_elec=6,
                      basis='6-311(2+,2+)G**',
                      ras_roots=10,
                      ras_do_hole=False,
                      ras_do_part=False)


# get data from Q-Chem calculation
output, electronic_structure = get_output_from_qchem(qc_input,
                                                     processors=4,
                                                     force_recalculation=False,
                                                     read_fchk=True,
                                                     parser=parser_rasci,
                                                     store_full_output=True)


# get plane from coordinates to determine orientation using the plane of the molecule
coordinates_mat = np.array(electronic_structure['structure'].get_coordinates())
center, normal = get_plane(coordinates_mat)


symmetry_measures = get_state_symmetry(electronic_structure,
                                       output['excited_states'],
                                       center=center,  # if not provided use center of electronic charge
                                       orientation=normal,  # if not provided auto-search best orientation
                                       group=point_group,
                                       extra_print=False
                                       )

energies = [state['excitation_energy'] for state in output['excited_states']]

print('Symmetry of RASCI excited states\n--------------------------------')
for i, (energy, state) in enumerate(zip(energies, symmetry_measures)):
    print('{:2} {:6} {:5.3f}'.format(i+1, state[0], state[1]) + '  {:5.3f} eV'.format(energy))


# J. Phys. Chem. 1992, 96, 26, 10756–10768
qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      basis='6-311(2+,2+)G**',
                      cis_n_roots=10)

# get data from Q-Chem calculation
output, electronic_structure = get_output_from_qchem(qc_input,
                                                     processors=4,
                                                     force_recalculation=False,
                                                     read_fchk=True,
                                                     parser=basic_cis,
                                                     store_full_output=True)

symmetry_measures = get_state_symmetry(electronic_structure,
                                       output['excited_states'],
                                       center=center,  # if not provided use center of electronic charge
                                       orientation=normal,  # if not provided auto-search best orientation
                                       group=point_group,
                                       extra_print=False
                                       )

energies = [state['excitation_energy'] for state in output['excited_states']]

print('\nSymmetry of CIS excited states\n--------------------------------')
for i, (energy, state) in enumerate(zip(energies, symmetry_measures)):
    print('{:2} {:6} {:5.3f}'.format(i+1, state[0], state[1]) + '  {:5.3f} eV'.format(energy))

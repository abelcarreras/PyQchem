# Simple single point calculation
from pyqchem import get_output_from_qchem, QchemInput, Structure
from pyqchem.parsers.basic import basic_parser_qchem
import numpy as np

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
                             force_recalculation=True,
                             parser=basic_parser_qchem,
                             )


print('scf energy', data['scf_energy'], 'H')

print('Orbital energies (H)')
print('  {:^10} {:^10}'.format('alpha', 'beta'))
for a, b in zip(data['orbital_energies']['alpha'], data['orbital_energies']['beta']):
    print('{:10.3f} {:10.3f}'.format(a, b))

print('dipole moment:', data['multipole']['dipole_moment'], data['multipole']['dipole_units'])
print('quadrupole moment ({}):\n'.format(data['multipole']['quadrupole_units']),
      np.array(data['multipole']['quadrupole_moment']))
print('octopole moment ({}):\n'.format(data['multipole']['octopole_units']),
      np.array(data['multipole']['octopole_moment']))

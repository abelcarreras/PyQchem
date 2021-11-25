# Example of use of solvation (PCM) in pyqchem
from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input
from pyqchem.parsers.basic import basic_parser_qchem
from pyqchem.structure import Structure
import numpy as np


# define molecule
coordinates = [[ 0.00000, 1.40272, 0.00000],
               [ 0.00000, 2.49029, 0.00000],
               [-1.21479, 0.70136, 0.00000],
               [-2.15666, 1.24515, 0.00000],
               [-1.21479, -0.70136, 0.00000],
               [-2.15666, -1.24515, 0.00000],
               [ 0.00000, -1.40272, 0.00000],
               [ 0.00000, -2.49029, 0.00000],
               [ 1.21479, -0.70136, 0.00000],
               [ 2.15666, -1.24515, 0.00000],
               [ 1.21479, 0.70136, 0.00000],
               [ 2.15666, 1.24515, 0.00000]]

symbols = ['C', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'C', 'H']

molecule = Structure(coordinates=coordinates,
                     symbols=symbols,
                     charge=0,
                     multiplicity=1)

print('Initial structure')
print(molecule)

# optimization
qc_input = create_qchem_input(molecule,
                              jobtype='sp',
                              exchange='hf',
                              basis='sto-3g',
                              unrestricted=True,
                              solvent_method='pcm',
                              solvent_params={'Dielectric': 8.93},  # Cl2CH2
                              pcm_params={'Theory': 'CPCM',
                                          'Method': 'SWIG',
                                          'Solver': 'Inversion',
                                          'Radii': 'Bondi'}
                              )


data, electronic_structure = get_output_from_qchem(qc_input,
                                                   force_recalculation=True,
                                                   processors=4,
                                                   parser=basic_parser_qchem,
                                                   read_fchk=True)


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

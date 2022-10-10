# Example of the calculation of the absorption spectrum for methyl peroxy radical
from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.parsers.parser_frequencies import basic_frequencies
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.tools.duschinsky import get_duschinsky
from pyqchem import Structure
import numpy as np
import matplotlib.pyplot as plt
from pyqchem.units import AU_TO_EV


molecule = Structure(coordinates=[[ 1.004123,  -0.180454,   0.000000],
                                  [-0.246002,   0.596152,   0.000000],
                                  [-1.312366,  -0.230256,   0.000000],
                                  [ 1.810765,   0.567203,   0.000000],
                                  [ 1.036648,  -0.805445,  -0.904798],
                                  [ 1.036648,  -0.805445,   0.904798]],
                     symbols=['C', 'O', 'O', 'H', 'H', 'H'],
                     multiplicity=2)

basis_set = '6-31+G*'
n_state = 1  # excited state number

# Optimization of ground state geometry
qc_input = QchemInput(molecule,
                      jobtype='opt',
                      exchange='b3lyp',
                      basis=basis_set,
                      )

output = get_output_from_qchem(qc_input,
                               processors=4,
                               force_recalculation=False,
                               parser=basic_optimization,
                               )

print('GS optimized')
gs_energy = output['energy']
print(output['optimized_molecule'])


# Ground state frequencies calculation
qc_input = QchemInput(output['optimized_molecule'],
                      jobtype='freq',
                      exchange='b3lyp',
                      basis=basis_set,
                      )


gs_output = get_output_from_qchem(qc_input,
                                  processors=4,
                                  force_recalculation=False,
                                  parser=basic_frequencies,
                                  store_full_output=True,
                                  )

print('GS frequencies done')
print('frequency         displacements')
for mode in gs_output['modes']:
    print('{:10}'.format(mode['frequency']), mode['displacement'])

# Optimization of the first excited state state geometry
qc_input = QchemInput(molecule,
                      jobtype='opt',
                      exchange='b3lyp',
                      basis=basis_set,
                      cis_n_roots=10,
                      cis_singlets=True,
                      cis_triplets=False,
                      #cis_convergence=6,
                      cis_state_deriv=n_state,
                      extra_rem_keywords={'XC_GRID': '000075000302'}
                      )


output = get_output_from_qchem(qc_input,
                               processors=4,
                               force_recalculation=False,
                               store_full_output=True,
                               parser=basic_optimization,
                               )
print('ES optimized')
es_energy = output['energy']
print(output['optimized_molecule'])

# Excited state frequencies calculation
qc_input = QchemInput(output['optimized_molecule'],
                      jobtype='freq',
                      exchange='b3lyp',
                      basis=basis_set,
                      cis_n_roots=10,
                      cis_singlets=True,
                      cis_triplets=False,
                      #cis_convergence=6,
                      cis_state_deriv=n_state,
                      mem_static=200,
                      extra_rem_keywords={'XC_GRID': '000075000302'}
                      )

es_output = get_output_from_qchem(qc_input,
                                  processors=4,
                                  force_recalculation=False,
                                  parser=basic_frequencies,
                                  store_full_output=True,
                                  )

print('ES frequencies done')
print('frequency         displacements')
for mode in es_output['modes']:
    print('{:10}'.format(mode['frequency']), mode['displacement'])


excitation_energy = (es_energy - gs_energy) * AU_TO_EV
print('\nexcitation energy (state {}): {}'.format(n_state, excitation_energy))


# build duschinsky object
duschinsky = get_duschinsky(gs_output, es_output, n_max_modes=6)

# align structures along principal axis of inertia
duschinsky.align_coordinates()

# compute and print data
transitions = duschinsky.get_transitions(max_vib_origin=1, max_vib_target=2, excitation_energy=excitation_energy)

temperature = 300

print('\n {:^10} {:^10} {:^10}  '.format('Energy', 'FCF', 'Intensity'), 'Transition')
for s in transitions:
    print('{:10.5f} {:10.5f} {:10.5f}  '.format(s.energy, s.get_intensity(temperature), s.fcf), s.get_label())


def gaussian(x, s, m):
    return 1/np.sqrt(2*np.pi*s**2)*np.exp(-0.5*((x-m)/s)**2)


# plot spectrum
sigma = 0.01
min = transitions[0].energy
max = transitions[-1].energy

energies = np.linspace(min, max, 500)
intensities = np.zeros_like(energies)
for trans in transitions:
    if trans.get_intensity(1) > 0.001:
        for i, e in enumerate(energies):
            intensities[i] += gaussian(e, sigma, trans.energy)*trans.get_intensity(temperature)
    if trans.get_intensity(1) > 0.01:
        height = gaussian(0, sigma, 0) * trans.get_intensity(temperature)
        plt.text(trans.energy, height, trans.get_label(), rotation=80)

plt.xlabel('Energy [eV]')
plt.yticks([], [])
plt.plot(energies, intensities)
plt.show()

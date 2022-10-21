# Example of the calculation of the emission/absorption spectrum for methyl peroxy radical
# This includes the calculation of the Franck Condon Weithed Density (FCWD)
from pyqchem import get_output_from_qchem, Structure, QchemInput
from pyqchem.parsers.parser_frequencies import basic_frequencies
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.parsers.parser_cis import basic_cis
from pyqchem.parsers.basic import basic_parser_qchem
from pyqchem.tools.duschinsky import get_duschinsky
from pyqchem.tools.gaussian import gaussian
from pyqchem.tools.spectrum import get_fcwd
from pyqchem.units import AU_TO_EV, KB_EV
import numpy as np
import matplotlib.pyplot as plt


molecule = Structure(coordinates=[[ 1.004123,  -0.180454,   0.000000],
                                  [-0.246002,   0.596152,   0.000000],
                                  [-1.312366,  -0.230256,   0.000000],
                                  [ 1.810765,   0.567203,   0.000000],
                                  [ 1.036648,  -0.805445,  -0.904798],
                                  [ 1.036648,  -0.805445,   0.904798]],
                     symbols=['C', 'O', 'O', 'H', 'H', 'H'],
                     multiplicity=2)

basis_set = '6-31+G*'
n_state = 1  # excited state number (target)
temperature = 50
plot_marcus = False

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
print('frequency(cm-1)         displacements')
for mode in gs_output['modes']:
    print('{:10}'.format(mode['frequency']), mode['displacement'])

# Excited state with ground state geometry
qc_input = QchemInput(output['optimized_molecule'],
                      jobtype='sp',
                      exchange='b3lyp',
                      basis=basis_set,
                      cis_n_roots=10,
                      cis_singlets=True,
                      cis_triplets=False,
                      # cis_convergence=6,
                      extra_rem_keywords={'XC_GRID': '000075000302'}
                      )


output_es_gs = get_output_from_qchem(qc_input,
                                     processors=4,
                                     force_recalculation=False,
                                     store_full_output=True,
                                     parser=basic_cis
                                     )

es_gs_energy = output_es_gs['excited_states'][n_state-1]['total_energy']

# Optimization of the first excited state state geometry
qc_input = QchemInput(molecule,
                      jobtype='opt',
                      exchange='b3lyp',
                      basis=basis_set,
                      cis_n_roots=10,
                      cis_singlets=True,
                      cis_triplets=False,
                      # cis_convergence=6,
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
                      # cis_convergence=6,
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
print('frequency(cm-1)         displacements')
for mode in es_output['modes']:
    print('{:10}'.format(mode['frequency']), mode['displacement'])

# Ground state with excited state geometry
qc_input = QchemInput(output['optimized_molecule'],
                      jobtype='sp',
                      exchange='b3lyp',
                      basis=basis_set,
                      )

output_gs_es = get_output_from_qchem(qc_input,
                                     processors=4,
                                     force_recalculation=False,
                                     parser=basic_parser_qchem
                                     )


gs_es_energy = output_gs_es['scf_energy']


excitation_energy = (es_energy - gs_energy) * AU_TO_EV
print('\nexcitation energy (state {}): {:7.4f} eV'.format(n_state, excitation_energy))

# reorganization energy
reorganization = ((es_gs_energy - gs_energy) + (gs_es_energy - es_energy)) * AU_TO_EV
print('Reorganization energy (total): {:7.4f} eV'.format(reorganization))

# build Duschinsky object
duschinsky = get_duschinsky(gs_output, es_output, n_max_modes=6)

# align structures along principal axis of inertia
duschinsky.align_coordinates_pmi()

# compute and print vibrational transitions
transitions = duschinsky.get_transitions(max_vib_origin=3, max_vib_target=3,
                                         excitation_energy=excitation_energy,
                                         reorganization_energy=reorganization,
                                         add_hot_bands=False)

print('\n {:^10} {:^10}  {:^10} {:^10}'.format('Energy', 'FCF', 'Int. abs.',  'Int. em.', 'Transition'))
for s in transitions:
    print('{:10.5f} {:10.5f} {:10.5f} {:10.5f}   {}'.format(s.energy,
                                                            s.fcf,
                                                            s.get_intensity_absorption(temperature),
                                                            s.get_intensity_emission(temperature),
                                                            s.get_label()))

# plot spectrum
sigma = np.sqrt(2*KB_EV*temperature*reorganization)  # Marcus model for band amplitude

cutoff = 0.001
cutoff_labels = 0.02
n_points = 500

min = transitions[0].energy_emission-sigma
max = transitions[-1].energy_absorption+sigma

energies = np.linspace(min, max, n_points)
intensities_abs = np.zeros_like(energies)
intensities_em = np.zeros_like(energies)

for trans in transitions:
    # absorption
    if trans.get_intensity_absorption(1) > cutoff:
        for i, e in enumerate(energies):
            intensities_abs[i] += gaussian(e, sigma, trans.energy_absorption) * trans.get_intensity_absorption(temperature)
    if trans.get_intensity_absorption(1) > cutoff_labels:
        height = gaussian(0, sigma, 0) * trans.get_intensity_absorption(temperature)
        plt.text(trans.energy_absorption, height, trans.get_label(sign=1), rotation=80, color='C0')
        plt.stem(trans.energy_absorption, height*0.9, markerfmt='', linefmt='C0-')
    # emission
    if trans.get_intensity_emission(1) > cutoff:  # cutoff for plot
        for i, e in enumerate(energies):
            intensities_em[i] += gaussian(e, sigma, trans.energy_emission) * trans.get_intensity_emission(temperature)
    if trans.get_intensity_emission(1) > cutoff_labels:
        height = gaussian(0, sigma, 0) * trans.get_intensity_emission(temperature)
        plt.text(trans.energy_emission, height, trans.get_label(sign=0), rotation=80, color='C1')
        plt.stem(trans.energy_emission, height*0.9, markerfmt='', linefmt='C1-')

plt.xlabel('Energy [eV]')
plt.ylabel('Intensity [eV-1]')
plt.yticks([], [])
plt.title('Spectrum at {} K'.format(temperature))
plt.plot(energies, intensities_abs, label='absorption', color='C0')
plt.plot(energies, intensities_em, label='emission', color='C1')

# marcus
def marcus_emission(e, de, lmb):
    prefactor = 1/np.sqrt(np.pi*4*KB_EV*temperature*lmb)
    return prefactor * np.exp(-(de - e + lmb)**2/(4*KB_EV*temperature*lmb))

def marcus_absorption(e, de, lmb):
    prefactor = 1/np.sqrt(np.pi*4*KB_EV*temperature*lmb)
    return prefactor * np.exp(-(-de + e + lmb)**2/(4*KB_EV*temperature*lmb))

if plot_marcus:
    marcus_em = np.zeros_like(energies)
    marcus_abs = np.zeros_like(energies)
    for i, e in enumerate(energies):
        marcus_em[i] = marcus_emission(e, excitation_energy, reorganization)
        marcus_abs[i] = marcus_absorption(e, excitation_energy, reorganization)

    plt.plot(energies, marcus_em, '--', label='marcus emission', color='C0')
    plt.plot(energies, marcus_abs, '--', label='marcus absorption', color='C1')


print('integral absorption: ', np.trapz(intensities_abs, energies))  # should be clode to 1
print('integral emission: ', np.trapz(intensities_em, energies))  # should be clode to 1

if plot_marcus:
    print('integral marcus absorption: ', np.trapz(marcus_abs, energies))  # should be clode to 1
    print('integral marcus emission: ', np.trapz(marcus_em, energies))  # should be clode to 1

print('FCWD: ', get_fcwd(transitions, temperature))
if plot_marcus:
    print('FCWD Marcus: {}'.format(np.trapz(marcus_em * marcus_abs, energies)))

plt.legend()
plt.show()

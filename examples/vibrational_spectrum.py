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
temperature = 100
plot_marcus = False  # plot Marcus model
plot_lj = True  # plot Levich–Jortner model

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
reorganization = ((es_gs_energy - es_energy) + (gs_es_energy - gs_energy)) * AU_TO_EV
print('Reorganization energy (total): {:7.4f} eV'.format(reorganization))

# build Duschinsky object
duschinsky = get_duschinsky(gs_output, es_output, n_max_modes=6)

# align structures along principal axis of inertia
duschinsky.align_coordinates_pmi()

# print modes included in the calculation
m1, m2 = duschinsky.get_restricted_modes()
print('Normal modes included:')
print('GS: ', m1)
print('ES: ', m2)

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
sigma = np.sqrt(2*KB_EV*temperature*reorganization/2)  # Marcus model for peaks amplitude

cutoff = 0.001
cutoff_labels = 0.02
n_points = 500

max = transitions[0].energy_emission+sigma
min = transitions[-1].energy_absorption-sigma

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

print('\nintegral absorption: {:7.4f}'.format(np.trapz(intensities_abs, energies)))  # should be close to 1
print('integral emission: {:7.4f}'.format(np.trapz(intensities_em, energies)))  # should be close to 1
print('FCWD: {:7.4f} eV^-1'.format(get_fcwd(transitions, temperature)))


# Marcus
def marcus_absorption(e, de, lmb):
    prefactor = 1/np.sqrt(np.pi*4*KB_EV*temperature*lmb)
    return prefactor * np.exp(-(de - e + lmb)**2/(4*KB_EV*temperature*lmb))

def marcus_emission(e, de, lmb):
    prefactor = 1/np.sqrt(np.pi*4*KB_EV*temperature*lmb)
    return prefactor * np.exp(-(-de + e + lmb)**2/(4*KB_EV*temperature*lmb))

if plot_marcus:
    marcus_em = np.zeros_like(energies)
    marcus_abs = np.zeros_like(energies)
    for i, e in enumerate(energies):
        marcus_em[i] = marcus_emission(e, excitation_energy, reorganization/2)
        marcus_abs[i] = marcus_absorption(e, excitation_energy, reorganization/2)

    plt.plot(energies, marcus_em, '--', label='marcus emission', color='C0')
    plt.plot(energies, marcus_abs, '--', label='marcus absorption', color='C1')

    print('\nintegral marcus absorption: {:7.4f}'.format(np.trapz(marcus_abs, energies)))  # should be close to 1
    print('integral marcus emission: {:7.4f}'.format(np.trapz(marcus_em, energies)))  # should be close to 1
    print('FCWD Marcus: {:7.4f} eV^-1'.format(np.trapz(marcus_em * marcus_abs, energies)))


# Levich–Jortner
def lj_function(e, de, l_cl, huang_rhys, frequencies):
    sign = np.sign(de)

    CM_TO_NS = 29.9792558
    HBAR_PLANCK = 6.58210893126952e-7  # eV * ns

    angl_freqs = np.array(frequencies) * CM_TO_NS * 2 * np.pi  # cm-1 (ordinary) -> ns-1 (angular)

    s_eff = np.sum(huang_rhys)
    w_eff = np.sum(np.multiply(huang_rhys, angl_freqs)) / s_eff  # angular frequency

    e = np.array(e, dtype=float)
    fcwd_term = np.zeros_like(e)
    for m in range(10):
        fcwd_term += s_eff**m / np.math.factorial(m) * np.exp(-s_eff) * np.exp(
            -(de - e * sign + l_cl + m * HBAR_PLANCK * w_eff)**2 / (4 * KB_EV * temperature * l_cl))

    return 1.0 / (np.sqrt(4 * np.pi * KB_EV * temperature * l_cl)) * fcwd_term

def lj_absorption(e, de, lmb, huang_rhys, frequencies):
    return lj_function(e, de, lmb, huang_rhys, frequencies)
def lj_emission(e, de, lmb, huang_rhys, frequencies):
    return lj_function(e, -de, lmb, huang_rhys, frequencies)


if plot_lj:
    s_i, s_f = duschinsky.get_huang_rys()
    l_i, l_f = duschinsky.get_reorganization_energies()  # eV
    freq_i = duschinsky._modes_initial.get_frequencies()  # freq in cm-1
    freq_f = duschinsky._modes_final.get_frequencies()  # freq in cm-1

    print('\nHuang Rhys')
    for i, (ssi, ssf) in enumerate(zip(s_i, s_f)):
        print('{} {:8.4f} {:8.4f}'.format(i, ssi, ssf))

    print('\nReorganization energies (meV)')
    for i, (lli, llf) in enumerate(zip(l_i, l_f)):
        print('{} {:8.4f} {:8.4f}'.format(i, lli * 1000, llf * 1000))

    print('sum: {:8.4f} eV {:8.4f} eV\ntotal: {:8.4f} eV'.format(np.sum(l_i), np.sum(l_f), np.sum(l_i) + np.sum(l_f)))

    lj_em = np.zeros_like(energies)
    lj_abs = np.zeros_like(energies)

    for i, e in enumerate(energies):
        lj_em[i] = lj_emission(e, excitation_energy, reorganization/2, s_f, freq_f)
        lj_abs[i] = lj_absorption(e, excitation_energy, reorganization/2, s_i, freq_i)

    print('\nintegral marcus absorption: {:7.4f}'.format(np.trapz(lj_em, energies)))  # should be close to 1
    print('integral marcus emission: {:7.4f}'.format(np.trapz(lj_abs, energies)))  # should be close to 1
    print('FCWD LJ: {:7.4f} eV^-1'.format(np.trapz(lj_em * lj_abs, energies)))

    plt.plot(energies, lj_em, ':', label='LJ emission', color='C0')
    plt.plot(energies, lj_abs, ':', label='LJ absorption', color='C1')

plt.legend()
plt.show()

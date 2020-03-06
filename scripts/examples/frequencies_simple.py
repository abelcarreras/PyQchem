from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input
from pyqchem.parsers.parser_frequencies import basic_frequencies
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.structure import Structure

import numpy as np
import matplotlib.pyplot as plt


# define molecule
coordinates = [[0.0,  1.0,  0.0],
               [0.0,  0.0,  1.0],
               [0.0,  0.0, -1.0]]

symbols = ['O', 'H', 'H']

molecule = Structure(coordinates=coordinates,
                     atomic_elements=symbols,
                     charge=0,
                     multiplicity=1)

print('Initial structure')
print(molecule)

# optimization
qc_input = create_qchem_input(molecule,
                              jobtype='opt',
                              exchange='hf',
                              basis='sto-3g',
                              geom_opt_tol_gradient=300,
                              geom_opt_tol_energy=100,
                              geom_opt_coords=-1,
                              geom_opt_tol_displacement=1200)


parsed_data, electronic_structure = get_output_from_qchem(qc_input,
                                                          processors=4,
                                                          parser=basic_optimization,
                                                          force_recalculation=False,
                                                          read_fchk=True)


opt_molecule = parsed_data['optimized_molecule']

print('Optimized structure')
print(opt_molecule)

# frequencies calculation
qc_input = create_qchem_input(opt_molecule,
                              jobtype='freq',
                              exchange='hf',
                              basis='sto-3g',
                              scf_guess=electronic_structure['coefficients'])

parsed_data = get_output_from_qchem(qc_input,
                                         processors=4,
                                         force_recalculation=False,
                                         parser=basic_frequencies)


# print results
print('Normal modes\n')

for mode, freq in enumerate(parsed_data['frequencies']):

    force_constants = parsed_data['force_constants'][mode]

    print('mode:                      {}'.format(mode+1))
    print('frequency (cm-1):          {:10.2f}'.format(freq))
    print('force constant (mdyne/A):  {:10.5f}\n'.format(force_constants))


# Thermodynamics
def get_thermodynamics(T):
    h_bar = 6.62607015e-34  # Plank constant J * s
    kb = 1.38064852e-23  # Boltzmann constant J / K
    cm_to_hz = 29979245800  # Hz/cm-1
    na = 6.02214076e23  # Avogadro constant
    joule_to_cal = 0.239005

    free_energy = 0
    entropy = 0
    total_energy = 0
    for mode, freq in enumerate(parsed_data['frequencies']):
        freq_hz = freq * cm_to_hz
        fact = np.exp(-h_bar * freq_hz / (kb * T))
        free_energy += freq_hz * h_bar / 2 + kb * T * np.log(1.0 - fact)
        entropy += -kb * np.log(1.0 - fact) + h_bar * freq_hz / T * fact/(1-fact)
        total_energy += freq_hz * h_bar/2 + h_bar*freq_hz * fact/(1-fact)

    free_energy *= na * joule_to_cal  # cal/mol
    entropy *= na * joule_to_cal  # cal/mol.K
    total_energy *= na * joule_to_cal  # cal/mol

    return free_energy, entropy, total_energy


T = 298.15  # Kelvin
free_energy, entropy, total_energy = get_thermodynamics(T)

# U = F + TS
print('Thermodynamics at {} K'.format(T))
print('Vibrational free energy (F): {:.4f} cal/mol'.format(free_energy))
print('Vibrational entropy (S): {:.4f} cal/mol.K'.format(entropy))
print('Total vibrational energy (U): {:.4f} cal/mol'.format(total_energy))

free_energy = []
entropy = []
total_energy = []
temp_range = range(1, 300)
for T in temp_range:
    fe, en, te = get_thermodynamics(T)
    free_energy.append(fe)
    entropy.append(en)
    total_energy.append(te)

plt.figure(1)
plt.plot(entropy, label='Entropy')
plt.xlabel('Temperature [K]')
plt.ylabel('cal/mol.K')
plt.legend()

plt.figure(2)
plt.plot(temp_range, free_energy, label='Free energy')
plt.plot(temp_range, total_energy, label='Total energy')
plt.xlabel('Temperature [K]')
plt.ylabel('cal/mol')
plt.legend()

plt.show()

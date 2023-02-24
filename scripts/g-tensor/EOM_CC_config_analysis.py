#####################################
#          MODULES SELECTION
#####################################
# module load Python/3.7.4-Anaconda3-2019.10
import numpy as np
from tabulate import tabulate

#####################################
#            EOM COMPARISON
#####################################
EOM_input = '../Sven_EOM-MP2/vo_h2o5_2+_eomea_mp2_def2-TZVP.in.out'
threshold_excitation_energy = 10 # in eV, energy difference with respect to ground state

from g_take_EOM_states import get_EOM_type, get_SCF_energy, get_irreps_energies, \
    get_maximum_amplitude_orbitals, get_EOM_SOCC_values, second_smallest_number

EOM_type = get_EOM_type(EOM_input)
SCF_energy = get_SCF_energy(EOM_input)

print("------------------------")
print("    EOM-CC ANALYSIS ")
print("------------------------")
print("SCF   reference energy: ", SCF_energy)

if (EOM_type == 'EOMIP'):
    irreps, optimized_state_index, energies, excit_energies = get_irreps_energies(EOM_input)
    EOM_SOCC = get_EOM_SOCC_values(EOM_input, optimized_state_index)
    orbitals = get_maximum_amplitude_orbitals(EOM_input, EOM_type)
    presentation_list = []
    presentation_list.append(['Symmetry', 'Transition', 'Energy (eV)', 'Excitation energy (eV)', 'SOCC (cm-1)', \
                              'Occ. 1', 'Occ. 2', 'Virtual 1'])

elif (EOM_type == 'EOMEA'):
    irreps, optimized_state_index, energies, excit_energies = get_irreps_energies(EOM_input)
    EOM_SOCC = get_EOM_SOCC_values(EOM_input, optimized_state_index)
    orbitals = get_maximum_amplitude_orbitals(EOM_input, EOM_type)
    presentation_list = []
    presentation_list.append(['Symmetry', 'Transition', 'Energy (eV)', 'Excitation energy (eV)', 'SOCC (cm-1)', \
                              'Occ. 1', 'Virtual 1', 'Virtual 2'])

EOM_state = 0
transition = 0
EOM_excitation_energies_eV_list = []

for i in range(0, len(irreps)):
    symmetry = irreps[i][1]
    # transition = irreps[i][0]
    total_energy = np.round(float(energies[i]), 3)
    excitation_energy = np.round(float(excit_energies[i]), 2)
    SOCC = np.round(float(EOM_SOCC[i]), 2)
    # SOCC = 0

    if (excitation_energy <= threshold_excitation_energy):
        EOM_excitation_energies_eV_list.append(excit_energies[i])
        transition += 1
        presentation_list.append([symmetry, transition, total_energy, excitation_energy, SOCC, \
                                  orbitals[i*3], orbitals[i*3+1], orbitals[i*3+2]])
        EOM_state += 1

    if (i < len(irreps)-1) and (irreps[i][1] != irreps[i+1][1]):
        presentation_list.append(['---'])

second_smallest_excit_energy = second_smallest_number(excit_energies)
if (second_smallest_excit_energy > threshold_excitation_energy):
    print('Increase the threshold_excitation_energy')
    exit()

EOM_excitation_energies_eV = np.array(EOM_excitation_energies_eV_list, dtype=float)
EOM_excitation_energies = EOM_excitation_energies_eV / 27.211399
#
# print(tabulate(presentation_list, headers='firstrow'))

print('\n'.join([','.join(['{:^20}'.format(item) for item in row]) \
                 for row in (presentation_list)]))

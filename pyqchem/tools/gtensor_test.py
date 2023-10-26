__author__ = 'Antonio Cebreiro-Gallardo'

import numpy as np
from gtensor_parser import complete_procedure

#####################################
#            INPUT VALUES
#####################################
qchem_output = 'h2o_triplet.out'

ppm = 0  # 0) g-shifts in ppt 1) g-shifts in ppm
state_selection = 0  # 0) use "nstates" 1) use all states in "qchem_output"
nstates = [1, 2, 3]  # States to be included when state_selection = 0

#####################################
#     G-MATRIX CALCULATION
#####################################
g_shifts, totalstates = complete_procedure(qchem_output, nstates, state_selection, ppm)

print("---------------------")
print("  INPUT SECTION")
print("---------------------")
print("File selected: ", qchem_output)
print("Number of states in output: ", totalstates)
print("Selected states selected: ", nstates)

print(" ")
print("---------------------")
print(" RAS-CI RESULTS")
print("---------------------")
print('g-factor (x y z dimensions):')
print(np.round(g_shifts.real[0], 3), np.round(g_shifts.real[1], 3), np.round(g_shifts.real[2], 3))
print('')

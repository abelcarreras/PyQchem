# Obtain 2e integrals and compute J and K externally
from pyqchem import get_output_from_qchem, QchemInput, Structure
import numpy as np

molecule = Structure(coordinates=[[0.0, 0.0, 0.0],
                                  [0.0, 0.0, 0.9]],
                     symbols=['H', 'F'],
                     charge=0,
                     multiplicity=1)

# create qchem input
qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      basis='sto-3g',
                      extra_rem_keywords={'SAVE_AO_INTEGRALS': 1,
                                          'USE_NEW_PATH2': False,
                                          'scf_final_print': 2},
                      symmetry=False,
                      )

# calculate and parse qchem output
data, ee = get_output_from_qchem(qc_input,
                                 processors=1,
                                 force_recalculation=False,
                                 return_electronic_structure=True,
                                 )

print(data)


def compute_JK(density_matrix, Vee):

    n_basis_functions = len(density_matrix)
    J = np.zeros((n_basis_functions, n_basis_functions))
    K = np.zeros((n_basis_functions, n_basis_functions))

    def J_element(i, j):
        J_sum = 0.0
        K_sum = 0.0

        for k in range(n_basis_functions):
            for l in range(n_basis_functions):
                density = density_matrix[k, l]
                J = Vee[i, j, k, l]
                K = Vee[i, l, k, j]
                J_sum += density * J
                K_sum += density * K

        return J_sum/2, K_sum/2  # handle double counting

    for i in range(n_basis_functions):
        for j in range(n_basis_functions):
            J[i, j], K[i, j] = J_element(i, j)

    return J, K


# extract data from Q-Chem
density_matrix = np.array(ee['total_scf_density'])
Vee = np.array(ee['ao_integrals'])

# get J & K matrices
J, K = compute_JK(density_matrix, Vee)

# compute energies
coulomb_energy = np.sum(density_matrix * J)
exchange_energy = -0.5 * np.sum(density_matrix * K)

print('Computed Properties\n')
print('Total Coulomb:', coulomb_energy)
print('HF Exchange:', exchange_energy)

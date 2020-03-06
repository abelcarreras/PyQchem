from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.structure import Structure
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
print('Final energy:', parsed_data['energy'])

energies = [step['energy'] for step in  parsed_data['optimization_steps']]
gradients = [step['gradient'] for step in parsed_data['optimization_steps']]
displacements = [step['displacement'] for step in parsed_data['optimization_steps']]

# Plot data
fig, axs = plt.subplots(3, sharex=True, gridspec_kw={'hspace': 0})
fig.suptitle('Optimization steps')
axs[0].set(ylabel='Energy')
axs[0].plot(energies, color='b')
axs[1].set(ylabel='Gradient')
axs[1].plot(gradients, color='r')
axs[2].set(ylabel='Displacement', xlabel='Steps')
axs[2].plot(displacements, color='g')

plt.show()
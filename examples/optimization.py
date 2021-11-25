# Geomtry optimization of water molecule
from pyqchem import get_output_from_qchem, Structure, QchemInput
from pyqchem.parsers.parser_optimization import basic_optimization
import matplotlib.pyplot as plt
from pyqchem.errors import OutputError

# define molecule
coordinates = [[0.0,  1.0,  0.0],
               [0.0,  0.0,  1.0],
               [0.0,  0.0, -1.0]]

symbols = ['O', 'H', 'H']

molecule = Structure(coordinates=coordinates,
                     symbols=symbols,
                     charge=0,
                     multiplicity=1)

print('Initial structure')
print(molecule)

# Optimization constrains definition
constrains = {'stre': [{'atoms': [1, 2], 'value': 0.9},
                       {'atoms': [1, 3], 'value': 0.9}],
              'bend': [{'atoms': [2, 1, 3], 'value': 110.0}]}

qc_input = QchemInput(molecule,
                      jobtype='opt',
                      exchange='hf',
                      basis='sto-3g',
                      geom_opt_tol_gradient=300,
                      geom_opt_tol_energy=100,
                      geom_opt_tol_displacement=1200,
                      geom_opt_max_cycles=50,  # reduce this number to test not convergence case
                      geom_opt_constrains=constrains  # comment/uncomment this line to use or not constrains
                      )

try:
    parsed_data = get_output_from_qchem(qc_input,
                                        processors=4,
                                        parser=basic_optimization,
                                        force_recalculation=False)

    opt_molecule = parsed_data['optimized_molecule']

    print('Optimized structure')
    print(opt_molecule)
    print('Final energy:', parsed_data['energy'])

except OutputError as e:
    print(e.error_lines)
    parsed_data = basic_optimization(e.full_output)
    print('Last structure')
    print(parsed_data['optimization_steps'][-1]['molecule'])
    print('Last energy:', parsed_data['optimization_steps'][-1]['energy'])


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
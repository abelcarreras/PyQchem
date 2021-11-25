# Computation of transition state optimization and IRC
from pyqchem import get_output_from_qchem, QchemInput
from pyqchem.parsers.parser_frequencies import basic_frequencies
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.parsers.parser_irc import basic_irc
from pyqchem.structure import Structure
from pyqchem.file_io import write_structure_to_xyz
import matplotlib.pyplot as plt


# define molecule
coordinates = [[ 0.6268743917,   -0.1366254266,   -0.0000000000],
               [ 0.3711605704,    1.0377672206,    0.0000000000],
               [-0.5903438458,   -0.0311449516,    0.0000000000]]

symbols = ['C', 'H', 'N']

molecule = Structure(coordinates=coordinates,
                     symbols=symbols,
                     charge=0,
                     multiplicity=1)

# Transition state optimization
qc_input = QchemInput(molecule,
                      jobtype='ts',
                      exchange='hf',
                      basis='sto-3g',
                      geom_opt_tol_gradient=300,
                      geom_opt_tol_energy=100,
                      geom_opt_tol_displacement=1200,
                      geom_opt_max_cycles=50,  # reduce this number to test not convergence case
                      )

opt_data, ee = get_output_from_qchem(qc_input,
                                     processors=4,
                                     parser=basic_optimization,
                                     force_recalculation=False,
                                     read_fchk=True,
                                     store_full_output=True)


print('Transition state')
print(opt_data['optimized_molecule'])

# frequencies calculation
qc_input = QchemInput(opt_data['optimized_molecule'],
                      jobtype='freq',
                      exchange='hf',
                      basis='sto-3g',
                      sym_ignore=True,
                      symmetry=False,
                      scf_guess=ee['coefficients'])

freq_data = get_output_from_qchem(qc_input,
                                  processors=4,
                                  # force_recalculation=True,
                                  parser=basic_frequencies,
                                  store_full_output=True)

# IRC calculation
qc_input = QchemInput(opt_data['optimized_molecule'],
                      jobtype='rpath',
                      exchange='hf',
                      basis='sto-3g',
                      rpath_max_cycles=30,
                      scf_guess=ee['coefficients'],
                      hessian=freq_data['hessian']
                      )

irc_data = get_output_from_qchem(qc_input,
                                 processors=4,
                                 parser=basic_irc,
                                 store_full_output=True
                                 )

# write coordinates into file
write_structure_to_xyz([step['molecule'] for step in irc_data['irc_forward']][::-1] +
                       [step['molecule'] for step in irc_data['irc_backward']], 'irc.xyz')

# plot SCF energies
energies_f = [step['energy'] for step in irc_data['irc_forward']]
energies_b = [step['energy'] for step in irc_data['irc_backward']]

plt.plot(energies_f[::-1] + energies_b)
plt.xticks([], [])
plt.xlabel('Intrinsic reaction coordinate')
plt.ylabel('Energy (Hartree)')
plt.show()

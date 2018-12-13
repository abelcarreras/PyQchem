from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input
from pyqchem.parsers.parser_frequencies import basic_frequencies
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.structure import Structure


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
print(molecule.get_xyz())

# optimization
txt_input = create_qchem_input(molecule,
                               jobtype='opt',
                               exchange='hf',
                               basis='sto-3g',
                               geom_opt_tol_gradient=300,
                               geom_opt_tol_energy=100,
                               geom_opt_coords=-1,
                               geom_opt_tol_displacement=1200)

parsed_data = get_output_from_qchem(txt_input,
                                    processors=4,
                                    parser=basic_optimization)


opt_molecule = parsed_data['optimized_molecule']

print('Optimized structure')
print(opt_molecule.get_xyz())

# frequencies calculation
txt_input = create_qchem_input(opt_molecule,
                               jobtype='freq',
                               exchange='hf',
                               basis='sto-3g',)

parsed_data = get_output_from_qchem(txt_input,
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




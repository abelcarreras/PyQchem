from pyqchem.structure import Structure
import numpy as np


def basic_optimization(output, print_data=False):

    # Molecule
    n = output.find('$molecule')
    n2 = output[n:].find('$end')

    molecule_region = output[n:n+n2-1].replace('\t', ' ').split('\n')[1:]
    charge, multiplicity = [int(num) for num in molecule_region[0].split()]
    coordinates = np.array([np.array(line.split()[1:4], dtype=float) for line in molecule_region[1:]])
    symbols = [line.split()[0].capitalize() for line in molecule_region[1:]]
    n_atoms = len(coordinates)

    # Optimization
    n = output.find('**  OPTIMIZATION CONVERGED  **')
    ne = output[n-200:n].find('Final energy')

    final_energy = float(output[ne+n-200: n].split()[3])
    optimization_section = output[n:]
    coordinates_section = optimization_section.split('\n')
    coordinates_final = [line.split()[2:5] for line in coordinates_section[5:5+n_atoms]]

    optimized_molecule = Structure(coordinates=np.array(coordinates_final, dtype=float).tolist(),
                                   atomic_elements=symbols,
                                   charge=charge,
                                   multiplicity=multiplicity)

    return {'optimized_molecule': optimized_molecule,
            'energy': final_energy}



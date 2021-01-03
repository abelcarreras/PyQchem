from pyqchem.structure import Structure
import numpy as np
import re


def basic_optimization(output, print_data=False):

    data_dict = {}
    # Molecule
    n = output.find('$molecule')
    n2 = output[n:].find('$end')

    molecule_region = output[n:n+n2-1].replace('\t', ' ').split('\n')[1:]
    charge, multiplicity = [int(num) for num in molecule_region[0].split()]
    coordinates = np.array([np.array(line.split()[1:4], dtype=float) for line in molecule_region[1:]])
    symbols = [line.split()[0].capitalize() for line in molecule_region[1:]]
    n_atoms = len(coordinates)

    step_s2 = None
    # Optimization steps
    optimization_steps = []
    list_iterations = [l.end() for l in re.finditer('Optimization Cycle', output)]
    for ini, fin in zip(list_iterations, list_iterations[1:] + [len(output)]):
        step_section = output[ini:fin]
        enum = step_section.find('Coordinates (Angstroms)')
        atoms_list = step_section[enum:].split('\n')[2:n_atoms+2]
        coordinates_step = np.array([atom.split()[2:] for atom in atoms_list], dtype=float).tolist()

        step_molecule = Structure(coordinates=coordinates_step,
                                  symbols=symbols,
                                  charge=charge,
                                  multiplicity=multiplicity)

        enum = step_section.find('Energy is')
        step_energy = float(step_section[enum: enum+50].split()[2])
        enum = step_section.find('      Gradient')
        step_gradient = float(step_section[enum: enum+50].split()[1])
        enum = step_section.find('      Displacement')
        step_displacement = float(step_section[enum: enum+50].split()[1])

        enum = step_section.find('<S^2>')
        if enum > 0:
            step_s2 = float(step_section[enum: enum+50].split()[2])

        optimization_steps.append({'molecule': step_molecule,
                                   'energy': step_energy,
                                   'gradient': step_gradient,
                                   'displacement': step_displacement})

    data_dict['optimization_steps'] = optimization_steps

    # Optimization Convergence
    enum = output.find('**  OPTIMIZATION CONVERGED  **')
    if enum > 0:
        ne = output[enum-200:enum].find('Final energy')

        final_energy = float(output[ne+enum-200: enum].split()[3])
        optimization_section = output[enum:]
        coordinates_section = optimization_section.split('\n')
        coordinates_final = [line.split()[2:5] for line in coordinates_section[5:5+n_atoms]]

        optimized_molecule = Structure(coordinates=np.array(coordinates_final, dtype=float).tolist(),
                                       symbols=symbols,
                                       charge=charge,
                                       multiplicity=multiplicity)

        data_dict['optimized_molecule'] = optimized_molecule
        data_dict['energy'] = final_energy
        data_dict['s2'] = step_s2

    return data_dict



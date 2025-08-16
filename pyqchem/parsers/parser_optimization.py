from pyqchem.structure import Structure
from pyqchem.parsers.common import read_input_structure
import numpy as np
import re


def parse_molecule(opt_section, charge=0, multiplicity=0):

    num = opt_section.find('Coordinates')
    if num < 0:
        num = opt_section.find('Standard Nuclear Orientation')
        coordinates_section = opt_section[num:].split('\n')[3:]
    else:
        coordinates_section = opt_section[num:].split('\n')[2:]

    # check number of atoms
    n_atoms = 0
    for i in coordinates_section:
        if len(i.split()) != 5:
            break
        n_atoms += 1

    coordinates_section = coordinates_section[:n_atoms]

    symbols = [line.split()[1] for line in coordinates_section]
    coordinates = [line.split()[2:] for line in coordinates_section]

    molecule = Structure(coordinates=np.array(coordinates, dtype=float).tolist(),
                         symbols=symbols,
                         charge=charge,
                         multiplicity=multiplicity)

    return molecule


def basic_optimization(output, print_data=False):

    data_dict = {}

    # Molecule
    structure = read_input_structure(output)

    step_s2 = None
    # Optimization steps
    optimization_steps = []
    list_iterations = [l.end() for l in re.finditer('Optimization Cycle', output, re.IGNORECASE)]
    for ini, fin in zip(list_iterations, list_iterations[1:] + [len(output)]):
        step_section = output[ini:fin]
        enum = step_section.find('Coordinates (Angstroms)')

        step_molecule = parse_molecule(step_section, structure.charge, structure.multiplicity)

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

        optimized_molecule = parse_molecule(optimization_section, structure.charge, structure.multiplicity)

        data_dict['optimized_molecule'] = optimized_molecule
        data_dict['energy'] = final_energy
        data_dict['s2'] = step_s2

    # Optimization Convergence
    enum = output.find('** TRANSITION STATE CONVERGED  **')
    if enum > 0:
        ne = output[enum - 200:enum].find('Final energy')

        final_energy = float(output[ne + enum - 200: enum].split()[3])
        optimization_section = output[enum:]

        optimized_molecule = parse_molecule(optimization_section, structure.charge, structure.multiplicity)

        data_dict['optimized_molecule'] = optimized_molecule
        data_dict['energy'] = final_energy
        data_dict['s2'] = step_s2

    return data_dict



from pyqchem.structure import Structure
import numpy as np
import re


def basic_irc(output, print_data=False):

    branch_mark = True
    data_dict = {}
    # Molecule
    n = output.find('$molecule')
    n2 = output[n:].find('$end')

    molecule_region = output[n:n+n2-1].replace('\t', ' ').split('\n')[1:]
    charge, multiplicity = [int(num) for num in molecule_region[0].split()]
    coordinates = np.array([np.array(line.split()[1:4], dtype=float) for line in molecule_region[1:]])
    symbols = [line.split()[0].capitalize() for line in molecule_region[1:]]
    n_atoms = len(coordinates)

    forward_steps = []
    backward_steps = []

    list_iterations = [l.end() for l in re.finditer('Reaction path following', output)]
    for ini, fin in zip(list_iterations, list_iterations[1:] + [len(output)]):
        step_section = output[ini:fin]
        enum = step_section.find('Standard Nuclear Orientation')
        atoms_list = step_section[enum:].split('\n')[3:n_atoms+3]
        coordinates_step = np.array([atom.split()[2:] for atom in atoms_list], dtype=float).tolist()

        step_energy = None
        for l in re.finditer('Total energy in the final basis set', step_section):
            step_energy = float(step_section[l.end(): l.end()+50].split()[1])

        step_molecule = Structure(coordinates=coordinates_step,
                                  symbols=symbols,
                                  charge=charge,
                                  multiplicity=multiplicity)

        if (step_section.find('IRC -- maximum number of cycles reached') > 0
                or step_section.find('IRC -- convergence criterion reached') > 0):
            if branch_mark:
                #forward_steps = forward_steps[::-1]
                branch_mark = False
            else:
                break

        if branch_mark:
            forward_steps.append({'molecule': step_molecule,
                                  'energy': step_energy,
                                  })
        else:
            backward_steps.append({'molecule': step_molecule,
                                   'energy': step_energy,
                                   })

    data_dict['irc_forward'] = forward_steps
    data_dict['irc_backward'] = forward_steps

    return data_dict
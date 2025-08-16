from pyqchem.structure import Structure
from pyqchem.parsers.common import read_input_structure, read_scf_info
import numpy as np
import re


def basic_irc(output, print_data=False):

    branch_mark = True
    data_dict = {}

    # Molecule
    structure = read_input_structure(output)
    n_atoms = structure.get_number_of_atoms()

    forward_steps = []
    backward_steps = []

    list_iterations = [l.end() for l in re.finditer('Reaction path following', output)]
    for ini, fin in zip(list_iterations, list_iterations[1:] + [len(output)]):
        step_section = output[ini:fin]
        enum = step_section.find('Standard Nuclear Orientation')
        atoms_list = step_section[enum:].split('\n')[3:n_atoms+3]
        coordinates_step = np.array([atom.split()[2:] for atom in atoms_list], dtype=float).tolist()

        if step_section.find('Total energy') < 0:
            break

        step_energy = read_scf_info(step_section)['scf_energy']

        step_molecule = Structure(coordinates=coordinates_step,
                                  symbols=structure.get_symbols(),
                                  charge=structure.charge,
                                  multiplicity=structure.multiplicity)

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
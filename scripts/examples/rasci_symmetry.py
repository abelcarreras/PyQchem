import numpy as np

from pyqchem.symmetry import get_wf_symmetry
from pyqchem.utils import get_plane, set_zero_coefficients
from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.structure import Structure
from pyqchem.file_io import build_fchk
from copy import deepcopy
from pyqchem.parsers.parser_rasci import rasci as rasci_parser
from pyqchem.utils import reorder_coefficients


def get_occupied_electrons(configuration, structure):
    alpha_e = np.sum([int(c) for c in configuration['alpha']])
    beta_e = np.sum([int(c) for c in configuration['beta']])
    return (structure.number_of_electrons + structure.charge - (alpha_e + beta_e))//2


def get_state_symmetry(parsed_fchk,
                       rasci_states,
                       center=None,
                       orientation=(0, 0, 1),
                       group='C2h',
                       extra_print=False
                       ):

    sym_states = {}
    for i, state in enumerate(rasci_states):

        structure = parsed_fchk['structure']
        total_orbitals = len(parsed_fchk['coefficients']['alpha'])

        if extra_print:
            print('HF/KS Ground state\n---------------------------')
            molsym = get_wf_symmetry(structure,
                                     parsed_fchk['basis'],
                                     parsed_fchk['coefficients'],
                                     center=center,
                                     orientation=orientation,
                                     group=group)

            molsym.print_alpha_mo_IRD()
            molsym.print_beta_mo_IRD()
            molsym.print_wf_mo_IRD()

            print('\nRASCI Exited state {}\n---------------------------'.format(i+1))

        # print(configurations)
        occupations_list = []
        for configuration in state['configurations']:
            occupied_orbitals = get_occupied_electrons(configuration, structure)
            n_extra = total_orbitals - occupied_orbitals - len(configuration['alpha'])
            vector_alpha = [1] * occupied_orbitals + [int(c) for c in configuration['alpha']] + [0] * n_extra

            n_extra = total_orbitals - occupied_orbitals - len(configuration['beta'])
            vector_beta = [1] * occupied_orbitals + [int(c) for c in configuration['beta']] + [0] * n_extra

            occupations_list.append({'alpha': vector_alpha, 'beta': vector_beta})

            if extra_print:
                print(configuration['alpha'], configuration['beta'], configuration['amplitude'])



        state_symmetry_list = []
        for occupations in occupations_list:
            # print('occupations', occupations)

            reordered_coefficients = reorder_coefficients(occupations, parsed_fchk['coefficients'])

            molsym = get_wf_symmetry(structure,
                                     parsed_fchk['basis'],
                                     reordered_coefficients,
                                     center=center,
                                     orientation=orientation,
                                     group=group)

            if extra_print:
                molsym.print_alpha_mo_IRD()
                molsym.print_beta_mo_IRD()
                molsym.print_wf_mo_IRD()

            state_symmetry = {}
            for label, measure in zip(molsym.IRLab, molsym.wf_IRd):
                state_symmetry[label] = measure

            # return maximum symmetry and value
            if extra_print:
                print([molsym.IRLab[np.argmax(molsym.wf_IRd)], np.max(molsym.wf_IRd)])
            state_symmetry_list.append([molsym.IRLab[np.argmax(molsym.wf_IRd)], np.max(molsym.wf_IRd)])

        # Make sure symmetry of all configurations is the same
        assert len(np.unique([a[0] for a in state_symmetry_list])) == 1

        symmetry_label = state_symmetry_list[0][0]
        average_measure = np.average([a[1] for a in state_symmetry_list])
        sym_states['state {}'.format(i+1)] = [symmetry_label, average_measure]
    return sym_states


ethene = [[0.0,  0.0000,   0.65750],
          [0.0,  0.0000,  -0.65750],
          [0.0,  0.92281,  1.22792],
          [0.0, -0.92281,  1.22792],
          [0.0, -0.92281, -1.22792],
          [0.0,  0.92281, -1.22792]]

symbols = ['C', 'C', 'H', 'H', 'H', 'H']

range_f1 = range(0, 6)

# create molecule
molecule = Structure(coordinates=ethene,
                     atomic_elements=symbols,
                     charge=0,
                     multiplicity=1)

# create Q-Chem input
qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      correlation='rasci',
                      ras_act=6,
                      ras_elec=6,
                      basis='sto-3g',
                      ras_roots=10,
                      ras_do_hole=False,
                      ras_do_part=False)

print(qc_input.get_txt())


# get data from Q-Chem calculation
output, err, electronic_structure = get_output_from_qchem(qc_input,
                                                          processors=4,
                                                          force_recalculation=False,
                                                          read_fchk=True,
                                                          parser=rasci_parser,
                                                          store_full_output=True)


# store original fchk info in file
open('test.fchk', 'w').write(build_fchk(electronic_structure))

# print(electronic_structure['nato_coefficients'])


# get plane from coordinates
coordinates_f1 = np.array(electronic_structure['structure'].get_coordinates())[range_f1]
center_f1, normal_f1 = get_plane(coordinates_f1)


wfsym = get_state_symmetry(electronic_structure,
                           output['excited states rasci'],
                           center=center_f1,
                           orientation=normal_f1,
                           group='D2h',
                           extra_print=False
                           )


for state in wfsym.items():
    print('{:10} '.format(state[0]) + ' {:5} {:5.3f}'.format(*state[1]))


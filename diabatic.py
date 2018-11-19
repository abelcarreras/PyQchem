#!/usr/bin/env python3

from structure import Structure
from qchem_core import get_output_from_qchem, create_qchem_input
import numpy as np
import pickle
import argparse

from parsers.basic import basic_parser_qchem
from parsers.parser_diabatic import analyze_diabatic  # parser for 4 states only


# Argument parser
parser = argparse.ArgumentParser(description='diabatic contributions')
parser.add_argument('filename', metavar='filename', type=str,
                    help='filename for input')
parser.add_argument('-d', metavar='distance', type=float, default=4.7,
                    help='intermolecular distance between monomers')

parser.add_argument('--yrange', metavar='yrange', type=float, nargs=3,
                    default=[0, 3, 0.5],
                    help='range: initial, final, step')

parser.add_argument('--zrange', metavar='zrange', type=float, nargs=3,
                    default=[0, 3, 0.5],
                    help='range: initial, final, step')

parser.add_argument('-n_states', metavar='N', type=int, default=4,
                    help='number of states')


args = parser.parse_args()


# get list of interesting states
def get_frontier_states(cis_data, n_states=4):
    interesting_transitions = []
    for i, state in enumerate(cis_data):
        for transition in state['transitions']:
            if transition['origin'] > ocup_orb - n_states // 2 and transition['target'] <= n_states // 2:
                interesting_transitions.append([i + 1,
                                                transition['origin'],
                                                transition['target'],
                                                transition['amplitude'] ** 2])

    interesting_transitions = np.array(interesting_transitions)

    list_states = interesting_transitions[np.array(interesting_transitions)[:, 3].argsort()]
    list_states = np.array(list_states[:, 0][::-1], dtype=int)
    indexes = np.unique(list_states, return_index=True)[1]
    list_states = np.sort([list_states[index] for index in sorted(indexes)][:n_states])
    return list_states


# monkey patch analyze_diabatic function defaults (change default behavior)
# analyze_diabatic.__setattr__('__defaults__', (False, 0.2, 6))

# common qchem input parameters
parameters = {'jobtype': 'sp',
              'exchange': 'hf',
              'basis': '6-31G',
              # cis
              'cis_n_roots': 20,
              'cis_convergence': 8,
              'cis_singlets': True,
              'cis_triplets': False,
              'cis_ampl_anal': True,
              # other
              'RPA': False,
              'gui': 0}

# grid parameters
distance = args.d
range_y = np.arange(args.yrange[0], args.yrange[1]+args.yrange[2], args.yrange[2]).tolist()
range_z = np.arange(args.zrange[0], args.zrange[1]+args.zrange[2], args.zrange[2]).tolist()

# Start parsing
total_data = {}
for slide_y in range_y:
    for slide_z in range_z:

        coordinates = [[0.0000000,  0.0000000,  0.6660120],
                       [0.0000000,  0.0000000, -0.6660120],
                       [0.0000000,  0.9228100,  1.2279200],
                       [0.0000000, -0.9228100,  1.2279200],
                       [0.0000000, -0.9228100, -1.2279200],
                       [0.0000000,  0.9228100, -1.2279200],
                       [distance,   0.0000000,  0.6660120],
                       [distance,   0.0000000, -0.6660120],
                       [distance,   0.9228100,  1.2279200],
                       [distance,  -0.9228100,  1.2279200],
                       [distance,  -0.9228100, -1.2279200],
                       [distance,   0.9228100, -1.2279200]]

        coordinates = np.array(coordinates)

        coordinates[:, 1] = coordinates[:, 1] + slide_y
        coordinates[:, 2] = coordinates[:, 2] + slide_z

        molecule = Structure(coordinates=coordinates,
                             atomic_elements=['C', 'C', 'H', 'H', 'H', 'H', 'C', 'C', 'H', 'H', 'H', 'H'],
                             charge=0)

        txt_input = create_qchem_input(molecule, **parameters)
        # print(txt_input)

        # parse CIS data
        data = get_output_from_qchem(txt_input, processors=4, force_recalculation=False, parser=basic_parser_qchem)

        # if closed shell
        ocup_orb = (np.sum(molecule.get_atomic_numbers()) - molecule.charge)//2

        # get interesting states
        list_states = get_frontier_states(data['excited states cis'], n_states=args.n_states)

        # store transition moment info
        states_info = [data['excited states cis'][state] for state in list_states]

        # update parameters
        parameters.update({'loc_cis_ov_separate': False,
                           'er_cis_numstate': args.n_states,
                           'cis_diabath_decompose': True,
                           'localized_diabatization': list_states})
        txt_input = create_qchem_input(molecule, **parameters)

        try:
            # parse adiabatic/diabatic data
            data = get_output_from_qchem(txt_input, processors=4, force_recalculation=True,
                                         parser=analyze_diabatic)
            data.update({'states_info': states_info})
            total_data['{}_{}'.format(slide_y, slide_z)] = data
            print(data)
        except:
            print('Failed!')

        print(data)

total_data.update({'range_y': range_y, 'range_z': range_z})

# store data on disk in pickle format
with open(args.filename, 'wb') as f:
   pickle.dump(total_data, f, pickle.HIGHEST_PROTOCOL)

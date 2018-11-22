#!/usr/bin/env python3

from structure import Structure
from qchem_core import get_output_from_qchem, create_qchem_input
import numpy as np
import pickle
import argparse
import molecules

from parsers.basic import basic_parser_qchem
from parsers.parser_diabatic import analyze_diabatic  # parser for 4 states only


# Argument parser
parser = argparse.ArgumentParser(description='diabatic contributions')
parser.add_argument('filename', metavar='filename', type=str,
                    help='filename for input')
parser.add_argument('-d', metavar='distance', type=float, default=4.7,
                    help='intermolecular distance between monomers')

parser.add_argument('-mol', metavar='name', type=str, default='dimer_ethene',
                    help='molecule name (in molecules.py)', )

parser.add_argument('--yrange', metavar='yrange', type=float, nargs=3,
                    default=[0, 3, 0.5],
                    help='range: initial, final, step')

parser.add_argument('--zrange', metavar='zrange', type=float, nargs=3,
                    default=[0, 3, 0.5],
                    help='range: initial, final, step')

parser.add_argument('--drange', metavar='drange', type=float, nargs=3,
                    default=None,
                    help='distance range: initial, final, step')

parser.add_argument('--force_recalculation', action='store_true',
                   help='force reacalculation of previously parsed data points')

parser.add_argument('-n_states', metavar='N', type=int, default=4,
                    help='number of states')

parser.add_argument('--origin', metavar='distance', type=int, default=None, nargs='*',
                    help='intermolecular distance between monomers')

parser.add_argument('--target', metavar='distance', type=int, default=None, nargs='*',
                    help='intermolecular distance between monomers')

args = parser.parse_args()

# get list of interesting states from frontier orbitals
def get_frontier_states(cis_data, ocup_orb, n_states=4):
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


# get list of interesting states from selected orbitals
def get_states_from_orbitals(cis_data, origin_orbitals, target_orbitals, n_states=4):
    interesting_transitions = []
    for i, state in enumerate(cis_data):
        for transition in state['transitions']:
            if transition['origin'] in origin_orbitals and transition['target'] in target_orbitals:
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
distance = [args.d]
range_y = np.arange(args.yrange[0], args.yrange[1]+args.yrange[2], args.yrange[2]).tolist()
range_z = np.arange(args.zrange[0], args.zrange[1]+args.zrange[2], args.zrange[2]).tolist()

# set distance range
if args.drange is not None:
    distance = np.arange(args.drange[0], args.drange[1] + args.drange[2], args.drange[2]).tolist()
    range_y = [0]
    range_z = [0]


# molecule selection
get_molecule = getattr(molecules, args.mol)

# Start calculation
total_data = {}

for slide_d in distance:
    for slide_y in range_y:
        for slide_z in range_z:

            # print('12')

            molecule, parser_info = get_molecule(slide_d, slide_y, slide_z)

            print(molecule.get_number_of_atoms())
            print('dist: {}  y: {}  z: {}'.format(slide_d, slide_y, slide_z))
            for s, c in zip(molecule.get_atomic_elements(), molecule.get_coordinates()):
                print('{:2} '.format(s) + '{:10.8f} {:10.8f} {:10.8f}'.format(*c))

            txt_input = create_qchem_input(molecule, **parameters)
            # print(txt_input)

            # parse CIS data
            data = get_output_from_qchem(txt_input, processors=4, force_recalculation=args.force_recalculation,
                                         parser=basic_parser_qchem)

            # get interesting states
            if args.origin is not None and args.target is not None:
                print('Using defined Orgin-Target orbitals with {} states'.format(args.n_states))
                list_states = get_states_from_orbitals(data['excited states cis'], n_states=args.n_states,
                                                       origin_orbitals=args.origin, target_orbitals=args.target)
            else:
                # assuming closed shell
                ocup_orb = (np.sum(molecule.get_atomic_numbers()) - molecule.charge) // 2
                print('Using {} states from frontier orbitals. Occupied orbitals: {}'.format(args.n_states, ocup_orb))
                list_states = get_frontier_states(data['excited states cis'], ocup_orb, n_states=args.n_states)

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
                data = get_output_from_qchem(txt_input, processors=4, force_recalculation=args.force_recalculation,
                                             parser=analyze_diabatic, parser_parameters=parser_info)
                data.update({'states_info': states_info})
                if args.drange is None:
                    total_data['{}_{}'.format(slide_y, slide_z)] = data
                else:
                    total_data['{}'.format(slide_d)] = data
                print(data)
            except:
                print('Failed!')


total_data.update({'range_y': range_y, 'range_z': range_z, 'distance': distance})

# store data on disk in pickle format
with open(args.filename, 'wb') as f:
   pickle.dump(total_data, f, pickle.HIGHEST_PROTOCOL)

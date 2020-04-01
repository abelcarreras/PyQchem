from pyqchem.structure import Structure
from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input
from pyqchem.parsers.basic import basic_parser_qchem
import numpy as np


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

distance = 4.0

slide_y = 0
slide_z = 0

range_z = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4]
range_y = [0]
for slide_y in range_y:
    for slide_z in range_z:
        print('y:{}  z:{}'.format(slide_y, slide_z))
        coordinates = [[0.0000000, 0.0000000, 0.6660120],
                       [0.0000000, 0.0000000, -0.6660120],
                       [0.0000000, 0.9228100, 1.2279200],
                       [0.0000000, -0.9228100, 1.2279200],
                       [0.0000000, -0.9228100, -1.2279200],
                       [0.0000000, 0.9228100, -1.2279200],
                       [distance, 0.0000000, 0.6660120],
                       [distance, 0.0000000, -0.6660120],
                       [distance, 0.9228100, 1.2279200],
                       [distance, -0.9228100, 1.2279200],
                       [distance, -0.9228100, -1.2279200],
                       [distance, 0.9228100, -1.2279200]]

        coordinates = np.array(coordinates)

        coordinates[6:, 1] = coordinates[6:, 1] + slide_y
        coordinates[6:, 2] = coordinates[6:, 2] + slide_z

        molecule = Structure(coordinates=coordinates,
                             symbols=['C', 'C', 'H', 'H', 'H', 'H', 'C', 'C', 'H', 'H', 'H', 'H'],
                             charge=0)

        txt_input = create_qchem_input(molecule, **parameters)
        # print(txt_input)

        data = get_output_from_qchem(txt_input, processors=4, force_recalculation=False, parser=basic_parser_qchem)

        # if closed shell
        ocup_orb = (np.sum(molecule.get_atomic_numbers()) - molecule.charge) // 2
        n_states = 2

        # get list of interesting states
        interesting_transitions = []
        for i, state in enumerate(data['excited states cis']):
            for transition in state['transitions']:
                if transition['origin'] > ocup_orb - n_states and transition['target'] <= n_states:
                    interesting_transitions.append(
                        [i + 1, transition['origin'], transition['target'], transition['amplitude'] ** 2])

        interesting_transitions = np.array(interesting_transitions)
        list_states = interesting_transitions[np.array(interesting_transitions)[:, 3].argsort()]
        list_states = np.array(list_states[:, 0][::-1], dtype=int)
        indexes = np.unique(list_states, return_index=True)[1]
        list_states = np.sort([list_states[index] for index in sorted(indexes)][:n_states * 2])
        # print(list_states)
        trans_moments = []
        for state in list_states:
            trans_moments.append(data['excited states cis'][state]['transition moment'])
        print(trans_moments)

        parameters.update({'loc_cis_ov_separate': False,
                           'er_cis_numstate': n_states * 2,
                           'cis_diabath_decompose': True,
                           'localized_diabatization': list_states})
        txt_input = create_qchem_input(molecule, **parameters)
        # print(txt_input)
        data = get_output_from_qchem(txt_input, force_recalculation=True, parser=None)
        with open('qchemout_{}_{}.out'.format(slide_y, slide_z), 'w') as f:
            f.write(data)

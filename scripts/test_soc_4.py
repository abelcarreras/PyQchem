from pyqchem.qchem_core import get_output_from_qchem, redefine_calculation_data_filename
from pyqchem.qc_input import QchemInput
from pyqchem.structure import Structure
from pyqchem.parsers.parser_rasci import parser_rasci as parser_rasci
from pyqchem.basis import get_basis_from_ccRepo, trucate_basis, basis_to_txt
import numpy as np
from pyqchem.errors import OutputError


redefine_calculation_data_filename('test_soc3s.pkl')

as_c_triplet = [[3, 1], 4, 1]
as_c_singlet = [[2, 2], 4, 1]
atom_c = Structure(coordinates=[[0.0, 0.0, 0.0]],
                   symbols=['C'],
                   charge=-1,
                   multiplicity=4,
                   name='C')

as_o_triplet = [[4, 2], 4, 1]
as_o_singlet = [[3, 3], 4, 1]
atom_o = Structure(coordinates=[[0.0, 0.0, 0.0]],
                   symbols=['O'],
                   charge=-2,
                   multiplicity=1,
                   name='O')

as_si_triplet = [[7, 5], 8, 1]
as_si_singlet = [[6, 6], 8, 1]
atom_si = Structure(coordinates=[[0.0, 0.0, 0.0]],
                    symbols=['Si'],
                    charge=-1,
                    multiplicity=4,
                    name='Si')

as_s_triplet = [[8, 6], 8, 1]
as_s_singlet = [[7, 7], 8, 1]
atom_s = Structure(coordinates=[[0.0, 0.0, 0.0]],
                   symbols=['S'],
                   charge=-2,
                   multiplicity=1,
                   name='S')


experimental = {'C': [10192.657, 21648.030],
                'O' : [15867.862, 33792.583],
                'Si': [6298.850,  15394.370],
                'S': [9238.609,  22179.954]}

print('\nExperimental values (cm-1)')
for e, v in experimental.items():
    print('{}: {}'.format(e, v))

use_experimental = False

for atom, active_spaces in [(atom_c, {'singlet': as_c_singlet, 'triplet': as_c_triplet}),
                            (atom_o, {'singlet': as_o_singlet, 'triplet': as_o_triplet}),
                            (atom_si, {'singlet': as_si_singlet, 'triplet': as_si_triplet}),
                            (atom_s, {'singlet': as_s_singlet, 'triplet': as_s_triplet})]:
    for basis_name in ['cc-pCVTZ']:
        basis_custom_repo = get_basis_from_ccRepo(atom, basis_name)

        basis_custom_repo = trucate_basis(basis_custom_repo,
                                          shells=['G', 'H', 'I', 'J', 'K'])

        qc_input_s = QchemInput(atom,
                                jobtype='sp',
                                exchange='hf',
                                correlation='rasci',
                                unrestricted='False',
                                thresh=14,
                                scf_convergence=11,
                                cis_convergence=9,
                                max_cis_cycles=150,
                                max_scf_cycles=100,
                                basis=basis_custom_repo,
                                ras_elec_alpha=active_spaces['singlet'][0][0],
                                ras_elec_beta=active_spaces['singlet'][0][1],
                                ras_act=active_spaces['singlet'][1],
                                ras_occ=active_spaces['singlet'][2],
                                ras_spin_mult=0,
                                ras_roots=10,
                                calc_soc=1,
                                n_frozen_core=0,
                                # ras_do_part=False,
                                set_iter=1000,
                                mem_total=15000,
                                mem_static=900
                                )

        qc_input_t = QchemInput(atom,
                                jobtype='sp',
                                exchange='hf',
                                correlation='rasci',
                                unrestricted='False',
                                thresh=14,
                                scf_convergence=11,
                                cis_convergence=9,
                                max_cis_cycles=150,
                                max_scf_cycles=100,
                                basis=basis_custom_repo,
                                ras_elec_alpha=active_spaces['triplet'][0][0],
                                ras_elec_beta=active_spaces['triplet'][0][1],
                                ras_act=active_spaces['triplet'][1],
                                ras_occ=active_spaces['triplet'][2],
                                ras_spin_mult=0,
                                ras_roots=10,
                                calc_soc=1,
                                n_frozen_core=0,
                                set_iter=1000,
                                mem_total=15000,
                                mem_static=900
                                )

        output_s = get_output_from_qchem(qc_input_s,
                                         processors=4,
                                         force_recalculation=False,
                                         # parser=parser_rasci,
                                         store_full_output=True
                                         )
        # print(output_s)
        with open('qchem_output_{}.out'.format(atom.name), 'w') as f:
            f.write(output_s)

        output_s = parser_rasci(output_s)

        output_t = get_output_from_qchem(qc_input_t,
                                         processors=4,
                                         force_recalculation=False,
                                         parser=parser_rasci,
                                         store_full_output=True
                                         )

        print('\n---------------------------------------------')
        print('Atom: {}'.format(atom.name))
        print('Basis: {}'.format(basis_name))
        print('Active space singlet (ele, act, occ): {}'.format(active_spaces['singlet']))
        print('---------------------------------------------')

        for i, state in enumerate(output_s['excited_states']):
            print('Energy state {} ({}):  {: 18.12f} au'.format(i+1, state['multiplicity'], state['total_energy']))
            print(' Alpha  Beta   Amplitude')
            for j, conf in enumerate(state['configurations']):
                print('  {}  {} {:8.3f}'.format(conf['alpha'], conf['beta'], conf['amplitude']))

        print('\nsoc_tot (cm-1) [imaginary part]')
        print('   ' + ''.join(['{:^18}'.format(n) for n in range(1, 9)]))
        for i in range(1, 10):
            line = '{:3}'.format(i)
            for j in range(1, 10):
                try:
                    line += '{:18.12f}'.format(np.array(output_s['interstate_properties'][(i, j)]['total_soc_mat'])[0, 0].imag)
                except KeyError:
                    line += '        -         '
            print(line)

        print('Frozen occupied: {} / virtual: {}'.format(output_s['rasci_dimensions']['frozen_occupied'],
                                                         output_s['rasci_dimensions']['frozen_virtual']))

        print('\n---------------------------------------------')
        print('Atom: {}'.format(atom.name))
        print('Basis: {}'.format(basis_name))
        print('Active space triplet (ele, act, occ): {}'.format(active_spaces['triplet']))
        print('---------------------------------------------')

        for i, state in enumerate(output_t['excited_states']):
            print('Energy state {} ({}):  {: 18.12f} au'.format(i+1, state['multiplicity'], state['total_energy']))
            print(' Alpha  Beta   Amplitude')
            for j, conf in enumerate(state['configurations']):
                print('  {}  {} {:8.3f}'.format(conf['alpha'], conf['beta'], conf['amplitude']))

        print('\nsoc_tot (cm-1) [imaginary part]')
        print('   ' + ''.join(['{:^18}'.format(n) for n in range(1, 9)]))
        for i in range(1, 10):
            line = '{:3}'.format(i)
            for j in range(1, 10):
                try:
                    line += '{:18.12f}'.format(np.array(output_t['interstate_properties'][(i, j)]['total_soc_mat'])[0, 0].imag)
                except KeyError:
                    line += '        -         '
            print(line)

        print('Frozen occupied: {} / virtual: {}'.format(output_t['rasci_dimensions']['frozen_occupied'],
                                                         output_t['rasci_dimensions']['frozen_virtual']))

        np.set_printoptions(precision=6, suppress=True)
        coupling = np.zeros((15, 15), dtype=complex)
        for i in range(3):
            for j in range(3):
                if i != j:
                    for ms_i in range(3):
                        for ms_j in range(3):
                            # print(output_t['interstate_properties'])
                            coupling[i*3+ms_i, j*3+ms_j] = np.array(output_t['interstate_properties'][(i+1, j+1)]['total_soc_mat'])[ms_i, ms_j]

        for i in range(3):
            for j in range(6):
                    for ms_i in range(3):
                        coupling[i*3+ms_i, j+9] = np.array(output_s['interstate_properties'][(i+1, j+4)]['total_soc_mat'])[0, ms_i]
                        coupling[j+9, i*3+ms_i] = np.array(output_s['interstate_properties'][(j+4, i+1)]['total_soc_mat'])[ms_i, 0]

        if use_experimental:
            for i in range(5):
                coupling[i + 9, i + 9] = experimental[atom.name][0]
            coupling[14, 14] = experimental[atom.name][1]

        else:
            for i in range(3):
                for j in range(3):
                    coupling[3*i+j, 3*i+j] = output_s['excited_states'][i]['excitation_energy']*8065.54429 # eV to cm-1

            for i in range(6):
                coupling[i + 9, i + 9] = output_s['excited_states'][i+3]['excitation_energy']*8065.54429 # eV to cm-1

        print('\n Couplings matrix [12 states] (cm-1)')
        t_labels = ' '.join(np.array([['T{:2}(ms={:2})              '.format(j, i) for i in [-1, 0, 1]] for j in ['1', '2', '3']]).flatten())
        s_labels = ' '.join(['S{:2}                     '.format(i) for i in range(4, 10)])

        print('            {}    {}'.format(t_labels, s_labels))

        print(np.array2string(coupling, max_line_width=10000))

        print('\n Diagonalizated couplings matrix [12 states] (cm-1)')

        eval = np.linalg.eigvalsh(coupling)
        # eval, evec = np.linalg.eigh(coupling)

        # for val, u in zip(eval, evec.T):
        #     print(np.round(np.linalg.norm(val * u - np.dot(coupling, u)), decimals=8))

        print(np.array2string(np.diag(eval), max_line_width=10000))

        coupling = coupling[:9, :9]
        print('\n Couplings matrix [9 states] (cm-1)')
        t_labels = ' '.join(np.array([['T{:2}(ms={:2})           '.format(j, i) for i in [-1, 0, 1]] for j in ['1', '2', '3']]).flatten())
        print('           {}'.format(t_labels))

        print(np.array2string(coupling, max_line_width=10000))

        print('\n Diagonalizated couplings matrix [9 states] (cm-1)')

        eval = np.linalg.eigvalsh(coupling)

        print(np.array2string(np.diag(eval), max_line_width=10000))

        print('\n***********************************************************************\n')

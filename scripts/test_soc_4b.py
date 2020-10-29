from pyqchem.qchem_core import get_output_from_qchem, redefine_calculation_data_filename
from pyqchem.qc_input import QchemInput
from pyqchem.structure import Structure
from pyqchem.parsers.parser_rasci import parser_rasci as parser_rasci
from pyqchem.basis import get_basis_from_ccRepo, trucate_basis, basis_to_txt
import numpy as np
from pyqchem.errors import OutputError


redefine_calculation_data_filename('test_soc3s_b.pkl')

as_f_doblet = [[4, 3], 4, 1]
atom_f = Structure(coordinates=[[0.0, 0.0, 0.0]],
                   symbols=['F'],
                   charge=-1,
                   multiplicity=1,
                   name='F')

as_cl_doblet = [[8, 7], 8, 1]
atom_cl = Structure(coordinates=[[0.0, 0.0, 0.0]],
                    symbols=['Cl'],
                    charge=-1,
                    multiplicity=1,
                    name='Cl')


experimental = {'F': [0, 0],
                'Cl': [0,0]}

use_experimental = False

for atom, active_spaces in [(atom_cl, as_cl_doblet), (atom_f, as_f_doblet)]:
    for basis_name in ['cc-pCVTZ']:
        basis_custom_repo = get_basis_from_ccRepo(atom, basis_name)

        basis_custom_repo = trucate_basis(basis_custom_repo,
                                          shells=['G', 'H', 'I', 'J', 'K'])

        qc_input_d = QchemInput(atom,
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
                                ras_elec_alpha=active_spaces[0][0],
                                ras_elec_beta=active_spaces[0][1],
                                ras_act=active_spaces[1],
                                ras_occ=active_spaces[2],
                                ras_spin_mult=0,
                                ras_roots=10,
                                calc_soc=1,
                                n_frozen_core=0,
                                # ras_do_part=False,
                                set_iter=1000,
                                mem_total=15000,
                                mem_static=900
                                )

        output_d = get_output_from_qchem(qc_input_d,
                                         processors=4,
                                         force_recalculation=False,
                                         # parser=parser_rasci,
                                         store_full_output=True
                                         )
        # print(output_s)
        with open('qchem_output_{}.out'.format(atom.name), 'w') as f:
            f.write(output_d)

        output_d = parser_rasci(output_d)

        print('\n---------------------------------------------')
        print('Atom: {}'.format(atom.name))
        print('Basis: {}'.format(basis_name))
        print('Active space (ele, act, occ): {}'.format(active_spaces))
        print('---------------------------------------------')

        for i, state in enumerate(output_d['excited_states']):
            print('Energy state {} ({}):  {: 18.12f} au'.format(i+1, state['multiplicity'], state['total_energy']))
            print(' Alpha  Beta   Amplitude')
            for j, conf in enumerate(state['configurations']):
                print('  {}  {} {:8.3f}'.format(conf['alpha'], conf['beta'], conf['amplitude']))

        print('\nsoc_tot (cm-1) [imaginary part]')
        print('   ' + ''.join(['{:^18}'.format(n) for n in range(1, 11)]))
        for i in range(1, 11):
            line = '{:3}'.format(i)
            for j in range(1, 11):
                try:
                    line += '{:18.12f}'.format(np.array(output_d['interstate_properties'][(i, j)]['total_soc_mat'])[0, 0].imag)
                except KeyError:
                    line += '        -         '
            print(line)

        print('Frozen occupied: {} / virtual: {}'.format(output_d['rasci_dimensions']['frozen_occupied'],
                                                         output_d['rasci_dimensions']['frozen_virtual']))


        np.set_printoptions(precision=6, suppress=True)
        coupling = np.zeros((6, 6), dtype=complex)
        for i in range(3):
            for j in range(3):
                if i != j:
                    for ms_i in range(2):
                        for ms_j in range(2):
                            # print(output_t['interstate_properties'])
                            # print(np.array(output_d['interstate_properties'][(i+1, j+1)]['total_soc_mat']).shape, '  :  ', [ms_i, ms_j])
                            coupling[i*2+ms_i, j*2+ms_j] = np.array(output_d['interstate_properties'][(i+1, j+1)]['total_soc_mat'])[ms_i, ms_j]

        if use_experimental:
            for i in range(5):
                coupling[i + 9, i + 9] = experimental[atom.name][0]
            coupling[14, 14] = experimental[atom.name][1]


        print('\n Couplings matrix (cm-1)')
        t_labels = ' '.join(np.array([['D{:2}(ms={:2})           '.format(j, i) for i in [1/2, -1/2]] for j in ['1', '2', '3']]).flatten())
        print('           {}'.format(t_labels))

        print(np.array2string(coupling, max_line_width=10000))

        print('\n Diagonalizated couplings matrix (cm-1)')

        eval = np.linalg.eigvalsh(coupling)

        print(np.array2string(np.diag(eval), max_line_width=10000))

        print('\n***********************************************************************\n')

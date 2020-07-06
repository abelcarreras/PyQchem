from pyqchem.qchem_core import get_output_from_qchem, redefine_calculation_data_filename
from pyqchem.qc_input import QchemInput
from pyqchem.structure import Structure
from pyqchem.parsers.parser_rasci import parser_rasci as parser_rasci
from pyqchem.basis import get_basis_from_ccRepo, trucate_basis, basis_to_txt
import numpy as np
from pyqchem.errors import OutputError


redefine_calculation_data_filename('test_soc3.pkl')

as_c = [[[3, 1], 4, 1]]
atom_c = Structure(coordinates=[[0.0, 0.0, 0.0]],
                   symbols=['C'],
                   charge=-1,
                   multiplicity=4,
                   name='C')

as_o = [[[4, 2], 4, 1]]
atom_o = Structure(coordinates=[[0.0, 0.0, 0.0]],
                   symbols=['O'],
                   charge=-2,
                   multiplicity=1,
                   name='O')

as_si = [[[3, 1], 4, 5], [[6, 4], 7, 2]]
atom_si = Structure(coordinates=[[0.0, 0.0, 0.0]],
                    symbols=['Si'],
                    charge=-1,
                    multiplicity=4,
                    name='Si')

as_s = [[[3, 1], 4, 5], [[7, 5], 7, 2], [[8, 6], 8, 1]]
atom_s = Structure(coordinates=[[0.0, 0.0, 0.0]],
                   symbols=['S'],
                   charge=-2,
                   multiplicity=1,
                   name='S')

as_f = [[[4, 3], 4, 1]]
atom_f = Structure(coordinates=[[0.0, 0.0, 0.0]],
                   symbols=['F'],
                   charge=-1,
                   multiplicity=1,
                   name='F')

as_cl = [[[4, 3], 4, 5], [[7, 6], 7, 2], [[8, 7], 8, 1]]
atom_cl = Structure(coordinates=[[0.0, 0.0, 0.0]],
                    symbols=['Cl'],
                    charge=-1,
                    multiplicity=1,
                    name='Cl')


for atom, active_space_list in [(atom_s, as_s)]:#, (atom_si, as_si), (atom_c, as_c), (atom_o, as_o),  (atom_f, as_f), (atom_cl, as_cl)]:
    for basis_name in ['cc-pVDZ', 'cc-pCVTZ', 'cc-pCVQZ']:
        for active_space in active_space_list:
            basis_custom_repo = get_basis_from_ccRepo(atom, basis_name)

            basis_custom_repo = trucate_basis(basis_custom_repo,
                                              shells=['G', 'H', 'I', 'J', 'K'])

            qc_input = QchemInput(atom,
                                  jobtype='sp',
                                  exchange='hf',
                                  correlation='rasci',
                                  unrestricted='False',
                                  thresh=14,
                                  scf_convergence=8,
                                  max_cis_cycles=150,
                                  max_scf_cycles=100,
                                  basis=basis_custom_repo,
                                  # ras_elec=active_space[0],
                                  ras_elec_alpha=active_space[0][0],
                                  ras_elec_beta=active_space[0][1],
                                  ras_act=active_space[1],
                                  ras_occ=active_space[2],
                                  ras_spin_mult=0,
                                  ras_roots=8,  # calculate 5 states
                                  calc_soc=1,
                                  # n_frozen_core=0,
                                  set_iter=1000,
                                  mem_total=15000,
                                  mem_static=900
                                  )

            try:
                output = get_output_from_qchem(qc_input,
                                               processors=4,
                                               force_recalculation=False,
                                               parser=parser_rasci,
                                               store_full_output=True
                                               )
            except OutputError as e:
                print('---------------------------------------------')
                print('Atom: {}'.format(atom.name))
                print('basis: {}'.format(basis_name))
                print('Active space ([alpha, beta], act, occ): {}'.format(active_space))
                print(e.error_lines)
                print('calculation_failed')
                continue

            print('\n---------------------------------------------')
            print('Atom: {}'.format(atom.name))
            print('Basis: {}'.format(basis_name))
            print('Active space (ele, act, occ): {}'.format(active_space))
            print('---------------------------------------------')

            for i, state in enumerate(output['excited_states']):
                print('Energy state {} ({}):  {: 18.12f} au'.format(i+1, state['multiplicity'], state['total_energy']))
                print(' Alpha  Beta   Amplitude')
                for j, conf in enumerate(state['configurations']):
                    print('  {}  {} {:8.3f}'.format(conf['alpha'], conf['beta'], conf['amplitude']))

            print('\ngamma_tot')
            print('   ' + ''.join(['{:^18}'.format(n) for n in range(1, 9)]))
            for i in range(1, 9):
                line = '{:3}'.format(i)
                for j in range(1, 9):
                    try:
                        line += '{:18.12f}'.format(output['interstate_properties'][(i, j)]['gamma_total'])
                    except KeyError:
                        line += '        -         '
                print(line)

            print('\nsoc_1e (cm-1) [imaginary part]')
            print('   ' + ''.join(['{:^18}'.format(n) for n in range(1, 9)]))
            for i in range(1, 9):
                line = '{:3}'.format(i)
                for j in range(1, 9):
                    try:
                        line += '{:18.12f}'.format(np.array(output['interstate_properties'][(i, j)]['1e_soc_mat'])[0, 0].imag)
                    except KeyError:
                        line += '        -         '
                print(line)

            print('\nsoc_2e (cm-1) [imaginary part]')
            print('   ' + ''.join(['{:^18}'.format(n) for n in range(1, 9)]))
            for i in range(1, 9):
                line = '{:3}'.format(i)
                for j in range(1, 9):
                    try:
                        line += '{:18.12f}'.format(np.array(output['interstate_properties'][(i, j)]['2e_soc_mat'])[0, 0].imag)
                    except KeyError:
                        line += '        -         '
                print(line)

            print('\nsoc_tot (cm-1) [imaginary part]')
            print('   ' + ''.join(['{:^18}'.format(n) for n in range(1, 9)]))
            for i in range(1, 9):
                line = '{:3}'.format(i)
                for j in range(1, 9):
                    try:
                        line += '{:18.12f}'.format(np.array(output['interstate_properties'][(i, j)]['total_soc_mat'])[0, 0].imag)
                    except KeyError:
                        line += '        -         '
                print(line)

            print('Frozen occupied: {} / virtual: {}'.format(output['rasci_dimensions']['frozen_occupied'],
                                                             output['rasci_dimensions']['frozen_virtual']))

            print('\nSOCC (cm-1)')
            print('   ' + ''.join(['{:^18}'.format(n) for n in range(1, 9)]))
            for i in range(1, 9):
                line = '{:3}'.format(i)
                for j in range(1, 9):
                    try:
                        line += '{:18.12f}'.format(output['interstate_properties'][(i, j)]['mf_socc'])
                    except KeyError:
                        line += '        -         '
                print(line)

            print('Frozen occupied: {} / virtual: {}'.format(output['rasci_dimensions']['frozen_occupied'],
                                                             output['rasci_dimensions']['frozen_virtual']))

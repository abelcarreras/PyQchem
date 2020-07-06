from pyqchem.qchem_core import get_output_from_qchem, redefine_calculation_data_filename
from pyqchem.qc_input import QchemInput
from pyqchem.structure import Structure
from pyqchem.parsers.parser_rasci import parser_rasci as parser_rasci
from pyqchem.basis import get_basis_from_ccRepo, trucate_basis, basis_to_txt
from pyqchem.errors import OutputError
import numpy as np

# import pyqchem.qchem_core
# pyqchem.qchem_core.__calculation_data_filename__ = 'test.pkl'
redefine_calculation_data_filename('test_soc2.pkl')

# create molecule
as_clo = [[3, 2, 11], [5, 3, 10], [11, 6, 7], [13, 7, 6], [19, 10, 3]]
mol_clo = Structure(coordinates=[[0.0, 0.0, 0.0000],
                                 [0.0, 0.0, 1.5696]],
                    symbols=['Cl', 'O'],
                    charge=-1,
                    multiplicity=1,
                    name='ClO')

as_bro = [[3, 2, 20], [9, 5, 17], [13, 7, 15], [23, 12, 10]]
mol_bro = Structure(coordinates=[[0.0, 0.0, 0.0000],
                                 [0.0, 0.0, 1.7172]],
                    symbols=['Br', 'O'],
                    charge=-1,
                    multiplicity=1,
                    name='BrO')

as_ncs = [[3, 2, 13], [7, 4, 11], [11, 6, 9], [15, 8, 7]]
mol_ncs = Structure(coordinates=[[0.0, 0.0, 1.159],
                                 [0.0, 0.0, 0.000],
                                 [0.0, 0.0,-1.631]],
                    symbols=['N', 'C', 'S'],
                    charge=-1,
                    multiplicity=1,
                    name='NCS')

as_n2o = [[3, 2, 9], [11, 6, 5], [15, 8, 3]]
mol_n2o = Structure(coordinates=[[0.0, 0.0, 1.154],
                                 [0.0, 0.0, 0.000],
                                 [0.0, 0.0,-1.182]],
                    symbols=['N', 'N', 'O'],
                    charge=0,
                    multiplicity=1,
                    name='N2O')

as_cco = [[3, 2, 9], [5, 3, 8], [11, 6, 5], [15, 8, 3]]
mol_cco = Structure(coordinates=[[0.0, 0.0, 1.308],
                                 [0.0, 0.0, 0.000],
                                 [0.0, 0.0,-1.221]],
                    symbols=['C', 'C', 'O'],
                    charge=-2,
                    multiplicity=1,
                    name='CCO')

as_ccf = [[3, 2, 9], [5, 3, 8], [11, 6, 5], [15, 8, 3]]
mol_ccf = Structure(coordinates=[[0.0, 0.0, 1.271],
                                 [0.0, 0.0, 0.000],
                                 [0.0, 0.0,-1.276]],
                    symbols=['C', 'C', 'F'],
                    charge=-1,
                    multiplicity=1,
                    name='CCF')

as_cccl = [[3, 2, 13], [5, 3, 12], [11, 6, 9], [15, 8, 7]]
mol_cccl = Structure(coordinates=[[0.0, 0.0, 1.2964],
                                  [0.0, 0.0, 0.000],
                                  [0.0, 0.0,-1.6224]],
                     symbols=['C', 'C', 'Cl'],
                     charge=-1,
                     multiplicity=1,
                     name='CCCl')


# Generate basis name list
basis_name_list = []

for molecule, active_space_list in [(mol_clo, as_clo), (mol_bro, as_bro), (mol_ncs, as_ncs), (mol_n2o, as_n2o),
                                    (mol_cco, as_cco), (mol_ccf, as_ccf), (mol_cccl, as_cccl)]:
    for basis_name in ['cc-pCVTZ', 'cc-pCVQZ']:
        for active_space in active_space_list:
            basis_custom_repo = get_basis_from_ccRepo(molecule,
                                                      basis_name,
                                                      if_missing=basis_name.replace('cc-pC', 'cc-p'))
            basis_custom_repo = trucate_basis(basis_custom_repo,
                                              shells=['G', 'H', 'I', 'J', 'K'])

            qc_input = QchemInput(molecule,
                                  jobtype='sp',
                                  exchange='hf',
                                  correlation='rasci',
                                  purecart='1111',
                                  unrestricted='False',
                                  thresh=14,
                                  scf_convergence=8,
                                  # max_cis_cycles=150,
                                  basis=basis_custom_repo,
                                  ras_elec=active_space[0],
                                  ras_act=active_space[1],
                                  ras_occ=active_space[2],
                                  ras_spin_mult=0,
                                  ras_roots=2,  # calculate 2 states
                                  calc_soc=1,
                                  set_iter=1000,
                                  mem_total=15000,
                                  mem_static=900,
                                  n_frozen_core=0
                                  )

            # print(qc_input.get_txt())
            try:
                output = get_output_from_qchem(qc_input,
                                               processors=4,
                                               force_recalculation=False,
                                               parser=parser_rasci,
                                               store_full_output=True
                                               )
            except OutputError as e:
                print('---------------------------------------------')
                print('Molecule: {}'.format(molecule.name))
                print('basis: {}'.format(basis_name))
                print('Active space (ele, act, occ): {}'.format(active_space))
                print(e.error_lines)
                print('calculation_failed')
                continue
            # print(_)
            # print(output['interstate_properties'])

            gamma_total = output['interstate_properties'][(1, 2)]['gamma_total']
            soc_1e = np.array(output['interstate_properties'][(1, 2)]['1e_soc_mat'])
            soc_2e = np.array(output['interstate_properties'][(1, 2)]['2e_soc_mat'])
            soc_tot = np.array(output['interstate_properties'][(1, 2)]['total_soc_mat'])
            socc = output['interstate_properties'][(1, 2)]['mf_socc']

            energy_1 = output['excited_states'][0]['total_energy']
            energy_2 = output['excited_states'][1]['total_energy']

            print('---------------------------------------------')
            print('molecule: {}'.format(molecule.name))
            print('basis: {}'.format(basis_name))
            print('Active space (ele, act, occ): {}'.format(active_space))
            print('gamma_tot: {}'.format(gamma_total))
            print('Energy state1:  {: 18.12f} au'.format(energy_1))
            print('Energy state2:  {: 18.12f} au'.format(energy_2))
            print('soc_1e  {0.real: 10.3f} + {0.imag: 10.8f} cm-1'.format(soc_1e[0, 0]))
            print('soc_2e  {0.real: 10.3f} + {0.imag: 10.8f} cm-1'.format(soc_2e[0, 0]))
            print('soc_tot {0.real: 10.3f} + {0.imag: 10.8f} cm-1'.format(soc_tot[0, 0]))
            # print('soc_1e (cm-1)')
            # print('{0.real: 8.3e} + {0.imag: 10.8e}i    {1.real: 8.3e} + {1.imag: 10.8e}i'.format(*soc_1e[0]))
            # print('{0.real: 8.3e} + {0.imag: 10.8e}i    {1.real: 8.3e} + {1.imag: 10.8e}i'.format(*soc_1e[1]))
            # print('soc_2e (cm-1)')
            # print('{0.real: 8.3e} + {0.imag: 10.8e}i    {1.real: 8.3e} + {1.imag: 10.8e}i'.format(*soc_2e[0]))
            # print('{0.real: 8.3e} + {0.imag: 10.8e}i    {1.real: 8.3e} + {1.imag: 10.8e}i'.format(*soc_2e[1]))
            # print('soc_tot (cm-1)')
            # print('{0.real: 8.3e} + {0.imag: 10.8e}i    {1.real: 8.3e} + {1.imag: 10.8e}i'.format(*soc_tot[0]))
            # print('{0.real: 8.3e} + {0.imag: 10.8e}i    {1.real: 8.3e} + {1.imag: 10.8e}i'.format(*soc_tot[1]))

            print('SOCC: {: 18.12f} cm-1'.format(socc))
            print('Frozen occupied: {} / virtual: {}'.format(output['rasci_dimensions']['frozen_occupied'],
                                                             output['rasci_dimensions']['frozen_virtual']))

            print('Frozen occupied: {} / virtual: {}'.format(output['rasci_dimensions']['frozen_occupied'],
                                                             output['rasci_dimensions']['frozen_virtual']))

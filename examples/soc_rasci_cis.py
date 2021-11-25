# Spin-orbit coupling calculation using TD-DFT and RASCI
from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.structure import Structure
from pyqchem.parsers.parser_rasci import parser_rasci
from pyqchem.parsers.parser_cis import basic_cis as parser_cis
import numpy as np

# Carbon atom defintion
atom_c = Structure(coordinates=[[0.0, 0.0, 0.0]],
                   symbols=['C'],
                   charge=-2,
                   multiplicity=1,
                   name='C')

basis_name = '6-31G**'
# RAS_CI calculation
qc_input_rasci = QchemInput(atom_c,
                            jobtype='sp',
                            exchange='hf',
                            correlation='rasci',
                            unrestricted=True,
                            thresh=14,
                            scf_convergence=8,
                            max_cis_cycles=150,
                            max_scf_cycles=100,
                            basis=basis_name,
                            ras_elec_alpha=3,
                            ras_elec_beta=1,
                            ras_act=4,
                            ras_occ=1,
                            ras_spin_mult=0,
                            ras_roots=8,
                            calc_soc=1,
                            n_frozen_core=0,
                            )

# CIS calculation
qc_input_cis = QchemInput(atom_c,
                          jobtype='sp',
                          exchange='b3lyp',
                          basis=basis_name,
                          cis_n_roots=8,
                          cis_singlets=True,
                          cis_triplets=True,
                          cis_convergence=8,
                          calc_soc=1,
                          n_frozen_core=0,
                          )

for qc_input, parser in [(qc_input_cis, parser_cis), (qc_input_rasci, parser_rasci)]:

    output = get_output_from_qchem(qc_input,
                                   processors=4,
                                   force_recalculation=True,
                                   parser=parser,
                                   )

    # print calculation data
    for i, state in enumerate(output['excited_states']):
        print('Energy state {} ({}):  {:18.12f} au'.format(i+1, state['multiplicity'], state['total_energy']))
        print('Transition energy  {:4.8} eV'.format(state['excitation_energy']))
        if qc_input is qc_input_cis:
            print(' Origin  Target   Amplitude')
        else:
            print(' Alpha  Beta   Amplitude')

        for j, conf in enumerate(state['configurations']):
            try:
                print('     {}     {}  {:8.3f}'.format(conf['origin'], conf['target'], conf['amplitude']))
            except KeyError:
                print('  {}  {} {:8.3f}'.format(conf['alpha'], conf['beta'], conf['amplitude']))

    print('\n1e_soc_mat (cm-1) [(0,0) element, imaginary part]')
    print('   ' + ''.join(['{:^18}'.format(n) for n in range(1, 9)]))
    for i in range(1, 9):
        line = '{:3}'.format(i)
        for j in range(1, 9):
            try:
                line += '{:18.12f}'.format(np.array(output['interstate_properties'][(i, j)]['1e_soc_mat'])[0, 0].imag)
            except KeyError:
                line += '        -         '
        print(line)


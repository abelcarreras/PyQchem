# Spin-orbit coupling calculation using TD-DFT and RASCI
from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.structure import Structure
from pyqchem.parsers.parser_rasci import parser_rasci


molecule = Structure(coordinates=[[  1.2452100000,    0.8399100000,    0.0000000000],
                                  [ -0.0000000000,    0.0000000000,    0.0000000000],
                                  [ -1.2452100000,    0.8399100000,   -0.0000000000],
                                  [ -0.0000000000,   -0.6635500000,    0.8647500000],
                                  [  0.0000000000,   -0.6635500000,   -0.8647500000],
                                  [  1.6893340004,    1.1394772932,    0.9278761732],
                                  [  1.6893332483,    1.1394767859,   -0.9278810935],
                                  [ -1.1704700000,    1.9087200000,    0.0000000000],
                                  [ -2.2082000000,    0.3702300000,    0.0000000000]],
                     symbols=['C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'],
                     charge=0,
                     multiplicity=3)


basis_name = 'cc-pVDZ'
qc_input_rasci = QchemInput(molecule,
                            jobtype='sp',
                            exchange='hf',
                            correlation='rasci',
                            unrestricted=False,
                            # thresh=14,
                            scf_convergence=8,
                            max_cis_cycles=100,
                            max_scf_cycles=100,
                            basis=basis_name,
                            ras_elec=6,
                            ras_act=5,
                            ras_occ=9,
                            ras_spin_mult=0,
                            ras_roots=4,
                            calc_soc=1,
                            n_frozen_core=0,
                            state_analysis=True,
                            set_iter=300,
                            )

output = get_output_from_qchem(qc_input_rasci,
                               processors=4,
                               force_recalculation=False,
                               parser=parser_rasci,
                               store_full_output=True,
                               )

# print calculation data
for i, state in enumerate(output['excited_states']):
    print('Energy state {} ({}):  {:18.12f} au'.format(i+1, state['multiplicity'], state['total_energy']))
    print('Transition energy  {:4.8} eV'.format(state['excitation_energy']))
    print(' Alpha  Beta   Amplitude')
    for j, conf in enumerate(state['configurations']):
       print('  {}  {} {:8.3f}'.format(conf['alpha'], conf['beta'], conf['amplitude']))

n_states = len(output['excited_states'])

print('\n1-elec SOCC (cm-1) [states x states]')
print('   ' + ''.join(['{:^18}'.format(n) for n in range(1, n_states+1)]))
for i in range(1, n_states+1):
    line = '{:3}'.format(i)
    for j in range(1, n_states+1):
        try:
            line += '{:18.12f}'.format(output['interstate_properties'][(i, j)]['1e_socc'])
        except KeyError:
            line += '        -         '
    print(line)


print('\nmean-Field SOCC (cm-1) [states x states]')
print('   ' + ''.join(['{:^18}'.format(n) for n in range(1, n_states+1)]))
for i in range(1, n_states+1):
    line = '{:3}'.format(i)
    for j in range(1, n_states+1):
        try:
            line += '{:18.12f}'.format(output['interstate_properties'][(i, j)]['mf_socc'])
        except KeyError:
            line += '        -         '
    print(line)


print('\nAngular Momentum Lx (imaginary component) [states x states]')
print('   ' + ''.join(['{:^18}'.format(n) for n in range(1, n_states+1)]))
for i in range(1, n_states+1):
    line = '{:3}'.format(i)
    for j in range(1, n_states+1):
        try:
            line += '{:18.12f}'.format(output['interstate_properties'][(i, j)]['angular_momentum'][0].imag)
        except KeyError:
            line += '        -         '
    print(line)

print('\nAngular Momentum Ly (imaginary component) [states x states]')
print('   ' + ''.join(['{:^18}'.format(n) for n in range(1, n_states+1)]))
for i in range(1, n_states+1):
    line = '{:3}'.format(i)
    for j in range(1, n_states+1):
        try:
            line += '{:18.12f}'.format(output['interstate_properties'][(i, j)]['angular_momentum'][1].imag)
        except KeyError:
            line += '        -         '
    print(line)

print('\nAngular Momentum Lz (imaginary component) [states x states]')
print('   ' + ''.join(['{:^18}'.format(n) for n in range(1, n_states+1)]))
for i in range(1, n_states+1):
    line = '{:3}'.format(i)
    for j in range(1, n_states+1):
        try:
            line += '{:18.12f}'.format(output['interstate_properties'][(i, j)]['angular_momentum'][2].imag)
        except KeyError:
            line += '        -         '
    print(line)

print('\nSpin matrix Sx (0,0 real component) [states x states]')
print('   ' + ''.join(['{:^18}'.format(n) for n in range(1, n_states+1)]))
for i in range(1, n_states+1):
    line = '{:3}'.format(i)
    for j in range(1, n_states+1):
        try:
            line += '{:18.12f}'.format(output['interstate_properties'][(i, j)]['spin_matrices'][0] [0][0].real)
        except KeyError:
            line += '        -         '
    print(line)

print('\nSpin matrix Sy (0,0 imag component) [states x states]')
print('   ' + ''.join(['{:^18}'.format(n) for n in range(1, n_states+1)]))
for i in range(1, n_states+1):
    line = '{:3}'.format(i)
    for j in range(1, n_states+1):
        try:
            line += '{:18.12f}'.format(output['interstate_properties'][(i, j)]['spin_matrices'][1] [0][0].imag)
        except KeyError:
            line += '        -         '
    print(line)

print('\nSpin matrix Sz (0,0 real component) [states x states]')
print('   ' + ''.join(['{:^18}'.format(n) for n in range(1, n_states+1)]))
for i in range(1, n_states+1):
    line = '{:3}'.format(i)
    for j in range(1, n_states+1):
        try:
            line += '{:18.12f}'.format(output['interstate_properties'][(i, j)]['spin_matrices'][2] [0][0].real)
        except KeyError:
            line += '        -         '
    print(line)


from structure import Structure
from qchem_core import get_output_from_qchem, create_qchem_input
from parsers.basic import basic_parser_qchem
import numpy as np
import matplotlib.pyplot as plt


energies = []
energies_es = []

distances = np.arange(0.2, 5.0, 0.2).tolist() + [10]
# distances = np.arange(0.2, 2.2, 0.2).tolist()
# distances = [2, 2.2, 2.5, 2.6, 2.8]

x_tot = []
for dist in distances:
    molecule = Structure(coordinates=[[0.0, 0.0, 0.0],
                                      [0.0, 0.0, dist]],
                         atomic_elements=['H', 'H'],
                         charge=0)

    txt_input = create_qchem_input(molecule,
                                   jobtype='sp',
                                   exchange='hf',
                                   correlation='rasci',
                                   basis='cc-pVTZ',
                                   #basis='6-31G(d,p)',
                                   ras_act=2,
                                   ras_elec=2,
                                   ras_spin_mult=0,
                                   ras_roots=1,
                                   ras_do_hole=True,
                                   ras_sts_tm=True,
                                   # rasci sr-dft
                                   ras_srdft=True,
                                   ras_omega=400,
                                   ras_srdft_cor='srpw92',
                                   ras_srdft_exc='srpbe',
                                   ras_natorb=False,
                                   ras_print=0,
                                   ras_srdft_spinpol=True,
                                   set_iter=30,
                                   ras_srdft_damp=0.5,
                                   gui=0)

    data = get_output_from_qchem(txt_input, processors=4, parser=basic_parser_qchem)
    #print(output)
    #exit()
    #data = parse_qchem_output(output)
    print(data)
    exit()
    if len([state['total energy'] for state in data['excited states']]) > 0:
        energies.append(data['scf energy'])
        energies_es.append([state['total energy'] for state in data['excited states']])
        x_tot.append(dist)
        print(energies_es[-1])
    else:
        print(dist)

x_tot = np.array(x_tot)
energies_es = np.array(energies_es)

print(energies_es)
print(energies)
plt.plot(x_tot, energies, label='HF')
plt.plot(x_tot, energies_es, label='RASCI-pol')
# np.savetxt('HF_ref.dat', np.vstack((x_tot, energies)).T, fmt='%10.8f')
#np.savetxt('HF_triplet_pol.dat', np.vstack((x_tot, energies_es.T)).T, fmt='%10.8f')

plt.legend()
plt.show()

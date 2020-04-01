from pyqchem.structure import Structure
from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input
from pyqchem.parsers.basic import basic_parser_qchem
import numpy as np
import matplotlib.pyplot as plt


parameters = {'jobtype': 'sp',
              'exchange': 'hf',
              'correlation': 'rasci',
              'basis': '6-31G(d,p)',
              'ras_act': 2,
              'ras_elec': 2,
              'ras_spin_mult': 3,
              'ras_roots': 1,
              'ras_do_hole': True,
              'ras_sts_tm': True,
              # rasci sr-dft
              'ras_srdft': True,
              'ras_srdft_cor': 'srpw92',
              'ras_srdft_exc': 'srpbe',
              'ras_natorb': False,
              'ras_print': 0,
              'ras_srdft_spinpol': True,
              'set_iter': 60,
              'ras_srdft_damp': 0.5,
              'gui': 0}

energies = []
energies_es = []

atoms = ['C', 'O', 'Si']

omegas = [pow(10, i) for i in np.arange(1, 2, 0.25)]

atom_data = {}
for atom in atoms:
    data_omega = []

    for omega in omegas:
        try:
            molecule = Structure(coordinates=[[0.0, 0.0, 0.0]],
                                 symbols=[atom],
                                 charge=0)

            parameters.update({'ras_omega': int(omega)})

            gap_list = []
            for ispin in [True, False]:
                print(ispin)
                parameters.update({'ras_srdft_spinpol': ispin})

                # singlet
                parameters.update({'ras_spin_mult': 1})
                txt_input = create_qchem_input(molecule, **parameters)
                data_singlet = get_output_from_qchem(txt_input, processors=4, parser=basic_parser_qchem)

                # triplet
                parameters.update({'ras_spin_mult': 3})
                txt_input = create_qchem_input(molecule, **parameters)
                data_triplet = get_output_from_qchem(txt_input, processors=4, parser=basic_parser_qchem)

                gap = data_triplet['excited states'][0]['total energy'] - data_singlet['excited states'][0]['total energy']
                gap_list.append(gap)

        except:
            continue

        data_omega.append(gap_list[0]-gap_list[1])

    print(atom, data_omega)
    plt.plot(omegas, data_omega, label=atom)
    plt.xscale('log')
    plt.show()
    atom_data[atom] = data_omega

#plt.legend()
#plt.show()
print(atom_data)
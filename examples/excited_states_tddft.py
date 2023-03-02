# Compute the excited states using RASCI method
from pyqchem import get_output_from_qchem, QchemInput, Structure
from pyqchem.parsers.parser_cis import basic_cis


molecule = Structure(coordinates=[[ -1.07076839,   -2.13175980,    0.03234382],
                                  [ -0.53741536,   -3.05918866,    0.04995793],
                                  [ -2.14073783,   -2.12969357,    0.04016267],
                                  [ -0.39112115,   -0.95974916,    0.00012984],
                                  [  0.67884827,   -0.96181542,   -0.00769025],
                                  [ -1.15875076,    0.37505495,   -0.02522296],
                                  [ -0.62213437,    1.30041753,   -0.05065831],
                                  [ -2.51391203,    0.37767199,   -0.01531698],
                                  [ -3.04726506,    1.30510083,   -0.03293196],
                                  [ -3.05052841,   -0.54769055,    0.01011971]],
                     symbols=['C', 'H', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'H'])


# create Q-Chem input
qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='b3lyp',
                      basis='sto-3g',
                      # unrestricted=True,
                      cis_n_roots=3,
                      cis_singlets=True,
                      cis_triplets=True,
                      )

data, ee = get_output_from_qchem(qc_input,
                                 processors=4,
                                 return_electronic_structure=True,
                                 parser=basic_cis
                                 )


structure = data['structure']
n_occ = structure.number_of_electrons/2


def get_orbital_label(orbital_number):
    """
    get label HOMO-n/LUMO+m labels from orbital number

    :param orbital_number: orbital number (int)
    :return: orbital label (str)
    """
    r_state = int(orbital_number - n_occ)
    if r_state > 0:
        label = 'LUMO +{}'.format(r_state - 1) if r_state - 1 > 0 else 'LUMO'
    else:
        label = 'HOMO {}'.format(r_state) if r_state < 0 else 'HOMO'
    return label


for i_state, state in enumerate(data['excited_states']):
    print('\nState {}'.format(i_state+ 1))
    print('Excitation energy: {} {}'.format(state['excitation_energy'], state['excitation_energy_units']))
    print('Multiplicity: {}'.format(state['multiplicity']))
    print('TDM {} ({:^8.4f})'.format(state['transition_moment'], state['strength']))

    for configuration in state['configurations']:

        print('{:6} {:8} -> {:8} {:8.4f}'.format(configuration['spin'],
                                                 get_orbital_label(configuration['origin']),
                                                 get_orbital_label(configuration['target']),
                                                 configuration['amplitude']))

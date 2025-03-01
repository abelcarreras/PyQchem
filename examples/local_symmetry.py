# Application of local symmetry analysis of RASCI WF for ethylene dimer
# This example requires posym (https://github.com/abelcarreras/posym)
from pyqchem import Structure, QchemInput, get_output_from_qchem
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.parsers.parser_rasci import parser_rasci
from pyqchem.symmetry import get_state_symmetry
from pyqchem.file_io import write_to_fchk
from pyqchem.utils import get_plane, crop_electronic_structure
from pyqchem.symmetry import get_orbitals_symmetry, get_wf_symmetry
import numpy as np


# define monomer
coor_monomer = [[ 0.6695,  0.0000,  0.0000],
                [-0.6695,  0.0000,  0.0000],
                [ 1.2321,  0.9289,  0.0000],
                [ 1.2321, -0.9289,  0.0000],
                [-1.2321,  0.9289,  0.0000],
                [-1.2321, -0.9289,  0.0000]]

symbols_monomer = ['C', 'C', 'H', 'H', 'H', 'H']

monomer = Structure(coordinates=coor_monomer,
                    symbols=symbols_monomer,
                    charge=0,
                    multiplicity=1)

# optimization qchem input
qc_input = QchemInput(monomer,
                      jobtype='opt',
                      exchange='hf',
                      basis='sto-3g',
                      geom_opt_tol_gradient=300,
                      geom_opt_tol_energy=100,
                      geom_opt_coords=-1,
                      geom_opt_tol_displacement=1200)

parsed_data = get_output_from_qchem(qc_input,
                                    processors=4,
                                    parser=basic_optimization)

opt_monomer = parsed_data['optimized_molecule']

print('Optimized monomer structure')
print(opt_monomer)

# Build dimer from monomer
coor_monomer2 = np.array(opt_monomer.get_coordinates())
coor_monomer2[:, 2] += 4.0  # molecules separation

coordinates = opt_monomer.get_coordinates() + coor_monomer2.tolist()
symbols_dimer = symbols_monomer * 2

dimer = Structure(coordinates=coordinates,
                  symbols=symbols_dimer,
                  charge=0,
                  multiplicity=1)

print('Dimer structure')
print(dimer)

# RASCI qchem input
qc_input = QchemInput(dimer,
                      jobtype='sp',
                      exchange='hf',
                      correlation='rasci',
                      basis='sto-3g',
                      ras_act=6,
                      ras_elec=4,
                      ras_spin_mult=1,  # singlets only
                      ras_roots=8,      # calculate 8 states
                      ras_do_hole=False,
                      ras_do_part=False,
                      set_iter=30)

parsed_data, ee = get_output_from_qchem(qc_input,
                                        processors=4,
                                        force_recalculation=False,
                                        parser=parser_rasci,
                                        return_electronic_structure=True,
                                        )

write_to_fchk(ee, filename='dimer.fchk')


sym = get_wf_symmetry(ee,
                      alpha_electrons=5,
                      beta_electrons=6)

states_symmetry = get_state_symmetry(ee,
                                     parsed_data['excited_states'],
                                     group='D2h',
                                     )


# Analysis of diabatic states to use in diabatization
print('\nSymmetry of the electronic states of the dimer')
list_diabatic = []
for i, state in enumerate(parsed_data['excited_states']):
    energy = state['excitation_energy']
    print('State {}: {}'.format(i+1, states_symmetry[i]))


range_frag_1 = list(range(0, 6))
range_frag_2 = list(range(6, 12))

orbitals_sym_dim = get_orbitals_symmetry(ee['structure'],
                                         ee['basis'],
                                         ee['coefficients']['alpha'],
                                         # orbital_numbers=[1, 2, 3]
                                         )

orbitals_sym_monomers = []
for i, range_frag in enumerate([range_frag_1, range_frag_2]):
    coordinates_frag = ee['structure'].get_coordinates(fragment=range_frag)
    center_frag, normal_frag = get_plane(coordinates_frag)

    electronic_structure_frag = crop_electronic_structure(ee, range_frag, renormalize=True)

    # save test fchk file with new coefficients
    write_to_fchk(electronic_structure_frag, filename='fragment_{}.fchk'.format(i))

    orbitals_sym = get_orbitals_symmetry(electronic_structure_frag['structure'],
                                         electronic_structure_frag['basis'],
                                         electronic_structure_frag['coefficients']['alpha'],
                                         center=center_frag,
                                         # orbital_numbers=[1, 2, 3]
                                         )

    orbitals_sym_monomers.append(orbitals_sym)

# Analysis of orbitals symmetry
print('\nSymmetry of the orbitals')

print(' '*15 + '{:5} {:10} {:10}'.format('dimer', 'monomer 1', 'monomer 2'))
for i in range(len(orbitals_sym_monomers[0])):
    # print(i+1, orbitals_sym_monomers[0], orbitals_sym_monomers[1])
    print('Orbital {:3}:     {}     {}        {}'.format(i+1,
                                                orbitals_sym_dim[i],
                                                orbitals_sym_monomers[0][i],
                                                orbitals_sym_monomers[1][i]))


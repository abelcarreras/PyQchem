# Application of local symmetry analysis of RASCI WF for ethylene dimer
# This example requires posym (https://github.com/abelcarreras/posym)
from pyqchem import Structure, QchemInput, get_output_from_qchem
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.parsers.parser_rasci import parser_rasci
from pyqchem.symmetry import get_state_symmetry
from pyqchem.file_io import write_to_fchk
from pyqchem.utils import get_plane, crop_electronic_structure
from pyqchem.symmetry import get_wf_symmetry
from posym import SymmetryBase
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


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

symmetry_measures = get_state_symmetry(ee,
                                       parsed_data['excited_states'],
                                       group='D2h',
                                       orientation=(1, 0, 0),
                                       orientation2=(0, 1, 0),
                                       )

# Analysis of diabatic states to use in diabatization
print('\nSymmetry of the electronic states of the dimer')
list_diabatic = []
for i, state in enumerate(parsed_data['excited_states']):
    sym_lab = symmetry_measures[i][0]
    energy = state['excitation_energy']
    print('State {}: {:3}'.format(i+1, sym_lab))

range_frag = list(range(0, 6))

coordinates_frag = ee['structure'].get_coordinates(fragment=range_frag)
center_frag, normal_frag = get_plane(coordinates_frag)

electronic_structure_frag = crop_electronic_structure(ee, range_frag, renormalize=True)

# save test fchk file with new coefficients
write_to_fchk(electronic_structure_frag, filename='fragment.fchk')

molsym = get_wf_symmetry(electronic_structure_frag['structure'],
                         electronic_structure_frag['basis'],
                         electronic_structure_frag['coefficients'],
                         center=center_frag,
                         orientation=(1, 0, 0),
                         orientation2=(0, 1, 0),
                         group='D2h')

molsym.print_alpha_mo_IRD()
molsym.print_beta_mo_IRD()
ir_labels = molsym.IRLab


def get_symmetry_wf(occupation_alpha, occupation_beta):
    state_wf = SymmetryBase(group='D2h', rep='Ag')
    for orbital_a, occupation in zip(molsym.mo_SOEVs_a, occupation_alpha):
        if abs(occupation) > 0.1:
            state_orb = SymmetryBase(group='D2h',
                                     rep=pd.Series(np.array(orbital_a),
                                                   index=[ "E", "C2", "C2'", "C2''", "sh", "i", "sv", "sd"]))
            state_wf = state_wf * state_orb
    for orbital_b, occupation in zip(molsym.mo_SOEVs_a, occupation_beta):
        if abs(occupation) > 0.1:
            state_orb = SymmetryBase(group='D2h',
                                     rep=pd.Series(np.array(orbital_b),
                                                   index=[ "E", "C2", "C2'", "C2''", "sh", "i", "sv", "sd"]))
            state_wf = state_wf * state_orb

    return state_wf


total_contributions = {}
for istate, state in enumerate(parsed_data['excited_states']):
    for conf in state['configurations']:
        amplitude = conf['amplitude']
        state_wf = get_symmetry_wf(conf['occupations']['alpha'], conf['occupations']['beta'])
        rep = state_wf.get_ir_representation().values
        ir_labels = state_wf.get_ir_representation().keys()

        if 'State: {}'.format(istate+1) in total_contributions:
            total_contributions['State: {}'.format(istate+1)] += amplitude**2 * rep
        else:
            total_contributions['State: {}'.format(istate+1)] = amplitude**2 * rep


print('\nLocal symmetry of electronic states on the monomer')
print('           ' + ' '.join(['{:5}'.format(l) for l in ir_labels]))
for istate, contributions in enumerate(total_contributions.items()):
    print(contributions[0], ' '.join(['{:5.2f}'.format(abs(c)) for c in contributions[1]]))

plot_data = pd.DataFrame(total_contributions, index=ir_labels)
plot_data.plot(kind="bar", width=1)
plt.title("Local symmetry of electronic states on the monomer")
plt.xlabel("Irreducible representations")
plt.ylabel("Contribution")
plt.show()


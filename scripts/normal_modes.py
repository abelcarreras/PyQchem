#!/usr/bin/env python3

from pyqchem.qchem_core import get_output_from_qchem, create_qchem_input
from pyqchem.parsers.parser_frequencies import basic_frequencies
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.structure import Structure

import numpy as np

import argparse


def get_rmsd(coor, coor2):
    coor = np.array(coor)
    coor2 = np.array(coor2)
    rms = np.average(np.square(coor - coor2))
    rmsd = np.sqrt(rms)
    return rmsd


# Argument parser
parser = argparse.ArgumentParser(description='normal modes')
parser.add_argument('filename', metavar='filename', type=str,
                    help='filename for output')
parser.add_argument('--show_plot', action='store_true',
                    help='perform test and show plots')
parser.add_argument('-t', metavar='temperature', type=float, default=300,
                    help='temperature to calculate normal mode amplitude', )


args = parser.parse_args()


coordinates = [[ 0.         ,  2.5593592,   -1.2291437],
               [ 0.         ,  1.2177072,   -0.7388171],
               [ 0.         ,  1.2539645,    0.6577044],
               [ 0.         ,  2.5457806,    1.0474186],
               [ 0.         ,  3.3773954,   -0.1013709],
               [ 0.         ,  0.       ,    1.5890149],
               [ 0.         , -1.2539645,    0.6577044],
               [ 0.         , -1.2177072,   -0.7388171],
               [ 0.         , -2.5593592,   -1.2291437],
               [ 0.         , -3.3773954,   -0.1013709],
               [ 0.         , -2.5457806,    1.0474186],
               [ 0.         , -2.9786512,   -2.6674481],
               [ 0.         ,  0.       ,   -1.4142026],
               [ 0.         , -2.9494503,    2.4837241],
               [ 0.         ,  2.9494503,    2.4837241],
               [ 0.         ,  2.9786512,   -2.6674481],
               [ 1.1494068  ,  0.       ,    2.3950298],
               [ -1.1494068 ,  0.       ,    2.3950298],
               [ 0.         ,  0.       ,   -2.5152032],
               [ 0.         ,  4.4782605,   -0.0840819],
               [ 0.         , -4.4782605,   -0.0840819],
               [ 0.         , -4.0761113,   -2.7556685],
               [ -0.8867106 , -2.6006386,   -3.199612],
               [ 0.8867106  , -2.6006386,   -3.199612],
               [ 0.8941856  , -2.5562337,    2.9920521],
               [ -0.8941856 , -2.5562337,    2.9920521],
               [ 0.         , -4.0461503,    2.5809354],
               [ -0.8867106 ,  2.6006386,   -3.199612],
               [ 0.         ,  4.0761113,   -2.7556685],
               [ 0.8867106  ,  2.6006386,   -3.199612],
               [ 0.         ,  4.0461503,    2.5809354],
               [ -0.8941856 ,  2.5562337,    2.9920521],
               [ 0.8941856  ,  2.5562337,    2.9920521]]

symbols = ['C', 'C', 'N', 'C', 'C', 'B', 'N', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
           'C', 'F', 'F', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
           'H', 'H', 'H']


molecule = Structure(coordinates=coordinates,
                     symbols=symbols,
                     charge=0,
                     multiplicity=1)

molecule.get_coordinates()

print('optimize molecule')

txt_input = create_qchem_input(molecule,
                               jobtype='opt',
                               method='b3lyp',
                               basis='cc-pVDZ',
                               geom_opt_tol_gradient=50,
                               geom_opt_coords=-1,
                               geom_opt_tol_displacement=800)


parsed_data = get_output_from_qchem(txt_input,
                                    processors=4,
                                    use_mpi=True,
                                    parser=basic_optimization,
                                    force_recalculation=False)


#parsed_data = basic_optimization(data)
opt_molecule = parsed_data['optimized_molecule']

print(opt_molecule)

print('calculate frequencies')
txt_input = create_qchem_input(opt_molecule,
                               jobtype='freq',
                               method='b3lyp',
                               basis='cc-pVDZ',)

parsed_data = get_output_from_qchem(txt_input,
                                    processors=4,
                                    use_mpi=True,
                                    force_recalculation=False,
                                    parser=basic_frequencies)

energy_i = parsed_data['energy']

f = open('{}.xyz'.format(args.filename), 'w')
f2 = open(args.filename, 'w')

# print data info

f2.write('OPTIMIZED STRUCTURE\n')
f2.write(opt_molecule)
f2.write('\n\n')
f2.write('DISPLACED STRUCTURES\n')

num_freq = len(parsed_data['frequencies'])

print('frequencies:', parsed_data['frequencies'])

for mode in range(num_freq):

    f2.write('----------------------------------\n')
    freq = parsed_data['frequencies'][mode]
    red_mass = parsed_data['red_mass'][mode]
    force_constants = parsed_data['force_constants'][mode]

    f2.write('mode:                      {}\n'.format(mode+1))
    f2.write('frequency (cm-1):          {}\n'.format(freq))
    f2.write('force constant (mdyne/A):  {}\n'.format(force_constants))

    if freq < 10:
        continue

    # prediction average displacement
    kb = 0.001987204  # kcal/molK
    max_a = np.sqrt(args.t * kb / (force_constants * 143.9325))
    f2.write('Average Kinetic at {:4} K: {:8.3f} Kcal/mol\n'.format(args.t, 0.5 * args.t * kb))

    data_for_plot = []
    for d in np.arange(-max_a, max_a+(max_a/11), max_a/10):
        mode_coor = parsed_data['coordinates'][mode]
        disp_coor = np.array(opt_molecule.get_coordinates()) + d * np.array(mode_coor)
        rmsd = get_rmsd(disp_coor, opt_molecule.get_coordinates())
        f2.write(''.format(rmsd))
        mol = Structure(coordinates=disp_coor,
                        symbols=symbols)

        hartree_to_kcalmol = 627.509474
        # energy2 = 0.5 * freq**2 * red_mass * d**2 * 8.516358951e-5  # kcal/mol
        energy2 = 0.5 * force_constants * d ** 2 * 143.9325  # kcal/mol

        if args.show_plot:
            txt_input = create_qchem_input(mol,
                                           jobtype='freq',
                                           method='b3lyp',
                                           basis='cc-pVDZ')

            parsed_data = get_output_from_qchem(txt_input, processors=4,
                                                force_recalculation=False,
                                                parser=basic_frequencies)
            energy = parsed_data['energy']

            data_for_plot.append([d, (energy - energy_i) * hartree_to_kcalmol, energy2])  # kcal/mol

        f.write(mol.get_xyz(title='{:8.5f} kcal/mol  rmsd: {:10.8f}\n'.format(energy2, rmsd)))
        f2.write(mol.get_xyz(title='{:8.5f} kcal/mol  rmsd: {:10.8f}\n'.format(energy2, rmsd)))

    if args.show_plot:
        import matplotlib.pyplot as plt
        plt.plot(np.array(data_for_plot)[:, 0], np.array(data_for_plot)[:, 1:])
        plt.show()
f.close()
f2.close()

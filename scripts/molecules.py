from pyqchem.structure import Structure
import numpy as np


# Ethene parallel position
def dimer_ethene(distance, slide_y, slide_z):

    coordinates = [[0.0000000,  0.0000000,  0.6660120],
                   [0.0000000,  0.0000000, -0.6660120],
                   [0.0000000,  0.9228100,  1.2279200],
                   [0.0000000, -0.9228100,  1.2279200],
                   [0.0000000, -0.9228100, -1.2279200],
                   [0.0000000,  0.9228100, -1.2279200],
                   [distance,   0.0000000,  0.6660120],
                   [distance,   0.0000000, -0.6660120],
                   [distance,   0.9228100,  1.2279200],
                   [distance,  -0.9228100,  1.2279200],
                   [distance,  -0.9228100, -1.2279200],
                   [distance,   0.9228100, -1.2279200]]

    coordinates = np.array(coordinates)

    coordinates[6:, 1] = coordinates[6:, 1] + slide_y
    coordinates[6:, 2] = coordinates[6:, 2] + slide_z
    symbols = ['C', 'C', 'H', 'H', 'H', 'H', 'C', 'C', 'H', 'H', 'H', 'H']

    molecule = Structure(coordinates=coordinates,
                         symbols=symbols,
                         charge=0)

    return molecule, {'state_threshold': 0.2,
                      'n_mon': 6}


# Tetracloroethene
def dimer_tetrafluoroethene(distance, slide_y, slide_z):

    monomer = [[ 0.6624670117,  0.0000000000, 0.0000000000],
               [-0.6624670117,  0.0000000000, 0.0000000000],
               [ 1.3834661472,  1.0993897934, 0.0000000000],
               [ 1.3834661472, -1.0993897934, 0.0000000000],
               [-1.3834661472, -1.0993897934, 0.0000000000],
               [-1.3834661472,  1.0993897934, 0.0000000000]]

    symbols = ['C', 'C', 'F', 'F', 'F', 'F']

    monomer2 = np.array(monomer)
    #monomer2 = np.dot(monomer, rotation_matrix([0, 1, 0], np.pi / 2))
    monomer2[:, 2] = monomer2[:, 2] + distance
    monomer2[:, 1] = monomer2[:, 1] + slide_y
    monomer2[:, 0] = monomer2[:, 0] + slide_z

    coordinates = np.vstack([monomer, monomer2])

    molecule = Structure(coordinates=coordinates,
                         symbols=symbols * 2,
                         charge=0)

    return molecule, {'state_threshold': 0.2,
                      'n_mon': len(monomer)}


# Tetracloroethene
def dimer_mix(distance, slide_y, slide_z):

    monomer1 = [[ 0.6660120,  0.0000000,  0.0000000,],
                [-0.6660120,  0.0000000,  0.0000000,],
                [ 1.2279200,  0.9228100,  0.0000000,],
                [ 1.2279200, -0.9228100,  0.0000000,],
                [-1.2279200, -0.9228100,  0.0000000,],
                [-1.2279200,  0.9228100,  0.0000000,]]

    symbols1 = ['C', 'C', 'H', 'H', 'H', 'H']

    monomer2 = [[ 0.6624670117,  0.0000000000, 0.0000000000],
                [-0.6624670117,  0.0000000000, 0.0000000000],
                [ 1.3834661472,  1.0993897934, 0.0000000000],
                [ 1.3834661472, -1.0993897934, 0.0000000000],
                [-1.3834661472, -1.0993897934, 0.0000000000],
                [-1.3834661472,  1.0993897934, 0.0000000000]]

    symbols2 = ['C', 'C', 'F', 'F', 'F', 'F']

    monomer2 = np.array(monomer2)

    monomer2[:, 2] = monomer2[:, 2] + distance
    monomer2[:, 1] = monomer2[:, 1] + slide_y
    monomer2[:, 0] = monomer2[:, 0] + slide_z

    coordinates = np.vstack([monomer1, monomer2])
    symbols = symbols1 + symbols2

    molecule = Structure(coordinates=coordinates,
                         symbols=symbols,
                         charge=0)

    return molecule, {'state_threshold': 0.4,
                      'n_mon': len(monomer1)}
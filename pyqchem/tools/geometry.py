import numpy as np


def unit_vector(vector):
    """
    Compute unit vector from general vector

    :param vector: the vector
    :return: the unit vector
    """
    return np.array(vector) / np.linalg.norm(vector)


def get_dihedral(coordinates, atoms):
    """
    Compute the dihedral angle between 4 atoms

    :param coordinates: list of coordinates of the molecule
    :param atoms: list of 4 atom indices to use to calculate the dihedral angle (from 1 to N)
    :return: the dihedral angle
    """

    b1 = np.array(coordinates[atoms[1]-1]) - np.array(coordinates[atoms[0]-1])
    b2 = np.array(coordinates[atoms[2]-1]) - np.array(coordinates[atoms[1]-1])
    b3 = np.array(coordinates[atoms[3]-1]) - np.array(coordinates[atoms[2]-1])

    n1 = unit_vector(np.cross(b1, b2))
    n2 = unit_vector(np.cross(b2, b3))

    b2 = unit_vector(b2)
    m = np.cross(n1, b2)

    assert np.abs(np.dot(b2, n2)) < 1e-8

    # sign given by convention https://en.wikipedia.org/wiki/Dihedral_angle
    return -np.rad2deg(np.arctan2(np.dot(m, n2), np.dot(n1, n2)))


def get_angle(coordinates, atoms):
    """
    Compute the angle between 3 atoms

    :param coordinates: list of coordinates of the molecule
    :param atoms: list of 3 atom indices to use to calculate the angle (from 1 to N)
    :return: the angle
    """

    b1 = np.array(coordinates[atoms[0]-1]) - np.array(coordinates[atoms[1]-1])
    b2 = np.array(coordinates[atoms[2]-1]) - np.array(coordinates[atoms[1]-1])

    b1 = unit_vector(b1)
    b2 = unit_vector(b2)

    dot = np.dot(b1, b2)

    return np.rad2deg(np.arccos(dot))


def get_distance(coordinates, atoms):
    """
    Compute the distance between 2 atoms

    :param coordinates: list of coordinates of the molecule
    :param atoms: list of 2 atom indices to use to calculate the distance (from 1 to N)
    :return: the distance
    """

    b1 = np.array(coordinates[atoms[0]-1]) - np.array(coordinates[atoms[1]-1])

    return np.linalg.norm(b1)


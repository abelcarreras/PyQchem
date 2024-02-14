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
    :return: the dihedral angle in degrees
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
    :return: the angle in degrees
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


def rotate_coordinates(coordinates, angle, axis, atoms_list=None, center=(0, 0, 0)):
    """
    Rotate the coordinates (or range of coordinates) with respect a given axis

    :param coordinates: coordinates to rotate
    :param angle: rotation angle in radians
    :param axis: rotation axis
    :param atoms_list: list of atoms to rotate (if None then rotate all)
    :return: rotated coordinates
    """

    axis = np.array(axis) / np.linalg.norm(axis)  # normalize axis
    coordinates = np.array(coordinates) - np.array(center)

    cos_term = 1 - np.cos(angle)
    rot_matrix = [[axis[0]**2*cos_term + np.cos(angle),              axis[0]*axis[1]*cos_term - axis[2]*np.sin(angle), axis[0]*axis[2]*cos_term + axis[1]*np.sin(angle)],
                  [axis[1]*axis[0]*cos_term + axis[2]*np.sin(angle), axis[1]**2*cos_term + np.cos(angle),              axis[1]*axis[2]*cos_term - axis[0]*np.sin(angle)],
                  [axis[2]*axis[0]*cos_term - axis[1]*np.sin(angle), axis[1]*axis[2]*cos_term + axis[0]*np.sin(angle), axis[2]**2*cos_term + np.cos(angle)]]

    if atoms_list is not None:
        coordinates[atoms_list] = np.dot(coordinates[atoms_list], rot_matrix)
    else:
        coordinates = np.dot(coordinates, rot_matrix) + np.array(center)

    return coordinates.tolist()

import numpy as np


def gaussian(x, s, m):
    """
    1-D Gaussian function

    :param x: coordinate
    :param s: standard deviation
    :param m: center of gaussian
    :return:
    """
    return 1/np.sqrt(2*np.pi*s**2)*np.exp(-0.5*((x-m)/s)**2)

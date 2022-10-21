import numpy as np
from pyqchem.tools.gaussian import gaussian
from pyqchem.units import AU_TO_EV, KB_EV


def get_fcwd(transitions, temperature=300):
    """
    Get Franck-Condon Weighted Density (FCWD) at a given temperature from a vibrational transitions list
    :param transitions: list of vibrational transitions
    :param temperature: temperature
    :return: FCWD
    """

    fcwd = 0
    for t1 in transitions:
        for t2 in transitions:
            reorganization = (t1.reorganization_energy + t2.reorganization_energy)/2
            sigma = np.sqrt(2 * KB_EV * temperature * reorganization)  # Marcus model for band amplitude
            intensity = t1.get_intensity_absorption(temperature) * t2.get_intensity_emission(temperature)
            fcwd += intensity * gaussian(t1.energy_absorption, sigma*np.sqrt(2), t2.energy_emission)
    return fcwd

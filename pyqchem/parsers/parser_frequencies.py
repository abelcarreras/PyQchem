import re
import numpy as np


# parser for frequencies calculations
def basic_frequencies(output, print_data=False):
    """
    Parser for frequencies calculations

    :param output:
    :param print_data:
    :return:
    """

    # Coordinates
    n = output.find('$molecule')
    n2 = output[n:].find('$end')
    molecule_region = output[n:n+n2-1].replace('\t', ' ').split('\n')[1:]
    coordinates = np.array([ np.array(line.split()[1:4], dtype=float) for line in molecule_region[1:]])
    symbols = [line.split()[0].capitalize() for line in molecule_region[1:]]
    n_atoms = len(coordinates)

    # Energy
    n = output.find('Total energy in the final basis set =')
    energy = float(output[n:n+70].split()[8])

    n = output.find('VIBRATIONAL ANALYSIS')
    vibration_section = output[n:]

    modes = []
    frequencies = []
    force_constants = []
    red_mass = []
    ir_active = []
    ir_intens = []
    raman_active = []
    for m in re.finditer('Mode:', vibration_section):
        end_line = vibration_section[m.end():].find('\n')
        modes += vibration_section[m.end():m.end()+end_line].split()[:3]

    for m in re.finditer('Frequency:', vibration_section):
        end_line = vibration_section[m.end():].find('\n')
        frequencies += vibration_section[m.end():m.end()+end_line].split()[:3]

    for m in re.finditer('Force Cnst:', vibration_section):
        end_line = vibration_section[m.end():].find('\n')
        force_constants += vibration_section[m.end():m.end()+end_line].split()[:3]

    for m in re.finditer('Red. Mass:', vibration_section):
        end_line = vibration_section[m.end():].find('\n')
        red_mass += vibration_section[m.end():m.end()+end_line].split()[:3]

    for m in re.finditer('IR Active:', vibration_section):
        end_line = vibration_section[m.end():].find('\n')
        ir_active += vibration_section[m.end():m.end()+end_line].split()[:3]

    for m in re.finditer('IR Intens:', vibration_section):
        end_line = vibration_section[m.end():].find('\n')
        ir_intens += vibration_section[m.end():m.end()+end_line].split()[:3]

    for m in re.finditer('Raman Active:', vibration_section):
        end_line = vibration_section[m.end():].find('\n')
        raman_active += vibration_section[m.end():m.end()+end_line].split()[:3]

    modes = [int(n) for n in modes]
    frequencies = [float(n) for n in frequencies]
    force_constants = [float(n) for n in force_constants]
    red_mass = [float(n) for n in red_mass]
    ir_active = [bool(n) for n in ir_active]
    ir_intens = [float(n) for n in ir_intens]
    raman_active = [bool(n) for n in raman_active]

    nm_coordinates = []
    for i, line in enumerate(vibration_section.split('\n')):
        if 'X      Y      Z' in line:
            nm_coordinate = []
            for j in range(n_atoms):
                coor_lines = vibration_section.split('\n')[j+ i+ 1]
                nm_coordinate.append(coor_lines.split()[1:])

            nm_coordinate = np.array(nm_coordinate, dtype=float)#.reshape(n_atoms, -1)
            # print(nm_coordinate.shape[1], nm_coordinate.shape[1]//3)

            nm_coordinates += [nm_coordinate[:, i*3:(i+1)*3].tolist() for i in range(nm_coordinate.shape[1]//3)]

    return {'modes': modes,
            'frequencies': frequencies,
            'force_constants': force_constants,
            'red_mass': red_mass,
            'ir_active': ir_active,
            'ir_intens': ir_intens,
            'raman_active': raman_active,
            'coordinates': nm_coordinates,
            'energy': energy}

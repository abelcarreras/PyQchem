import re
import numpy as np
from pyqchem.structure import Structure


# parser for frequencies calculations
def basic_frequencies(output, print_data=False):
    """
    Parser for frequencies calculations

    :param output:
    :param print_data:
    :return:
    """

    # Input moelcular data
    n = output.find('$molecule')
    n2 = output[n:].find('$end')
    molecule_region = output[n:n+n2-1].replace('\t', ' ').split('\n')[1:]
    charge, multiplicity = [int(num) for num in molecule_region[0].split()]
    coordinates_input = np.array([ np.array(line.split()[1:4], dtype=float) for line in molecule_region[1:]])
    symbols = [line.split()[0].capitalize() for line in molecule_region[1:]]
    n_atoms = len(coordinates_input)

    # Standard oriented structure
    enum = output.find('Standard Nuclear Orientation')
    section_structure = output[enum:enum + 200*n_atoms].split('\n')
    section_structure = section_structure[3:n_atoms+3]
    coordinates = [[float(num) for num in s.split()[2:]] for s in section_structure]

    structure = Structure(coordinates=coordinates,
                          symbols=symbols,
                          charge=charge,
                          multiplicity=multiplicity)

    # Energy
    n = output.find('Total energy in the final basis set =')
    energy = float(output[n:n+70].split()[8])

    n_hess = output.find('Hessian of the SCF Energy')
    n_van = output.find('VIBRATIONAL ANALYSIS')

    if n_hess > -1:
        # Hessian
        ncol = 6
        ndim = n_atoms * 3
        hessian_section = output[n_hess: n_van]
        hess_block = hessian_section.split('\n')[1:]


        hessian = []
        for i in range(ndim):
            line = []
            for block in range((ndim-1)//ncol + 1):
                line += hess_block[block*(ndim+1) + i +1].split()[1:]
            hessian.append(line)

        hessian = np.array(hessian, dtype=float).tolist()
    else:
        hessian = None

    # Vibration analysis
    vibration_section = output[n_van:]

    frequencies = []
    force_constants = []
    red_mass = []
    ir_active = []
    ir_intens = []
    raman_active = []

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

    frequencies = [float(n) for n in frequencies]
    force_constants = [float(n) for n in force_constants]
    red_mass = [float(n) for n in red_mass]
    ir_active = [bool(n) for n in ir_active]
    ir_intens = [float(n) for n in ir_intens]
    raman_active = [bool(n) for n in raman_active]

    displacements = []
    for i, line in enumerate(vibration_section.split('\n')):
        if 'X      Y      Z' in line:
            disp_coordinate = []
            for j in range(n_atoms):
                coor_lines = vibration_section.split('\n')[j+ i+ 1]
                disp_coordinate.append(coor_lines.split()[1:])

            disp_coordinate = np.array(disp_coordinate, dtype=float)#.reshape(n_atoms, -1)
            # print(nm_coordinate.shape[1], nm_coordinate.shape[1]//3)

            displacements += [disp_coordinate[:, i*3:(i+1)*3].tolist() for i in range(disp_coordinate.shape[1]//3)]

    modes = []
    for i in range(len(frequencies)):
        modes.append({'frequency': frequencies[i],
                      'frequency_units': 'cm-1',
                      'force_constant': force_constants[i],
                      'force_constant_units': 'mDyn/Angs',
                      'reduced_mass': red_mass[i],
                      'reduced_mass_units': 'AMU',
                      'ir_active': ir_active[i],
                      'ir_intensity': ir_intens[i],
                      'ir_intensity_units': 'KM/mol',
                      'raman_active': raman_active[i],
                      'displacement': displacements[i]})

    # thermochemistry
    def get_data_from_line(text):
        n_words = len(text.split())
        n = output.find(text)

        value, units = output[n:n + 70].split()[n_words: n_words + 2]

        return float(value), units

    zpe, zpe_units = get_data_from_line('Zero point vibrational energy')
    trans_enthalpy, trans_enthalpy_units = get_data_from_line('Translational Enthalpy')
    rot_enthalpy, rot_enthalpy_units = get_data_from_line('Rotational Enthalpy')
    vib_enthalpy, vib_enthalpy_units = get_data_from_line('Vibrational Enthalpy')
    rt, rt_units = get_data_from_line('gas constant (RT)')

    trans_entropy, trans_entropy_units = get_data_from_line('Translational Entropy')
    rot_entropy, rot_entropy_units = get_data_from_line('Rotational Entropy')
    vib_entropy, vib_entropy_units = get_data_from_line('Vibrational Entropy')

    tot_enthalpy, tot_enthalpy_unis = get_data_from_line('Total Enthalpy')
    tot_entropy, tot_entropy_units = get_data_from_line('Total Entropy')

    thermochemistry = {'conditions': {'temperature': 298.15, 'temperature_units': 'K',
                                      'pressure': 1.00, 'pressure_units': 'atm' },
                       'zpe': zpe, 'zpe_units': zpe_units,
                       'trans_enthalpy': trans_enthalpy, 'trans_enthalpy_units': trans_enthalpy_units,
                       'rot_enthalpy': rot_enthalpy, 'rot_enthalpy_units': rot_enthalpy_units,
                       'vib_enthalpy': vib_enthalpy, 'vib_enthalpy_units': vib_enthalpy_units,
                       'rt': rt, 'rt_units': rt_units,
                       'trans_entropy': trans_entropy, 'trans_entropy_units': trans_entropy_units,
                       'rot_entropy': rot_entropy, 'rot_entropy_units': rot_entropy_units,
                       'vib_entropy': vib_entropy, 'vib_entropy_units': vib_entropy_units,
                       'tot_enthalpy': tot_enthalpy, 'tot_enthalpy_unis': tot_enthalpy_unis,
                       'tot_entropy': tot_entropy, 'tot_entropy_units': tot_entropy_units,
                       }

    return {'modes': modes,
            'hessian': hessian,
            'thermochemistry': thermochemistry,
            'structure': structure,
            'scf_energy': energy}

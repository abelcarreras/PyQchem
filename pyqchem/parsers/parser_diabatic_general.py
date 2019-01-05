import re
import numpy as np


# parser for 4 states
def analyze_diabatic(output, print_data=False, state_threshold=0.2, n_mon=6):

    # Coordinates
    n = output.find('$molecule')
    n2 = output[n:].find('$end')
    molecule_region = output[n:n+n2-1].replace('\t', ' ').split('\n')[1:]
    coordinates = np.array([np.array(line.split()[1:4], dtype=float) for line in molecule_region[1:]])
    symbols = [line.split()[0].capitalize() for line in molecule_region[1:]]
    n_atoms = len(coordinates)

    # Diabatic states
    for i, line in enumerate(output.split('\n')):
        if 'adiabatic states' in line.lower():
            loc_diab = [int(num) for num in output.split('\n')[i+1].split()]

    if print_data:
        print('-------------------------------')
        print('Mulliken')

    # determine diabatic states type
    state_order = []
    for m in re.finditer('Mulliken analysis of TDA State', output):
        #print output[m.end():m.end()+800].split()
        data = output[m.end():m.end()+37*(n_atoms+3)].split()
        state_number = int(data[0])
        state_info = []
        for i in range(n_atoms):
            # print(data[7+i*4:7+i*4+4][1:4])
            dat = data[7+i*4:7+i*4+4][1:4]
            state_info.append([float(d) for d in dat])

        state_info = np.array(state_info)
        e_sum1 = np.sum(state_info[0:n_mon, 0])
        e_sum2 = np.sum(state_info[n_mon:n_atoms, 0])

        h_sum1 = np.sum(state_info[0:n_mon, 1])
        h_sum2 = np.sum(state_info[n_mon:n_atoms, 1])

        # print (state_info)
        # print ('e_sum', e_sum1, e_sum2)
        # print ('h_sum', h_sum1, h_sum2)

        label = ''
        if state_number in loc_diab:
            if ((np.abs(e_sum1) > 1-state_threshold and np.abs(e_sum2) < state_threshold) or
                (np.abs(e_sum1) < state_threshold and np.abs(e_sum2) > 1-state_threshold)) and \
                 np.abs(np.abs(e_sum1) - np.abs(h_sum1)) < state_threshold:

                if (e_sum1 - e_sum2) > 0:
                    label = '01'
                else:
                    label = '10'

                state_order.append(label)

            if ((np.abs(np.abs(e_sum1)) > 1-state_threshold and np.abs(h_sum1) < state_threshold) or
                (np.abs(np.abs(e_sum1)) < state_threshold and np.abs(h_sum1) > 1-state_threshold)) and \
                np.abs(np.abs(e_sum1) - np.abs(h_sum1)) > 1-state_threshold:

                if (e_sum1 - e_sum2) > 0:
                    label = 'CA'
                else:
                    label = 'AC'
                state_order.append(label)

        if print_data:
            print ('state: {}  {}'.format(state_number, label))

    # AdiabatH
    if len(state_order) < 4:
        raise Exception('Some states not found')

    if print_data:
        print ('-------------------------------')
    adiabatic_energies = {}
    for m in re.finditer('showmatrix adiabatH', output):
        #print output[m.end():m.end()+800].split()
        data = output[m.end():m.end()+30].split()
        indices = output[m.end()+1:m.end()+4].split(',')
        indices = [int(i) for i in indices]
        # print (indices)
        adiabat_h = float(data[2])
        # print (adiabat_h)
        if indices[0] == indices[1]:
            adiabatic_energies.update({'E_{}'.format(indices[0]+1): adiabat_h * 27.2114})

    if print_data:
        for item in sorted(adiabatic_energies):
            print('{:5} : {:10.5f}'.format(item, adiabatic_energies[item]))

        print('-------------------------------')

    if print_data:
        for item in sorted(adiabatic_energies):
            print('{:5} : {:10.5f}'.format(item, adiabatic_energies[item]))
        print('-------------------------------')

    diabatic_energies = {}
    for m in re.finditer('showmatrix diabatH', output):
        data = output[m.end():m.end()+30].split()
        indices = output[m.end()+1:m.end()+4].split(',')
        indices = [int(i) for i in indices]
        diabat_h = float(data[2])
        if indices[0] == indices[1]:
            # diabatic_energies.append(adiabat_h * 27.2114)
            if state_order[indices[0]] == '10':
                diabatic_energies.update({'E_10': diabat_h * 27.2114})
            if state_order[indices[0]] == '01':
                diabatic_energies.update({'E_01': diabat_h * 27.2114})
            if state_order[indices[0]] == 'AC':
                diabatic_energies.update({'E_AC': diabat_h * 27.2114})
            if state_order[indices[0]] == 'CA':
                diabatic_energies.update({'E_CA': diabat_h * 27.2114})

        diabatic_energies.update({'E_{}_{}'.format(*np.array(indices) + 1): diabat_h * 27.2114})

        # print ('test_indices:', indices, state_order, loc_diab)
        # print ('test_states:', [state_order[indices[0]], state_order[indices[1]]])

        if [state_order[indices[0]], state_order[indices[1]]] == ['10', '01']:
            diabatic_energies.update({'V_DC_1': diabat_h * 27.2114})

        if [state_order[indices[0]], state_order[indices[1]]] == ['01', '10']:
            diabatic_energies.update({'V_DC_2': diabat_h * 27.2114})

        if [state_order[indices[0]], state_order[indices[1]]] == ['AC', 'CA']:
            diabatic_energies.update({'V_CT_1': diabat_h * 27.2114})

        if [state_order[indices[0]], state_order[indices[1]]] == ['CA', 'AC']:
            diabatic_energies.update({'V_CT_2': diabat_h * 27.2114})

        if [state_order[indices[0]], state_order[indices[1]]] == ['10', 'CA']:
            diabatic_energies.update({'V_e_1': diabat_h * 27.2114})

        if [state_order[indices[0]], state_order[indices[1]]] == ['01', 'AC']:
            diabatic_energies.update({'V_e_2': diabat_h * 27.2114})

        if [state_order[indices[0]], state_order[indices[1]]] == ['10', 'AC']:
            diabatic_energies.update({'V_h_1': diabat_h * 27.2114})

        if [state_order[indices[0]], state_order[indices[1]]] == ['01', 'CA']:
            diabatic_energies.update({'V_h_2': diabat_h * 27.2114})

    diabatic_energies.update({'V_DC': [diabatic_energies['V_DC_1'],
                                       diabatic_energies['V_DC_2']]})
    del diabatic_energies['V_DC_1']
    del diabatic_energies['V_DC_2']

    diabatic_energies.update({'V_CT': [diabatic_energies['V_CT_1'],
                                       diabatic_energies['V_CT_2']]})
    del diabatic_energies['V_CT_1']
    del diabatic_energies['V_CT_2']

    # arrange diabatic energies in lists
    diabatic_energies.update({'E_LE': (diabatic_energies['E_10'] +
                                       diabatic_energies['E_01'])/2})

    diabatic_energies.update({'dE_LE': (diabatic_energies['E_10'] -
                                        diabatic_energies['E_01'])/2})

    diabatic_energies.update({'E_CT': np.average((diabatic_energies['E_AC'] +
                                                  diabatic_energies['E_CA'])/2)})

    diabatic_energies.update({'dE_CT': (diabatic_energies['E_AC'] -
                                        diabatic_energies['E_CA'])/2})

    diabatic_energies.update({'V_e': [diabatic_energies['V_e_1'],
                                      diabatic_energies['V_e_2']]})
    del diabatic_energies['V_e_1']
    del diabatic_energies['V_e_2']

    diabatic_energies.update({'V_h': [diabatic_energies['V_h_1'],
                                      diabatic_energies['V_h_2']]})
    del diabatic_energies['V_h_1']
    del diabatic_energies['V_h_2']


    if len(diabatic_energies) < 6:
        raise Exception('diabatic states not resolved!')

    # Showmatrix Decomposed diabatics
    for diabatic_type in ['One', 'K', 'J']:
        for m in re.finditer('showmatrix Decomposed_{}_diabatic'.format(diabatic_type), output):
            data = output[m.end():m.end()+30].split()
            indices = output[m.end()+1:m.end()+4].split(',')
            indices = [int(i) for i in indices]
            element_diabatic = float(data[2])
            if indices[0] == indices[1]:
                if state_order[indices[0]] == '10':
                    diabatic_energies.update({'{}_10'.format(diabatic_type): element_diabatic * 27.2114})
                if state_order[indices[0]] == '01':
                    diabatic_energies.update({'{}_01'.format(diabatic_type): element_diabatic * 27.2114})
                if state_order[indices[0]] == 'AC':
                    diabatic_energies.update({'{}_AC'.format(diabatic_type): element_diabatic * 27.2114})
                if state_order[indices[0]] == 'CA':
                    diabatic_energies.update({'{}_CA'.format(diabatic_type): element_diabatic * 27.2114})

            if [state_order[indices[0]], state_order[indices[1]]] == ['10', '01']:
                diabatic_energies.update({'{}_DC_1'.format(diabatic_type): element_diabatic * 27.2114})

            if [state_order[indices[0]], state_order[indices[1]]] == ['01', '10']:
                diabatic_energies.update({'{}_DC_2'.format(diabatic_type): element_diabatic * 27.2114})

            if [state_order[indices[0]], state_order[indices[1]]] == ['AC', 'CA']:
                diabatic_energies.update({'{}_CT_1'.format(diabatic_type): element_diabatic * 27.2114})

            if [state_order[indices[0]], state_order[indices[1]]] == ['CA', 'AC']:
                diabatic_energies.update({'{}_CT_2'.format(diabatic_type): element_diabatic * 27.2114})

            if [state_order[indices[0]], state_order[indices[1]]] == ['10', 'CA']:
                diabatic_energies.update({'{}_e_1'.format(diabatic_type): element_diabatic * 27.2114})

            if [state_order[indices[0]], state_order[indices[1]]] == ['01', 'AC']:
                diabatic_energies.update({'{}_e_2'.format(diabatic_type): element_diabatic * 27.2114})

            if [state_order[indices[0]], state_order[indices[1]]] == ['10', 'AC']:
                diabatic_energies.update({'{}_h_1'.format(diabatic_type): element_diabatic * 27.2114})

            if [state_order[indices[0]], state_order[indices[1]]] == ['01', 'CA']:
                diabatic_energies.update({'{}_h_2'.format(diabatic_type): element_diabatic * 27.2114})

        diabatic_energies.update({'{}_DC'.format(diabatic_type): [diabatic_energies['{}_DC_1'.format(diabatic_type)],
                                                                  diabatic_energies['{}_DC_2'.format(diabatic_type)]]})
        del diabatic_energies['{}_DC_1'.format(diabatic_type)]
        del diabatic_energies['{}_DC_2'.format(diabatic_type)]

        diabatic_energies.update({'{}_CT'.format(diabatic_type): [diabatic_energies['{}_CT_1'.format(diabatic_type)],
                                                                  diabatic_energies['{}_CT_2'.format(diabatic_type)]]})
        del diabatic_energies['{}_CT_1'.format(diabatic_type)]
        del diabatic_energies['{}_CT_2'.format(diabatic_type)]

        diabatic_energies.update({'{}_e'.format(diabatic_type): [diabatic_energies['{}_e_1'.format(diabatic_type)],
                                                                 diabatic_energies['{}_e_2'.format(diabatic_type)]]})
        del diabatic_energies['{}_e_1'.format(diabatic_type)]
        del diabatic_energies['{}_e_2'.format(diabatic_type)]

        diabatic_energies.update({'{}_h'.format(diabatic_type): [diabatic_energies['{}_h_1'.format(diabatic_type)],
                                                                 diabatic_energies['{}_h_2'.format(diabatic_type)]]})
        del diabatic_energies['{}_h_1'.format(diabatic_type)]
        del diabatic_energies['{}_h_2'.format(diabatic_type)]


    if print_data:
        for item in sorted(diabatic_energies):
            try:
                print('{:5} : {:10.5f}'.format(item, diabatic_energies[item]))
            except:
                print('{:5} : '.format(item) + ' '.join(['{:10.5f}'.format(num) for num in diabatic_energies[item]]))

        print ('-------------------------------')

    coefficients = {}
    for m in re.finditer('showmatrix final adiabatic -> diabatic RotMatrix', output):
        data = output[m.end():m.end()+30].split()
        # print data
        indices = output[m.end()+1:m.end()+4].split(',')
        indices = [int(i) for i in indices]
        adiabat_f = float(data[2])

        for e_index in range(4):
            if [indices[0], state_order[indices[1]]] == [e_index, '10']:
                coefficients.update({'S{}_C10'.format(e_index+1): adiabat_f})
            if [indices[0], state_order[indices[1]]] == [e_index, '01']:
                coefficients.update({'S{}_C01'.format(e_index+1): adiabat_f})
            if [indices[0], state_order[indices[1]]] == [e_index, 'CA']:
                coefficients.update({'S{}_CCA'.format(e_index+1): adiabat_f})
            if [indices[0], state_order[indices[1]]] == [e_index, 'AC']:
                coefficients.update({'S{}_CAC'.format(e_index+1): adiabat_f})

            coefficients.update({'S{}_C{}'.format(*np.array(indices) + 1): adiabat_f})

    coefficients.update({'S_10': [coefficients.pop('S{}_C10'.format(i+1)) for i in range(4)]})
    coefficients.update({'S_01': [coefficients.pop('S{}_C01'.format(i+1)) for i in range(4)]})
    coefficients.update({'S_CA': [coefficients.pop('S{}_CCA'.format(i+1)) for i in range(4)]})
    coefficients.update({'S_AC': [coefficients.pop('S{}_CAC'.format(i+1)) for i in range(4)]})

    if print_data:
        for item in sorted(coefficients):
            try:
                print('{:5} : {:10.5f}'.format(item, coefficients[item]))
            except:
                print('{:5} : '.format(item) + ' '.join(['{:10.5f}'.format(num) for num in coefficients[item]]))

    diabatic_contributions = {}
    w_dc = []
    w_ct = []
    w_e = []
    w_h = []
    for i in range(4):
        factor = [coefficients['S_10'][i] * coefficients['S_01'][i],
                  coefficients['S_01'][i] * coefficients['S_10'][i]]
        w_dc.append(np.average([diabatic_energies['V_DC'][0] * np.sign(factor[0]),
                                diabatic_energies['V_DC'][1] * np.sign(factor[1])]))

        factor = [coefficients['S_AC'][i] * coefficients['S_CA'][i],
                  coefficients['S_CA'][i] * coefficients['S_AC'][i]]
        w_ct.append(np.average([diabatic_energies['V_CT'][0] * np.sign(factor[0]),
                                diabatic_energies['V_CT'][1] * np.sign(factor[1])]))

        factor = [coefficients['S_10'][i] * coefficients['S_AC'][i],
                  coefficients['S_01'][i] * coefficients['S_CA'][i]]
        w_h.append(np.average([diabatic_energies['V_h'][0] * np.sign(factor[0]),
                               diabatic_energies['V_h'][1] * np.sign(factor[1])]))

        factor = [coefficients['S_10'][i] * coefficients['S_CA'][i],
                  coefficients['S_01'][i] * coefficients['S_AC'][i]]
        w_e.append(np.average([diabatic_energies['V_e'][0] * np.sign(factor[0]),
                               diabatic_energies['V_e'][1] * np.sign(factor[1])]))

    diabatic_contributions.update({'W_DC': w_dc})
    diabatic_contributions.update({'W_CT': w_ct})
    diabatic_contributions.update({'W_h': w_h})
    diabatic_contributions.update({'W_e': w_e})

    if print_data:
        print ('-------------------------------')
        for item in ['W_DC', 'W_CT', 'W_h', 'W_e']:
            print('{:5} : {:10.5f} {:10.5f} {:10.5f} {:10.5f}'.format(
                *[item] + diabatic_contributions[item] ))

    for diabatic_type in ['One', 'K', 'J']:
        element_dc = []
        element_ct = []
        element_e = []
        element_h = []
        for i in range(4):
            factor = [coefficients['S_10'][i] * coefficients['S_01'][i],
                      coefficients['S_01'][i] * coefficients['S_10'][i]]
            element_dc.append(np.average([diabatic_energies['{}_DC'.format(diabatic_type)][0] * np.sign(factor[0]),
                                          diabatic_energies['{}_DC'.format(diabatic_type)][1] * np.sign(factor[1])]))

            factor = [coefficients['S_AC'][i] * coefficients['S_CA'][i],
                      coefficients['S_CA'][i] * coefficients['S_AC'][i]]
            element_ct.append(np.average([diabatic_energies['{}_CT'.format(diabatic_type)][0] * np.sign(factor[0]),
                                          diabatic_energies['{}_CT'.format(diabatic_type)][1] * np.sign(factor[1])]))

            factor = [coefficients['S_10'][i] * coefficients['S_AC'][i],
                      coefficients['S_01'][i] * coefficients['S_CA'][i]]
            element_h.append(np.average([diabatic_energies['{}_h'.format(diabatic_type)][0] * np.sign(factor[0]),
                                         diabatic_energies['{}_h'.format(diabatic_type)][1] * np.sign(factor[1])]))

            factor = [coefficients['S_10'][i] * coefficients['S_CA'][i],
                      coefficients['S_01'][i] * coefficients['S_AC'][i]]
            element_e.append(np.average([diabatic_energies['{}_e'.format(diabatic_type)][0] * np.sign(factor[0]),
                                         diabatic_energies['{}_e'.format(diabatic_type)][1] * np.sign(factor[1])]))

        diabatic_contributions.update({'W_{}_DC'.format(diabatic_type): element_dc})
        diabatic_contributions.update({'W_{}_CT'.format(diabatic_type): element_ct})
        diabatic_contributions.update({'W_{}_h'.format(diabatic_type): element_h})
        diabatic_contributions.update({'W_{}_e'.format(diabatic_type): element_e})

        if print_data:
            print('-------------------------------')
            for item in ['W_{}_DC'.format(diabatic_type),
                         'W_{}_CT'.format(diabatic_type),
                         'W_{}_h'.format(diabatic_type),
                         'W_{}_e'.format(diabatic_type)]:
                print('{:5} : {:10.5f} {:10.5f} {:10.5f} {:10.5f}'.format(
                    *[item] + diabatic_contributions[item]))

    l2 = []
    d2le = []
    d2ct = []
    for i in range(4):
        l2.append(coefficients['S_AC'][i]**2 + coefficients['S_CA'][i]**2)
        d2le.append(coefficients['S_10'][i]**2 - coefficients['S_01'][i]**2)
        d2ct.append(coefficients['S_AC'][i]**2 - coefficients['S_CA'][i]**2)

    asymmetries = {}
    asymmetries.update({'lambda2': l2})
    asymmetries.update({'delta2_le': d2le})
    asymmetries.update({'delta2_ct': d2ct})


    if print_data:
        print ('-------------------------------')
        for i, _l in enumerate(l2):
            print ('{}: lambda^2: {:10.5f}  c^2_le: {:10.5f}  c^2_ct: {:10.5f}'.format(i+1, _l, d2le[i], d2ct[i]))

    omega_dc = []
    for i in range(4):
        omega_dc.append(np.sqrt((1 - l2[i]) ** 2 - d2le[i] ** 2) * diabatic_contributions['W_DC'][i])

    omega_ct = []
    for i in range(4):
        omega_ct.append(np.sqrt(l2[i] ** 2 - d2ct[i] ** 2) * diabatic_contributions['W_CT'][i])

    omega_e = []
    for i in range(4):
        omega_e.append(2 * (coefficients['S_10'][i] * coefficients['S_CA'][i] * diabatic_energies['V_e'][0] +
                            coefficients['S_01'][i] * coefficients['S_AC'][i] * diabatic_energies['V_e'][1]))

    omega_h = []
    for i in range(4):
        omega_h.append(2 * (coefficients['S_10'][i] * coefficients['S_AC'][i] * diabatic_energies['V_h'][0] +
                            coefficients['S_01'][i] * coefficients['S_CA'][i] * diabatic_energies['V_h'][1]))

    diabatic_contributions.update({'Omega_dc': omega_dc})
    diabatic_contributions.update({'Omega_ct': omega_ct})

    diabatic_contributions.update({'Omega_e': omega_e})
    diabatic_contributions.update({'Omega_h': omega_h})

    for diabatic_type in ['One', 'K', 'J']:

        omega_e_dc = []
        for i in range(4):
            omega_e_dc.append(np.sqrt((1 - l2[i]) ** 2 - d2le[i] ** 2) *
                            diabatic_contributions['W_{}_DC'.format(diabatic_type)][i])

        omega_e_ct = []
        for i in range(4):
            omega_e_ct.append(np.sqrt(l2[i] ** 2 - d2ct[i] ** 2) *
                            diabatic_contributions['W_{}_CT'.format(diabatic_type)][i])

        omega_e_e = []
        for i in range(4):
            omega_e_e.append(2 * (coefficients['S_10'][i] * coefficients['S_CA'][i] *
                                  diabatic_energies['{}_e'.format(diabatic_type)][0] +
                                  coefficients['S_01'][i] * coefficients['S_AC'][i] *
                                  diabatic_energies['{}_e'.format(diabatic_type)][1]))
        omega_e_h = []
        for i in range(4):
            omega_e_h.append(2 * (coefficients['S_10'][i] * coefficients['S_AC'][i] *
                                  diabatic_energies['{}_h'.format(diabatic_type)][0] +
                                  coefficients['S_01'][i] * coefficients['S_CA'][i] *
                                  diabatic_energies['{}_h'.format(diabatic_type)][1]))

        diabatic_contributions.update({'Omega_{}_dc'.format(diabatic_type): omega_e_dc})
        diabatic_contributions.update({'Omega_{}_ct'.format(diabatic_type): omega_e_ct})

        diabatic_contributions.update({'Omega_{}_e'.format(diabatic_type): omega_e_e})
        diabatic_contributions.update({'Omega_{}_h'.format(diabatic_type): omega_e_h})

    reordering = []
    for m in re.finditer('Reordering necessary!', output):
        reordering.append([int(num) for num in output[m.end():m.end()+7].split('->')])

    return {'loc_diab': loc_diab,
            'adiabatic_energies': adiabatic_energies,
            'diabatic_energies': diabatic_energies,
            'diabatic_contributions': diabatic_contributions,
            'coefficients': coefficients,
            'asymmetries': asymmetries,
            'reordering': reordering}

import re
import numpy as np
import os


def analyze_diabatic(output, print_data=False, state_threshold=0.2, n_mon=6):

    # Coordinates
    n = output.find('$molecule')
    n2 = output[n:].find('$end')
    molecule_region = output[n:n+n2-1].replace('\t', ' ').split('\n')[1:]
    coordinates = np.array([ np.array(line.split()[1:4], dtype=float) for line in molecule_region[1:]])
    symbols = [line.split()[0].capitalize() for line in molecule_region[1:]]

    n_atoms = len(coordinates)


    # Diabatic

    n = output.find('On the next line, list which excited')
    loc_diab = output[n:n+100].split()[13:17]

    if print_data:
        # Mulliken
        print ('-------------------------------')

        print ('Mulliken')
    state_order = []
    for m in re.finditer('Mulliken analysis of TDA State', output):
        #print output[m.end():m.end()+800].split()
        data = output[m.end():m.end()+800].split()
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
    if len(state_order) < len(loc_diab) or len(state_order) < 4:
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

        print ('-------------------------------')

    if print_data:
        for item in sorted(adiabatic_energies):
            print('{:5} : {:10.5f}'.format(item, adiabatic_energies[item]))
        print ('-------------------------------')

    diabatic_energies = {}
    for m in re.finditer('showmatrix diabatH', output):
        data = output[m.end():m.end()+30].split()
        indices = output[m.end()+1:m.end()+4].split(',')
        indices = [int(i) for i in indices]
        adiabat_h = float(data[2])
        if indices[0] == indices[1]:
            # diabatic_energies.append(adiabat_h * 27.2114)
            if indices[0] == 0:
                diabatic_energies.update({'E_LE_1': adiabat_h * 27.2114})
            if indices[0] == 1:
                diabatic_energies.update({'E_LE_2': adiabat_h * 27.2114})
            if indices[0] == 2:
                diabatic_energies.update({'E_CT_1': adiabat_h * 27.2114})
            if indices[0] == 3:
                diabatic_energies.update({'E_CT_2': adiabat_h * 27.2114})

        diabatic_energies.update({'E_{}_{}'.format(*np.array(indices) + 1): adiabat_h * 27.2114})

        # print ('test_indices:', indices, state_order, loc_diab)
        # print ('test_states:', [state_order[indices[0]], state_order[indices[1]]])

        if [state_order[indices[0]], state_order[indices[1]]] == ['01', '10']:
            diabatic_energies.update({'V_DC': adiabat_h * 27.2114})

        if [state_order[indices[0]], state_order[indices[1]]] == ['AC', 'CA']:
            diabatic_energies.update({'V_CT': adiabat_h * 27.2114})

        if [state_order[indices[0]], state_order[indices[1]]] == ['10', 'CA']:
            diabatic_energies.update({'V_e_1': adiabat_h * 27.2114})

        if [state_order[indices[0]], state_order[indices[1]]] == ['01', 'AC']:
            diabatic_energies.update({'V_e_2': adiabat_h * 27.2114})

        if [state_order[indices[0]], state_order[indices[1]]] == ['10', 'AC']:
            diabatic_energies.update({'V_h_1': adiabat_h * 27.2114})

        if [state_order[indices[0]], state_order[indices[1]]] == ['01', 'CA']:
            diabatic_energies.update({'V_h_2': adiabat_h * 27.2114})

    # arrange diabatic energies in lists
    diabatic_energies.update({'E_LE': [diabatic_energies['E_LE_1'],
                                       diabatic_energies['E_LE_2']]})
    del diabatic_energies['E_LE_1']
    del diabatic_energies['E_LE_2']

    diabatic_energies.update({'E_CT': [diabatic_energies['E_CT_1'],
                                       diabatic_energies['E_CT_2']]})
    del diabatic_energies['E_CT_1']
    del diabatic_energies['E_CT_2']

    diabatic_energies.update({'V_e': [diabatic_energies['V_e_1'],
                                      diabatic_energies['V_e_2']]})
    del diabatic_energies['V_e_1']
    del diabatic_energies['V_e_2']

    diabatic_energies.update({'V_h': [diabatic_energies['V_h_1'],
                                      diabatic_energies['V_h_2']]})
    del diabatic_energies['V_h_1']
    del diabatic_energies['V_h_2']

    if print_data:
        for item in sorted(diabatic_energies):
            try:
                print('{:5} : {:10.5f}'.format(item, diabatic_energies[item]))
            except:
                print('{:5} : '.format(item) + ' '.join(['{:10.5f}'.format(num) for num in diabatic_energies[item]]))

        print ('-------------------------------')

    if len(diabatic_energies) < 6:
        raise Exception('diabatic states not resolved!')

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
            print('{:5} : {:10.5f}'.format(item, coefficients[item]))

    diabatic_contributions = {}
    w_dc = []
    w_ct = []
    w_e = []
    w_h = []

    for i in range(4):
        factor = coefficients['S_10'][i] * coefficients['S_01'][i]
        w_dc.append(diabatic_energies['V_DC'] * np.sign(factor))

        factor = coefficients['S_CA'][i] * coefficients['S_AC'][i]
        w_ct.append(diabatic_energies['V_CT'] * np.sign(factor))

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
                *[item] + [diabatic_contributions[item+'_{}'.format(i+1)] for i in range(4)] ))

    l = []
    for i in range(4):
        l.append(np.sqrt(2) * np.abs(coefficients['S_AC'][i]))

    #l = np.sqrt(2) * coefficients['S1_C3']
    #l_prima = np.sqrt(2) * coefficients['S2_C3']

    if print_data:
        print ('-------------------------------')
        for i, _l in enumerate(l):
            print ('lambda {}: {:10.5f}'.format(i+1, _l))

    reordering = []
    for m in re.finditer('Reordering necessary!', output):
        reordering.append([int(num) for num in output[m.end():m.end()+7].split('->')])

    return {'loc_diab': loc_diab,
            'adiabatic_energies': adiabatic_energies,
            'diabatic_energies': diabatic_energies,
            'diabatic_contributions': diabatic_contributions,
            'coefficients': coefficients,
            'lambda': l,
            'reordering' : reordering}


# Start script -------------------------------------------------------------

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    dir_list = os.listdir('data/eclipsed')

    total_data = []
    for dir in dir_list:
        try:
            data = analyze_diabatic('data/eclipsed/{}/diabatic.out'.format(dir), print_data=True)
            print (data)
            if data is not None:
                total_data.append(data)
        except IOError:
            pass

    x = [data['c_coordinates'][0] for data in total_data]
    indices =  np.argsort(x)
    x = np.array(x)[indices]

    l = np.array([data['lambda'] for data in total_data])[indices].T


    plt.figure(1)
    for i, _l in enumerate(l):
        plt.plot(x, _l, label='lambda {}'.format(i+1))

    plt.xlabel('distance (A)')
    plt.ylabel('unidades')
    plt.legend()
    #plt.show()



    plt.figure(2)
    ll = []
    for i, lam in enumerate(data['lambda']):
        ll.append(np.array([2*lam * np.sqrt(1 - lam**2) *
                      (data['diabatic_contributions']['W_e_{}'.format(i+1)] + data['diabatic_contributions']['W_h_{}'.format(i+1)])
                      for data in total_data])[indices])

    ll = np.array(ll)

    for i, lam in enumerate(ll):
        plt.plot(x, lam, label='ll {}'.format(i+1))
    # plt.plot(x, ll2, label='ll2')
    plt.xlabel('distance (A)')
    plt.ylabel('unidades')
    plt.legend()
    #plt.show()


    plt.figure(3)
    for i in range(4):
        label = 'E_{}'.format(i+1)
        e = np.array([data['adiabatic_energies'][label] for data in total_data])[indices]
        plt.plot(x, e, label=label)

    plt.xlabel('distance (A)')
    plt.ylabel('unidades')
    plt.legend()
    #plt.show()

    plt.figure(4)

    d = {}
    for e in np.array([data['diabatic_energies'] for data in total_data])[indices]:
        for label in e.keys():
            try:
                d[label].append(e[label])
            except:
                d[label] = [e[label]]

    print(d)

    for label in ['E_LE', 'E_CT']:
        plt.plot(x, d[label] , label=label)
    plt.xlabel('distance (A)')
    plt.ylabel('unidades')
    plt.legend()
    #plt.show()

    plt.figure(5)
    d2 = {}
    for e in np.array([data['diabatic_contributions'] for data in total_data])[indices]:
        for label in e.keys():
            try:
                d2[label].append(e[label])
            except:
                d2[label] = [e[label]]

    print(2)

    for label in d2.keys():
        if label.endswith('_1'):
            plt.plot(x, d2[label] , label=label)
    plt.xlabel('distance (A)')
    plt.ylabel('unidades')
    plt.legend()


    plt.figure(6)
    for label in d2.keys():
        if label.endswith('_2'):
            plt.plot(x, d2[label] , label=label)
    plt.xlabel('distance (A)')
    plt.ylabel('unidades')
    plt.legend()


    plt.figure(7)

    for i, _l in enumerate(l):
        e2 = np.array(_l)**2 * (np.array(d['E_CT']) - np.array(d['E_LE']) +
                                 np.array(d2['W_CT_{}'.format(i+1)]) - np.array(d2['W_DC_{}'.format(i+1)]))
    #       plt.plot(x, e1_, label='E1-')
        plt.plot(x, np.array(e2), label='E2-{}'.format(i+1))

    plt.xlabel('distance (A)')
    plt.ylabel('unidades')
    plt.legend()

    plt.show()



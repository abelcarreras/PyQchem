from pyqchem.qchem_core import parse_output
import numpy as np
import matplotlib.pyplot as plt

from pyqchem.parsers.parser_diabatic import analyze_diabatic
from pyqchem.parsers.basic import basic_parser_qchem
# common qchem input parameters


def get_order_states(states, epsilon=1):

    print([state['total energy'] for state in states])
    print([state['transition moment'] for state in states])
    print([state['transition moment'][1] for state in states])

    order = np.argsort([state['total energy'] for state in states]).tolist()

    for i, tm in enumerate([state['transition moment'] for state in states]):
        if np.linalg.norm(tm) > epsilon:
            index = order.pop(i)
            order = [index] + order

    order = np.argsort([np.abs(state['transition moment'][1]) for state in states]).tolist()[::-1]
    order[1], order[3] = order[3], order[1]
#    order[1], order[0] = order[0], order[1]

    print('--------------order------------\n', order)
    return order


def correct_order(list, order):

    alist = np.array(list)
    orddered_list = alist[order]

    return orddered_list


@parse_output
def get_output(filename):
    with open(filename, 'r') as f:
        output = f.read()
    return output, ''

dif1 = []
dif2 = []
dif3 = []

ener = []
ener1 = []
ener2 = []
ener3 = []

kill = []
kill2 = []
kill3 = []
kill4 = []
kill5 = []

points = []
points_total = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4]
for num in points_total:
    file = 'test_outputs_3/qchemout_0_{}.out'.format(num)

    print('file: {}'.format(file))
    data1 = get_output(file, parser=basic_parser_qchem, force_recalculation=True)

    try:
        states = np.array(get_output(file).split('\n')[180].split(), dtype=int)-1
    except ValueError:
        states = np.array(get_output(file).split('\n')[171].split(), dtype=int)-1

    # get list of interesting states
    states_info = []
    for state in states:
        states_info.append(data1['excited states cis'][state])
    print(states_info)

    order = get_order_states(states_info)
    print(states)

    e = np.array([data1['excited states cis'][state]['excitation energy'] for state in states])[order]
    kill3.append(e)

    t = np.array([np.linalg.norm(data1['excited states cis'][state]['transition moment']) for state in states])[order]
    kill4.append(t)


    # get energy contributions
    try:
        data = get_output(file, parser=analyze_diabatic, force_recalculation=True)
    except Exception:
        print('Failed!')
        continue

    coefficients = np.zeros([4, 4])
    adiabatic_energies = np.zeros(4)
    diabatic_energies = np.zeros([4, 4])
    for i in range(4):
        for j in range(4):
            coefficients[i, j] = data['coefficients']['S{}_C{}'.format(i + 1, j + 1)]
            diabatic_energies[i, j ] = data['diabatic_energies']['E_{}_{}'.format(i+1, j+1)]

        adiabatic_energies[i] = data['adiabatic_energies']['E_{}'.format(i + 1)]

    print('coefficients')
    print(coefficients)
    print('adiabatic energies')
    print(adiabatic_energies)
    print('diabatic')
    print(diabatic_energies)

    print('recalculated adiabatic energies')
    calculated_adiabatic = np.diag(np.dot(coefficients, np.dot(diabatic_energies, coefficients.T)))
    print(calculated_adiabatic)

    calculated_adiabatic = correct_order(calculated_adiabatic, order)
    adiabatic_energies = correct_order(adiabatic_energies, order)

    print('diference adiabatic (original - calculated)')
    print(calculated_adiabatic - adiabatic_energies)

    dif1.append(calculated_adiabatic - adiabatic_energies)
    ener1.append(calculated_adiabatic)
    ener.append(adiabatic_energies)

    print('-------------------------------------')

    diabatic_energies_2 = np.zeros([4, 4])

    diabatic_energies_2[0, 0] = np.average(data['diabatic_energies']['E_LE'])
    diabatic_energies_2[2, 2] = np.average(data['diabatic_energies']['E_CT'])

    diabatic_energies_2[0, 1] = data['diabatic_energies']['V_DC']
    diabatic_energies_2[2, 3] = data['diabatic_energies']['V_CT']

    diabatic_energies_2[1, 2] = data['diabatic_energies']['V_e'][0]
    diabatic_energies_2[1, 3] = data['diabatic_energies']['V_h'][0]

    # prima
    diabatic_energies_2[1, 1] = np.average(data['diabatic_energies']['E_LE'])
    diabatic_energies_2[3, 3] = np.average(data['diabatic_energies']['E_CT'])

    diabatic_energies_2[0, 3] = data['diabatic_energies']['V_e'][1]
    diabatic_energies_2[0, 2] = data['diabatic_energies']['V_h'][1]

    # fill lower diagonal matrix
    for i in range(4):
        for j in range(i):
            diabatic_energies_2[i, j] = diabatic_energies_2[j, i]

    # Diabatic matrix
    #   E_LE   V_DC   V_h'   V_e'
    #   V_DC   E_LE'  V_e    V_h
    #   V_h'   V_e    E_CT   V_e
    #   V_e'   V_h    V_e    E_CT'

    print('diabatic_energies_2')
    print(diabatic_energies_2)

    coefficients_2 = [data['coefficients']['S_01'],
                      data['coefficients']['S_10'],
                      data['coefficients']['S_CA'],
                      data['coefficients']['S_AC']]

    coefficients_2 = np.array(coefficients_2).T

    print('coefficients_2')
    print(coefficients_2)

    calculated_adiabatic_2 = np.diag(np.dot(coefficients_2, np.dot(diabatic_energies_2, coefficients_2.T)))
    print(calculated_adiabatic_2)

    calculated_adiabatic_2 = correct_order(calculated_adiabatic_2, order)

    print('diference adiabatic_2 (original - calculated)')
    print(calculated_adiabatic_2 - adiabatic_energies)

    dif2.append(calculated_adiabatic_2 - adiabatic_energies)
    ener2.append(calculated_adiabatic_2)

    print('-------------------------------------')

    calculated_adiabatic_3 = []
    w_e_list = []
    w_h_list = []
    lambda_list = []

    e_le = np.average(data['diabatic_energies']['E_LE'])
    e_ct = np.average(data['diabatic_energies']['E_CT'])

    for i in range(4):

        lmb = data['lambda'][i]

        w_dc = data['diabatic_energies']['V_DC'] * np.sign(data['coefficients']['S_10'][i] *
                                                           data['coefficients']['S_01'][i])

        w_ct = data['diabatic_energies']['V_CT'] * np.sign(data['coefficients']['S_CA'][i] *
                                                           data['coefficients']['S_AC'][i])

        w_h1 = data['diabatic_energies']['V_h'][0] * np.sign(data['coefficients']['S_10'][i] *
                                                             data['coefficients']['S_AC'][i])

        w_h2 = data['diabatic_energies']['V_h'][1] * np.sign(data['coefficients']['S_01'][i] *
                                                             data['coefficients']['S_CA'][i])

        w_e1 = data['diabatic_energies']['V_e'][0] * np.sign(data['coefficients']['S_10'][i] *
                                                             data['coefficients']['S_CA'][i])

        w_e2 = data['diabatic_energies']['V_e'][1] * np.sign(data['coefficients']['S_01'][i] *
                                                             data['coefficients']['S_AC'][i])

        w_e = np.average([w_e1, w_e2])
        w_h = np.average([w_h1, w_h2])

        o_e = 2 * lmb * np.sqrt(1 - lmb ** 2) * w_e
        o_h = 2 * lmb * np.sqrt(1 - lmb ** 2) * w_h

        #w_dc = data['diabatic_contributions']['W_DC'][i]
        #w_ct = data['diabatic_contributions']['W_CT'][i]

        #w_e = data['diabatic_contributions']['W_e'][i]
        #w_h = data['diabatic_contributions']['W_h'][i]

        #o_e = data['diabatic_contributions']['Omega_e'][i]
        #o_h = data['diabatic_contributions']['Omega_h'][i]

        ###### manual fix ########
        manual_fix = False
        if manual_fix:
            if num < 0.9:
                kk2 = -1
            else:
                kk2 = 1

            if num > 2.2:
                kk3 = -1
            else:
                kk3 = 1

            if i == order[0]:
                lmb = lmb * kk2
                w_h = w_h * kk2
                w_e = w_e * kk2

            #if j == order[1]:
            #    w_h = w_h * kk2
            #    w_e = w_e * kk2

            #if j == order[2]:
            #    w_h = w_h * kk3
            #    w_e = w_e * kk3

            if i == order[1]:
                lmb = lmb * kk3
                w_h = w_h * kk3
                w_e = w_e * kk3
        ###########

        w_e_list.append(w_e)
        w_h_list.append(w_h)
        lambda_list.append(lmb)

        calculated_adiabatic_3.append(e_le + w_dc + 2 * lmb * np.sqrt(1 - lmb ** 2) * (w_e + w_h) + lmb ** 2 * (e_ct - e_le + w_ct - w_dc))

    calculated_adiabatic_3 = correct_order(calculated_adiabatic_3, order)
    w_e_list = correct_order(w_e_list, order)
    w_h_list = correct_order(w_h_list, order)
    lambda_list = correct_order(lambda_list, order)

    kill.append(w_e_list)
    kill2.append(w_h_list)
    kill5.append(lambda_list)

    ener3.append(calculated_adiabatic_3)
    dif3.append(calculated_adiabatic_3 - adiabatic_energies)

    points.append(num)

plt.figure()
plt.title('excitation energy')
for i, k in enumerate(np.array(kill3).T[:]):
    plt.plot(points_total, k, '-o', label=i)
plt.legend()

plt.figure()
plt.title('transition moment module')
for i, k in enumerate(np.array(kill4).T[:]):
    plt.plot(points_total, k, '-o', label=i)
plt.legend()
plt.show()


#exit()
plt.figure()
plt.title('lambda')

for i, k in enumerate(np.array(kill5).T[2:]):
    plt.plot(points, k, '-o', label=i)
plt.legend()

for i, k in enumerate(np.array(kill5).T[:2]):
    plt.plot(points, k, '--', label=i+2)
plt.legend()

plt.figure()
plt.title('W_e')

for i, k in enumerate(np.array(kill).T[:2]):
    plt.plot(points, k, '--', label=i)
plt.legend()

for i, k in enumerate(np.array(kill).T[2:]):
    plt.plot(points, k, '-o', label=i+2)
plt.legend()


plt.figure()
plt.title('W_h')

for i, k in enumerate(np.array(kill2).T[:2]):
    plt.plot(points, k, '--', label=i)
plt.legend()

for i, k in enumerate(np.array(kill2).T[2:]):
    plt.plot(points, k, '-o', label=i+2)
plt.legend()


plt.show()

#exit()

dif1 = np.array(dif1).T
dif2 = np.array(dif2).T
dif3 = np.array(dif3).T

ener1 = np.array(ener1).T
ener2 = np.array(ener2).T
ener3 = np.array(ener3).T

ener = np.array(ener).T


plt.title('diference adiabatic')
for i, p in enumerate(dif1):
    plt.plot(points, p, '-o',label='{}: {}'.format('Direct', i+1))

for i, p in enumerate(dif2):
    plt.plot(points, p, '-+', label='{}: {}'.format('Separated', i+1))

for i, p in enumerate(dif3):
    plt.plot(points, p, '-x', label='{}: {}'.format('Parts', i+1))

plt.legend()
plt.show()

#exit()

plt.title('Adiabatic energies')
for i, p in enumerate(ener1):
    plt.plot(points, p, '-o', label='{}: {}'.format('Direct', i+1))

for i, p in enumerate(ener2):
    plt.plot(points, p, '-+', label='{}: {}'.format('Separated', i+1))

for i, p in enumerate(ener3):
    plt.plot(points, p, '-x', label='{}: {}'.format('Parts', i+1))

for i, p in enumerate(ener):
    plt.plot(points, p, '--', label='{}: {}'.format('Reference', i+1))

plt.legend()
plt.show()



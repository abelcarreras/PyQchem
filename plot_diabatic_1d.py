import pickle
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl


def quadriplot(data1, data2, data3, data4, label1, label2, label3, label4,
               x_range, title=None, factor=False, range_y=(-1, 1), ylabel='Energy [eV]'):

    if factor:
        factor_h_to_ev = 27.2116

        data1 = factor_h_to_ev * np.array(data1)
        data2 = factor_h_to_ev * np.array(data2)
        data3 = factor_h_to_ev * np.array(data3)
        data4 = factor_h_to_ev * np.array(data4)

    f = plt.figure(1)

    plt.title(title)
    plt.plot(x_range, data1, label=label1)
    plt.plot(x_range, data2, label=label2)
    plt.plot(x_range, data3, label=label3)
    plt.plot(x_range, data4, label=label4)

    plt.xlim(3.5, 6.0)
    plt.ylim(*range_y)
    plt.legend()
    plt.xlabel('Distance X [Å]')
    plt.ylabel(ylabel)

    plt.show()
    return f


def multiplot(data , labels, x_range, title=None, factor=False, range_y=(-1, 1), ylabel='Energy [eV]'):

    if factor:
        factor_h_to_ev = 27.2116
        for i, dat in enumerate(data):
            data[i] = factor_h_to_ev * np.array(dat)

    f = plt.figure(1)

    plt.title(title)
    for dat, label in zip(data, labels):
        plt.plot(x_range, dat, label=label)

    plt.xlim(3.5, 6.0)
    plt.ylim(*range_y)
    plt.legend()
    plt.xlabel('Distance X [Å]')
    plt.ylabel(ylabel)

    plt.show()
    return f


def biplot(data1, data2, label1, label2, x_range, factor=False, range_y=(-1, 1), ylabel='Energy [eV]'):

    if factor:
        factor_h_to_ev = 27.2116

        data1 = factor_h_to_ev * np.array(data1)
        data2 = factor_h_to_ev * np.array(data2)

    f = plt.figure(1)

    plt.title(label1[:-2])
    plt.plot(x_range, data1, label=label1)
    plt.plot(x_range, data2, label=label2)
    plt.xlim(3.5, 6.0)
    if range_y is not None:
        plt.ylim(*range_y)
    plt.legend()
    plt.xlabel('Distance X [Å]')
    plt.ylabel(ylabel)

    plt.show()
    return f


folder = '1d_plot/'

with open('my_data_distance_2.pkl', 'rb') as input:
    calculation_data = pickle.load(input)
    print('Loaded data from calculation_data.pkl')
for slide_x in calculation_data['range_x']:
        if '{}'.format(slide_x) in calculation_data:
            print(calculation_data['{}'.format(slide_x)])
            data_i = calculation_data['{}'.format(slide_x)]


x_range = calculation_data['range_x']


total_data = []
full_range = []
i = 0
x_coordinate = []
for slide_x in x_range:
        print('{}'.format(slide_x))
        if '{}'.format(slide_x) in calculation_data:
            data = calculation_data['{}'.format(slide_x)]
            total_data.append(data)
            full_range.append(i)
            x_coordinate.append(slide_x)
            i += 1
# print(total_data[0])
x_range = np.array(x_coordinate)

###################
# Z, Y = np.meshgrid(z_range, y_range)

data_1 = []
for i, e in enumerate([data['diabatic_contributions']['W_DC_1'] for data in total_data]):
    data_1.append(e)

data_2 = []
for i, e in enumerate([data['diabatic_contributions']['W_DC_2'] for data in total_data]):
    data_2.append(e)


f = biplot(data_1, data_2, 'W_DC_1', 'W_DC_2', x_range)
f.savefig(folder + "W_DC.pdf", bbox_inches='tight')

wdc1 = data_1
wdc2 = data_2


#exit()

####################
data_1 = []
for i, e in enumerate([data['diabatic_contributions']['W_CT_1'] for data in total_data]):
    data_1.append(e)

data_2 = []
for i, e in enumerate([data['diabatic_contributions']['W_CT_2'] for data in total_data]):
    data_2.append(e)

f = biplot(data_1, data_2, 'W_CT_1', 'W_CT_2', x_range)
f.savefig(folder + "W_CT.pdf", bbox_inches='tight')
wct1 = data_1
wct2 = data_2


#####################

data_1 = []
for i, e in enumerate([data['diabatic_contributions']['W_e_1'] for data in total_data]):
    data_1.append(e)

data_2 = []
for i, e in enumerate([data['diabatic_contributions']['W_e_2'] for data in total_data]):
    data_2.append(e)

# correct sign:
for data in [data_1, data_2]:
    for i in range(len(data)):
        s = np.sign(np.average(data))
        if s * np.sign(data[i]) < 0:
            data[i] *= -1


f = biplot(data_1, data_2, 'W_e_1', 'W_e_2', x_range)
f.savefig(folder + "W_e.pdf", bbox_inches='tight')

we1 = data_1
we2 = data_2

########################

data_1 = []
for i, e in enumerate([data['diabatic_contributions']['W_h_1'] for data in total_data]):
    data_1.append(e)

data_2 = []
for i, e in enumerate([data['diabatic_contributions']['W_h_2'] for data in total_data]):
    data_2.append(e)

# correct order:
for i, lc in enumerate(zip(data_1, data_2)):
    l1, l2 = lc
    if l1 > 0:
        data_1[i] *= -1
    if l2 < 0:
        data_2[i] *= -1

f = biplot(data_1, data_2, 'W_h_1', 'W_h_2', x_range)
f.savefig(folder + "W_h.pdf", bbox_inches='tight')

wh1 = data_1
wh2 = data_2

######################
# W per state

f = quadriplot(wh1, we1, wct1, wdc1, 'W_h', 'W_e', 'W_CT', 'W_DC', x_range, title='State 1')
f.savefig(folder + "W_1.pdf", bbox_inches='tight')

f = quadriplot(wh2, we2, wct2, wdc2, 'W_h', 'W_e', 'W_CT', 'W_DC', x_range, title='State 2')
f.savefig(folder + "W_2.pdf", bbox_inches='tight')

#######################


data_1 = []
for i, e in enumerate(np.array([data['lambda'] for data in total_data]).T[0]):
    data_1.append(e)

data_2 = []
for i, e in enumerate(np.array([data['lambda'] for data in total_data]).T[1]):
    data_2.append(e)

f = biplot(data_1, data_2, 'lambda 1', 'lambda 2', x_range, factor=False, range_y=[0, 0.6], ylabel='')
f.savefig(folder + "lambda.pdf", bbox_inches='tight')

l1 = np.array(data_1)
l2 = np.array(data_2)

##########################

l = np.array([data['lambda'] for data in total_data]).T

d_energies = {}
for e in np.array([data['diabatic_energies'] for data in total_data]):
    for label in e.keys():
        try:
            d_energies[label].append(e[label])
        except:
            d_energies[label] = [e[label]]

#f = biplot(d_energies['E_CT'], data_2, 'E1-1', 'E1-2', x_range)




d_contributions = {}
for c in np.array([data['diabatic_contributions'] for data in total_data]):
    for label in c.keys():
        try:
            d_contributions[label].append(c[label])
        except:
            d_contributions[label] = [c[label]]

# data_1 = 2 * np.array(l[0]) * np.square(1 - l[0]**2) * (np.array(d_contributions['W_e_1']) - np.array(d_contributions['W_h_1']))
# data_2 = 2 * np.array(l[1]) * np.square(1 - l[1]**2) * (np.array(d_contributions['W_e_2']) - np.array(d_contributions['W_h_2']))

data_1 = 2 * l1 * np.sqrt(1 - l1**2) * (np.array(we1) + np.array(wh1))
data_2 = 2 * l2 * np.sqrt(1 - l2**2) * (np.array(we2) + np.array(wh2))
print(x_range[5])
print(we1[5])
print(wh1[5])
print(l1[5])
print(data_1[5])


f = biplot(data_1, data_2, 'E1-1', 'E1-2', x_range)
f.savefig(folder + "E1.pdf", bbox_inches='tight')
e1_1 = data_1
e1_2 = data_2

####
data_1 = l1**2 * (np.array(d_energies['E_CT']) - np.array(d_energies['E_LE']) +
                  np.array(d_contributions['W_CT_1']) - np.array(d_contributions['W_DC_1']))

data_2 = l2**2 * (np.array(d_energies['E_CT']) - np.array(d_energies['E_LE']) +
                  np.array(d_contributions['W_CT_2']) - np.array(d_contributions['W_DC_2']))


f = biplot(data_1, data_2, 'E2-1', 'E2-2', x_range)
f.savefig(folder + "E2.pdf", bbox_inches='tight')
e2_1 = data_1
e2_2 = data_2



###
data_1 = []
for i, e in enumerate([data['diabatic_energies']['V_DC'] for data in total_data]):
    data_1.append(e)

data_2 = []
for i, e in enumerate([data['diabatic_energies']['V_CT'] for data in total_data]):
    data_2.append(e)

data_3 = []
for i, e in enumerate([data['diabatic_energies']['V_e'] for data in total_data]):
    data_3.append(e)

data_4 = []
for i, e in enumerate([data['diabatic_energies']['V_h'] for data in total_data]):
    data_4.append(e)


# correct sign:
for data in [data_1, data_2, data_3, data_4]:
    for i in range(len(data)):
        s = np.sign(np.average(data))
        if s * np.sign(data[i]) < 0:
            data[i] *= -1

f = quadriplot(data_1, data_2, data_3, data_4, 'V_DC', 'V_CT', 'V_e', 'V_h', x_range, title='Diabatic contributions')

f.savefig(folder + "diabatic.pdf", bbox_inches='tight')


for state in ['10', '01', 'AC', 'CA']:

    data_1 = []
    for i, e in enumerate([data['coefficients']['S1_C{}'.format(state)] for data in total_data]):
        data_1.append(e)

    data_2 = []
    for i, e in enumerate([data['coefficients']['S2_C{}'.format(state)] for data in total_data]):
        data_2.append(e)

    data_3 = []
    for i, e in enumerate([data['coefficients']['S3_C{}'.format(state)] for data in total_data]):
        data_3.append(e)

    data_4 = []
    for i, e in enumerate([data['coefficients']['S4_C{}'.format(state)] for data in total_data]):
        data_4.append(e)


    # correct sign:
    for data in [data_1, data_2, data_3, data_4]:
        for i in range(len(data)):
            s = np.sign(data[1])
            if s * np.sign(data[i]) < 0:
                data[i] *= -1

    f = quadriplot(data_1, data_2, data_3, data_4, 'S1_{}'.format(state),
                   'S2_{}'.format(state), 'S2_{}'.format(state), 'S3_{}'.format(state),
                   x_range, title='diabatic RotMatrix {}'.format(state), ylabel='')

    f.savefig(folder + "adiabatic_S{}.pdf".format(state), bbox_inches='tight')

########

# gràfica amb W_DC(2), E1_2 i E2_2 en eV.
f = multiplot([wdc1, e1_1, e2_1], ['W_DC', 'E1', 'E2'], x_range, title='State 1')
f.savefig(folder + "contributions_1.pdf", bbox_inches='tight')

f = multiplot([wdc2, e1_2, e2_2], ['W_DC', 'E1', 'E2'], x_range, title='State 2')
f.savefig(folder + "contributions_2.pdf", bbox_inches='tight')



data_1 = []
for i, e in enumerate([data['adiabatic_energies']['E_1'] for data in total_data]):
    data_1.append(e)
data_2 = []
for i, e in enumerate([data['adiabatic_energies']['E_2'] for data in total_data]):
    data_2.append(e)
data_3 = []
for i, e in enumerate([data['adiabatic_energies']['E_3'] for data in total_data]):
    data_3.append(e)
data_4 = []
for i, e in enumerate([data['adiabatic_energies']['E_4'] for data in total_data]):
    data_4.append(e)


f = quadriplot(data_1, data_2, data_3, data_4, 'state_1', 'state_2', 'state_3', 'state_4', x_range,
               title='Adiabatic energies', range_y=[6, 13])
f.savefig(folder + "adiabatic_energies.pdf", bbox_inches='tight')


data_1 = []
for i, e in enumerate([data['diabatic_energies']['E_LE'] for data in total_data]):
    data_1.append(e)
data_2 = []
for i, e in enumerate([data['diabatic_energies']['E_CT'] for data in total_data]):
    data_2.append(e)


f = multiplot([data_1, data_2], ['E_LE', 'E_CT'], x_range, title='Diabatic energies', range_y=[6, 13])
f.savefig(folder + "diabatic_energies.pdf", bbox_inches='tight')

exit()


##############################################################
##############################################################


l = np.array([data['lambda'] for data in total_data]).T

plt.figure(1)
for i, _l in enumerate(l):
    plt.plot(full_range, _l, label='lambda {}'.format(i + 1))

plt.xlabel('distance (A)')
plt.ylabel('unidades')
plt.legend()
#plt.show()


plt.figure(2)
ll = []
for i, lam in enumerate(data['lambda']):
    ll.append(np.array([2*lam * np.sqrt(1 - lam**2) *
                  (data['diabatic_contributions']['W_e_{}'.format(i+1)] + data['diabatic_contributions']['W_h_{}'.format(i+1)])
                  for data in total_data]))

ll = np.array(ll)

for i, lam in enumerate(ll):
    plt.plot(full_range, lam, label='ll {}'.format(i + 1))
# plt.plot(x, ll2, label='ll2')
plt.xlabel('distance (A)')
plt.ylabel('unidades')
plt.legend()
#plt.show()


plt.figure(3)
for i in range(4):
    label = 'E_{}'.format(i+1)
    e = np.array([data['adiabatic_energies'][label] for data in total_data])
    plt.plot(full_range, e, label=label)

plt.xlabel('distance (A)')
plt.ylabel('unidades')
plt.legend()
#plt.show()

plt.figure(4)

d_energies = {}
for e in np.array([data['diabatic_energies'] for data in total_data]):
    for label in e.keys():
        try:
            d_energies[label].append(e[label])
        except:
            d_energies[label] = [e[label]]

print(d_energies)

for label in ['E_LE', 'E_CT']:
    plt.plot(full_range, d_energies[label], label=label)
plt.xlabel('distance (A)')
plt.ylabel('unidades')
plt.legend()
#plt.show()

plt.figure(5)
d_contributions = {}
for e in np.array([data['diabatic_contributions'] for data in total_data]):
    for label in e.keys():
        try:
            d_contributions[label].append(e[label])
        except:
            d_contributions[label] = [e[label]]

print(2)

for label in d_contributions.keys():
    if label.endswith('_1'):
        plt.plot(full_range, d_contributions[label], label=label)
plt.xlabel('distance (A)')
plt.ylabel('unidades')
plt.legend()


plt.figure(6)
for label in d_contributions.keys():
    if label.endswith('_2'):
        plt.plot(full_range, d_contributions[label], label=label)
plt.xlabel('distance (A)')
plt.ylabel('unidades')
plt.legend()


plt.figure(7)

for i, _l in enumerate(l):
    e2 = np.array(_l)**2 * (np.array(d_energies['E_CT']) - np.array(d_energies['E_LE']) +
                            np.array(d_contributions['W_CT_{}'.format(i + 1)]) - np.array(d_contributions['W_DC_{}'.format(i + 1)]))
#       plt.plot(x, e1_, label='E1-')
    plt.plot(full_range, np.array(e2), label='E2-{}'.format(i + 1))

plt.xlabel('distance (A)')
plt.ylabel('unidades')
plt.legend()

plt.show()






import pickle
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages


# very simple order (for 2 states)
def get_order_states(transition_moments, epsilon=0.1):
    order = []
    for i, tm in enumerate(transition_moments):
        if np.linalg.norm(tm) < epsilon:
            order.append(0)
        else:
            order.append(1)
    return order


def correct_order(data_1, data_2, order):
    data_1_o = []
    data_2_o = []

    for i, o in enumerate(order):
        print(o)
        if o[0] == 1 and o[1] == 0:
            data_1_o.append(data_2[i])
            data_2_o.append(data_1[i])
        else:
            data_1_o.append(data_1[i])
            data_2_o.append(data_2[i])

    return data_1_o, data_2_o


def interpolate_data(points, data , y_range, z_range):
    from scipy.interpolate import griddata

    Y, Z = np.meshgrid(y_range, z_range)

    grid_z2 = griddata(points, data, (Y, Z), method='cubic')

    return grid_z2.T.flatten()


def close(p1, p2, p3):
    x = np.array([1, 2, 3], dtype=float)
    y = np.array([p1, p2, p3], dtype=float)

    return np.sum((np.polyval(np.polyfit(x, y, 1), x) - y) ** 2)


def close2(p1_1, p1_2, p2_1, p2_2):
    a1 = np.abs(p1_1 - p1_2) - np.abs(p2_1 - p1_2)
    a2 = np.abs(p2_1 - p2_2) - np.abs(p1_1 - p2_2)

    return False


def correct_data(data_1, data_2, y_range, z_range):
    ny = len(y_range)
    nz = len(z_range)
    data1_rs = np.array(data_1).reshape(len(y_range), len(z_range))
    data2_rs = np.array(data_2).reshape(len(y_range), len(z_range))
    for i in range(ny):
        for j in range(nz):
            if i > 2:
                if close(data1_rs[i, j], data1_rs[i - 1, j], data1_rs[i - 2, j]) > close(data2_rs[i, j], data1_rs[i - 1, j], data1_rs[i - 2, j]):
                    a = data1_rs[i, j]
                    data1_rs[i, j] = data2_rs[i, j]
                    data2_rs[i, j] = a

            if i < ny - 2:
                if close(data1_rs[i, j], data1_rs[i + 1, j], data1_rs[i + 2, j]) > close(data2_rs[i, j],
                                                                                         data1_rs[i + 1, j],
                                                                                         data1_rs[i + 2, j]):
                    a = data1_rs[i, j]
                    data1_rs[i, j] = data2_rs[i, j]
                    data2_rs[i, j] = a

            if j > 2:
                if close(data1_rs[i, j], data1_rs[i, j-1], data1_rs[i, j - 2]) > close(data2_rs[i, j], data1_rs[i, j-1], data1_rs[i, j-2]):
                    a = data1_rs[i, j]
                    data1_rs[i, j] = data2_rs[i, j]
                    data2_rs[i, j] = a


            if j < nz -2:
                if close(data1_rs[i, j], data1_rs[i, j+1], data1_rs[i, j + 2]) > close(data2_rs[i, j], data1_rs[i, j+1], data1_rs[i, j+2]):
                    a = data1_rs[i, j]
                    data1_rs[i, j] = data2_rs[i, j]
                    data2_rs[i, j] = a


    return data1_rs.flatten(), data2_rs.flatten()


def correct_data2(data_1, data_2, points):
    # return data_1, data_2

    d = 0.25
    for _ in range(1):

        for k, p in enumerate(points):
            i, j = p

            if [i + d, j] in points and [i + 2*d, j] in points and i < 0:
                i1 = points.index([i + d, j])
                i2 = points.index([i + 2*d, j])
                if close(data_1[k], data_1[i1], data_1[i2]) > close(data_2[k], data_1[i1], data_1[i2]) and \
                    close2(data_1[k], data_1[i1], data_2[k], data_2[i1]):
                    a = data_1[k]
                    data_1[k] = data_2[k]
                    data_2[k] = a

            if [i - d, j] in points and [i - 2*d, j] in points and i < 0:
                i1 = points.index([i - d, j])
                i2 = points.index([i - 2*d, j])
                if close(data_1[k], data_1[i1], data_1[i2]) > close(data_2[k], data_1[i1], data_1[i2]) and \
                    close2(data_1[k], data_1[i1], data_2[k], data_2[i1]):
                    a = data_1[k]
                    data_1[k] = data_2[k]
                    data_2[k] = a

            if [i, j + d] in points and [i, j + 2*d] in points and j > 0:
                i1 = points.index([i, j + d])
                i2 = points.index([i, j + d])
                if close(data_1[k], data_1[i1], data_1[i2]) > close(data_2[k], data_1[i1], data_1[i2]) and \
                    close2(data_1[k], data_1[i1], data_2[k], data_2[i1]):
                    a = data_1[k]
                    data_1[k] = data_2[k]
                    data_2[k] = a

            if [i, j - d] in points and [i, j - 2*d] in points and j < 0:
                i1 = points.index([i, j - d])
                i2 = points.index([i, j - d])
                if close(data_1[k], data_1[i1], data_1[i2]) > close(data_2[k], data_1[i1], data_1[i2]) and \
                    close2(data_1[k], data_1[i1], data_2[k], data_2[i1]):
                    a = data_1[k]
                    data_1[k] = data_2[k]
                    data_2[k] = a

    return data_1, data_2


def triplot(data1, data2, label1, label2, y_range, z_range, wireframe=False, pdf=None,
            zlabel='Energy [eV]',zrange=np.arange(-0.2, 0.2+0.025, 0.025)):

    # from matplotlib.colors import LinearSegmentedColormap
    # from matplotlib.colors import BoundaryNorm
    # from matplotlib.ticker import MaxNLocator

    cmap = plt.get_cmap('PiYG')

    Y, Z = np.meshgrid(y_range, z_range)

    plt.figure(1)

    plt.title(label1)
    CS = plt.contourf(Y, Z, np.array(data1).reshape(len(y_range), len(z_range)), levels=zrange, cmap=cmap)
    CS2 = plt.contour(CS, levels=CS.levels[::1], colors='black')
    plt.clabel(CS2, inline=1, fontsize=10)
    plt.xlabel('distance Y [Å]')
    plt.ylabel('distance Z [Å]')
    plt.xlim([-4, 4])
    plt.ylim([-4, 4])

    cbar = plt.figure(1).colorbar(CS)

    cbar.ax.set_ylabel(zlabel)

    if pdf is not None:
        pdf.savefig()

    plt.figure(2)

    plt.title(label2)
    CS = plt.contourf(Y, Z, np.array(data2).reshape(len(y_range), len(z_range)), levels=zrange, cmap=cmap)
    CS2 = plt.contour(CS, levels=CS.levels[::1], colors='black')
    plt.clabel(CS2, inline=1, fontsize=10)
    plt.xlabel('distance Y [Å]')
    plt.ylabel('distance Z [Å]')
    plt.xlim([-4, 4])
    plt.ylim([-4, 4])

    cbar = plt.figure(2).colorbar(CS)
    cbar.ax.set_ylabel(zlabel)

    if pdf is not None:
        pdf.savefig()

    fig = plt.figure()

    ax = fig.add_subplot(111, projection='3d')
    if wireframe:
        ax.plot_wireframe(Y, Z, np.array(data1).reshape(len(y_range), len(z_range)), color='b')
        ax.plot_wireframe(Y, Z, np.array(data2).reshape(len(y_range), len(z_range)), color='r')
    else:
        ax.plot_surface(Y, Z, np.array(data1).reshape(len(y_range), len(z_range)), color='b')
        ax.plot_surface(Y, Z, np.array(data2).reshape(len(y_range), len(z_range)), color='r')

    fake2Dline = mpl.lines.Line2D([0], [0], linestyle="none", c='b', marker='o')
    fake2Dline2 = mpl.lines.Line2D([0], [0], linestyle="none", c='r', marker='o')
    ax.legend([fake2Dline, fake2Dline2], [label1, label2], numpoints=1)
    ax.set_xlabel('distance Y [Å]')
    ax.set_ylabel('distance Z [Å]')
    ax.set_zlabel(zlabel)

    if pdf is not None:
        pdf.savefig(fig)

    plt.show()


#############################
folder = '3d_plot/'
#############################

with open('my_data_4.2_2.pkl', 'rb') as input:
    calculation_data = pickle.load(input)
    print('Loaded data from calculation_data.pkl')

do_full = True
interpolate = True

for slide_y in calculation_data['range_y']:
    for slide_z in calculation_data['range_z']:
        print(slide_y, slide_z)
        if '{}_{}'.format(slide_y, slide_z) in calculation_data:
            print(calculation_data['{}_{}'.format(slide_y, slide_z)])
            data_i = calculation_data['{}_{}'.format(slide_y, slide_z)]
            if do_full:

                calculation_data.update({'{}_{}'.format(-slide_y, slide_z): data_i,
                                         '{}_{}'.format(slide_y, -slide_z): data_i,
                                         '{}_{}'.format(-slide_y, -slide_z): data_i})

if do_full:
    y_range = np.unique((-np.array(calculation_data['range_y'])).tolist() + calculation_data['range_y'])
    z_range = np.unique((-np.array(calculation_data['range_z'])).tolist() + calculation_data['range_z'])
else:
    y_range = calculation_data['range_y']
    z_range = calculation_data['range_z']

points = []
total_data = []
full_range = []
i = 0
for slide_y in y_range:
    for slide_z in z_range:
        print('{}_{}'.format(slide_y, slide_z))
        if '{}_{}'.format(slide_y, slide_z) in calculation_data:
            data = calculation_data['{}_{}'.format(slide_y, slide_z)]
            total_data.append(data)
            full_range.append(i)
            points.append([slide_y, slide_z])
            i += 1
# print(total_data[0])

states_orders = [get_order_states(data['transition_moments'][0:2]) for data in total_data]

###################
# Z, Y = np.meshgrid(z_range, y_range)


data_1 = []
for i, e in enumerate([data['diabatic_contributions']['W_DC_1'] for data in total_data]):
    data_1.append(e)

data_2 = []
for i, e in enumerate([data['diabatic_contributions']['W_DC_2'] for data in total_data]):
    data_2.append(e)

data_1, data_2 = correct_order(data_1, data_2, states_orders)

if interpolate:
    data_1 = interpolate_data(points, data_1, y_range, z_range)
    data_2 = interpolate_data(points, data_2, y_range, z_range)


wdc1 = np.array(data_1)
wdc2 = np.array(data_2)

pdf = PdfPages(folder + 'W_DC.pdf')
triplot(data_1, data_2, 'W_DC_1', 'W_DC_2', y_range, z_range, pdf=pdf, wireframe=True)
pdf.close()

#factor = coefficients['S{}_C10'.format(i + 1)] * coefficients['S{}_C01'.format(i + 1)]
#w_dc = factor / np.abs(factor) * diabatic_energies['V_DC']

####################
data_1 = []
data_2 = []
for diab, coeff in [[data['diabatic_energies'], data['coefficients']] for data in total_data]:
    #factor = coeff['S1_C10'] * coeff['S1_CCA'] + coeff['S1_C01'] * coeff['S1_CAC']
    #data_1.append(diab['V_e'] * factor)

    factor = coeff['S1_CCA'] * coeff['S1_CAC']
    data_1.append(diab['V_CT'] * np.sign(factor))

    factor = coeff['S2_CCA'] * coeff['S2_CAC']
    data_2.append(diab['V_CT'] * np.sign(factor))

    #factor = coeff['S1_C10'] * coeff['S1_C01']
    #data_1.append(diab['V_DC'] * factor)

data_1, data_2 = correct_order(data_1, data_2, states_orders)

if interpolate:
    data_1 = interpolate_data(points, data_1, y_range, z_range)
    data_2 = interpolate_data(points, data_2, y_range, z_range)

wct1 = np.array(data_1)
wct2 = np.array(data_2)

pdf = PdfPages(folder + 'W_CT.pdf')
triplot(data_1, data_2, 'W_CT_1', 'W_CT_2', y_range, z_range, pdf=pdf, wireframe=False)
pdf.close()


#################### W_CT
data_1 = []
for i, e in enumerate([data['diabatic_contributions']['W_CT_1'] for data in total_data]):
    #print(states_orders[i][0])

    #factor = total_data[i]['coefficients']['S{}_CCA'.format(0 + 1)] * total_data[i]['coefficients']['S{}_CAC'.format(0 + 1)]
    #w_ct = factor / np.abs(factor) * total_data[i]['diabatic_energies']['V_CT']
    #data_1.append(w_ct)

    data_1.append(e)

data_2 = []
for i, e in enumerate([data['diabatic_contributions']['W_CT_2'] for data in total_data]):
    #factor = total_data[i]['coefficients']['S{}_CCA'.format(1 + 1)] * total_data[i]['coefficients']['S{}_CAC'.format(1 + 1)]
    #w_ct = factor / np.abs(factor) * total_data[i]['diabatic_energies']['V_CT']
    #data_2.append(w_ct)

    data_2.append(e)

#data_1, data_2 = correct_order(data_1, data_2, states_orders)

if interpolate:
    data_1 = interpolate_data(points, data_1, y_range, z_range)
    data_2 = interpolate_data(points, data_2, y_range, z_range)

#pdf = PdfPages(folder + 'W_CT.pdf')
#triplot(data_1, data_2, 'W_CT_1', 'W_CT_2', y_range, z_range, pdf=pdf)
#pdf.close()


#####################

data_1 = []
for i, e in enumerate([data['diabatic_contributions']['W_e_1'] for data in total_data]):
    data_1.append(e)

data_2 = []
for i, e in enumerate([data['diabatic_contributions']['W_e_2'] for data in total_data]):
    data_2.append(e)


# correct sign:
for i, lc in enumerate(zip(data_1, data_2)):
    la, lb = lc
    if la < 0:
        data_1[i] *= -1
    if lb > 0:
        data_2[i] *= -1

#data_1, data_2 = correct_order(data_1, data_2, states_orders)

if interpolate:
    data_1 = interpolate_data(points, data_1, y_range, z_range)
    data_2 = interpolate_data(points, data_2, y_range, z_range)

pdf = PdfPages(folder + 'W_e.pdf')
triplot(data_1, data_2, 'W_e_1', 'W_e_2', y_range, z_range, pdf=pdf)
pdf.close()

we1 = data_1
we2 = data_2

#exit()

############
data_1 = []
data_2 = []
for diab, coeff in [[data['diabatic_energies'], data['coefficients']] for data in total_data]:

    factor = coeff['S1_C10'] * coeff['S1_CCA'] + coeff['S1_C01'] * coeff['S1_CAC']
    # factor = coeff['S1_C10'] * coeff['S1_CAC'] + coeff['S1_C01'] * coeff['S1_CCA']

    data_1.append(diab['V_e'] * np.sign(factor))

    #factor = coeff['S1_C10'] * coeff['S1_CAC'] + coeff['S1_C01'] * coeff['S1_CCA']
    #data_1.append(diab['V_h'] * factor)

    factor = coeff['S2_C10'] * coeff['S2_CCA'] + coeff['S2_C01'] * coeff['S2_CAC']
    # factor = coeff['S2_C10'] * coeff['S2_CAC'] + coeff['S2_C01'] * coeff['S2_CCA']

    data_2.append(diab['V_e'] * np.sign(factor))

    #factor = coeff['S2_C10'] * coeff['S2_CAC'] + coeff['S2_C01'] * coeff['S2_CCA']
    #data_2.append(diab['V_h'] * factor)


#data_1, data_2 = correct_order(data_1, data_2, states_orders)

# correct sign:
#for i, lc in enumerate(zip(data_1, data_2)):
#    la, lb = lc
#    if la > 0:
#        data_1[i] *= -1
#    if lb < 0:
#        data_2[i] *= -1

if interpolate:
    data_1 = interpolate_data(points, data_1, y_range, z_range)
    data_2 = interpolate_data(points, data_2, y_range, z_range)


pdf = PdfPages(folder + 'W_e.pdf')
triplot(data_1, data_2, 'W_e_1', 'W_e_2', y_range, z_range, pdf=pdf)
pdf.close()


########################

data_1 = []
for i, e in enumerate([data['diabatic_contributions']['W_h_1'] for data in total_data]):
    data_1.append(e)

data_2 = []
for i, e in enumerate([data['diabatic_contributions']['W_h_2'] for data in total_data]):
    data_2.append(e)

#data_1, data_2 = correct_order(data_1, data_2, states_orders)

# correct order:
for i, lc in enumerate(zip(data_1, data_2)):
    l1, l2 = lc
    if l1 > 0:
        data_1[i] *= -1
    if l2 < 0:
        data_2[i] *= -1

if interpolate:
    data_1 = interpolate_data(points, data_1, y_range, z_range)
    data_2 = interpolate_data(points, data_2, y_range, z_range)

with PdfPages(folder + 'W_h.pdf') as pdf:
    triplot(data_1, data_2, 'W_h_1', 'W_h_2', y_range, z_range, pdf=pdf)

wh1 = data_1
wh2 = data_2

##
data_1 = []
data_2 = []
for diab, coeff in [[data['diabatic_energies'], data['coefficients']] for data in total_data]:

    #factor = coeff['S1_C10'] * coeff['S1_CAC'] + coeff['S1_C01'] * coeff['S1_CCA']

    data_1.append(diab['V_h'] * np.sign(factor))

    #factor = coeff['S1_C10'] * coeff['S1_CAC'] + coeff['S1_C01'] * coeff['S1_CCA']
    #data_1.append(diab['V_h'] * factor)

    factor = coeff['S2_C10'] * coeff['S2_CAC'] + coeff['S2_C01'] * coeff['S2_CCA']

    data_2.append(diab['V_h'] * np.sign(factor))

    #factor = coeff['S2_C10'] * coeff['S2_CAC'] + coeff['S2_C01'] * coeff['S2_CCA']
    #data_2.append(diab['V_h'] * factor)


#data_1, data_2 = correct_order(data_1, data_2, states_orders)

# correct sign:
#for i, lc in enumerate(zip(data_1, data_2)):
#    la, lb = lc
#    if la > 0:
#        data_1[i] *= -1
#    if lb < 0:
#        data_2[i] *= -1

if interpolate:
    data_1 = interpolate_data(points, data_1, y_range, z_range)
    data_2 = interpolate_data(points, data_2, y_range, z_range)

with PdfPages(folder + 'W_h.pdf') as pdf:
    triplot(data_1, data_2, 'W_h_1', 'W_h_2', y_range, z_range, pdf=pdf)


#exit()
#######################


data_1 = []
for i, e in enumerate(np.array([data['lambda'] for data in total_data]).T[0]):
    data_1.append(e)

data_2 = []
for i, e in enumerate(np.array([data['lambda'] for data in total_data]).T[1]):
    data_2.append(e)

#data_1, data_2 = correct_data2(data_1, data_2, points)
data_1, data_2 = correct_order(data_1, data_2, states_orders)


if interpolate:
    data_1 = interpolate_data(points, data_1, y_range, z_range)
    data_2 = interpolate_data(points, data_2, y_range, z_range)


#data_1, data_2 = correct_data(data_1, data_2, y_range, z_range)

with PdfPages(folder + 'lambda.pdf') as pdf:
    triplot(data_1, data_2, 'lambda 1', 'lambda 2', y_range, z_range, wireframe=True, pdf=pdf, zlabel='')

l1 = np.array(data_1)
l2 = np.array(data_2)


##########################

d_energies = {}
for e in np.array([data['diabatic_energies'] for data in total_data]):
    for label in e.keys():
        try:
            d_energies[label].append(e[label])
        except:
            d_energies[label] = [e[label]]


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

with PdfPages(folder + 'E1.pdf') as pdf:
    triplot(data_1, data_2, 'E1-1', 'E1-2', y_range, z_range, wireframe=True, pdf=pdf)

#exit()

############

data_1 = []
for i, e in enumerate([data['diabatic_energies']['E_LE'] for data in total_data]):
    data_1.append(e)

data_2 = []
for i, e in enumerate([data['diabatic_energies']['E_CT'] for data in total_data]):
    data_2.append(e)

#data_1, data_2 = correct_order(data_1, data_2, states_orders)

# correct order:
#for i, lc in enumerate(zip(data_1, data_2)):
#    la, lb = lc
#    if la > 0:
#        data_1[i] *= -1
#    if lb < 0:
#        data_2[i] *= -1

#data_1, data_2 = correct_data(data_1, data_2, y_range, z_range)

if interpolate:
    data_1 = interpolate_data(points, data_1, y_range, z_range)
    data_2 = interpolate_data(points, data_2, y_range, z_range)


e_le = np.array(data_1)
e_ct = np.array(data_2)

with PdfPages(folder + 'diabatic_energies.pdf') as pdf:
    triplot(data_1, data_2, 'E_LE', 'E_CT', y_range, z_range, pdf=pdf, zrange=None, wireframe=False)


###
#l = np.array([data['lambda'] for data in total_data]).T

#data_1 = np.array(l[0])**2 * (np.array(d_energies['E_CT']) - np.array(d_energies['E_LE']) +
#                  np.array(d_contributions['W_CT_1']) - np.array(d_contributions['W_DC_1']))

#data_2 = np.array(l[1])**2 * (np.array(d_energies['E_CT']) - np.array(d_energies['E_LE']) +
#                  np.array(d_contributions['W_CT_2']) - np.array(d_contributions['W_DC_2']))

#if interpolate:
#    data_1 = interpolate_data(points, data_1, y_range, z_range)
#    data_2 = interpolate_data(points, data_2, y_range, z_range)

data_1 = l1**2 * (e_ct - e_le + wct1 - wdc1)
data_2 = l2**2 * (e_ct - e_le + wct2 - wdc2)

#data_1, data_2 = correct_data(data_1, data_2, y_range, z_range)

with PdfPages(folder + 'E2.pdf') as pdf:
    triplot(data_1, data_2, 'E2-1', 'E2-2', y_range, z_range, wireframe=True, pdf=pdf)



data_1 = []
for i, e in enumerate([data['adiabatic_energies']['E_1'] for data in total_data]):
    data_1.append(e)

data_2 = []
for i, e in enumerate([data['adiabatic_energies']['E_2'] for data in total_data]):
    data_2.append(e)

data_1, data_2 = correct_order(data_1, data_2, states_orders)

if interpolate:
    data_1 = interpolate_data(points, data_1, y_range, z_range)
    data_2 = interpolate_data(points, data_2, y_range, z_range)

with PdfPages(folder + 'adiabatic_energies.pdf') as pdf:
    triplot(data_1, data_2, 'E_1', 'E_2', y_range, z_range, pdf=pdf, zrange=None)







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




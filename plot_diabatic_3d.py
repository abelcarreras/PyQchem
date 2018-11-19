import pickle
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import argparse

# Argument parser
parser = argparse.ArgumentParser(description='Plot 3D data')
parser.add_argument('filename', metavar='filename', type=str,
                    help='filename for input')

parser.add_argument('--output_folder', metavar='distance', type=str, default='3d_plot/',
                    help='folder to store PDF plots')

parser.add_argument('--show_plots', action='store_true',
                   help='show plots')

parser.add_argument('--full', action='store_true',
                   help='apply periodicity')

parser.add_argument('--interpolate', action='store_true',
                   help='interpolate missing points')


args = parser.parse_args()


# very simple order (for 2 states)
def get_order_states2(transition_moments, epsilon=0.1):
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

        if o[0] == 1 and o[1] == 0:
            data_1_o.append(data_2[i])
            data_2_o.append(data_1[i])
        else:
            data_1_o.append(data_1[i])
            data_2_o.append(data_2[i])

    return data_1_o, data_2_o


def get_order_states(states, epsilon=0.1):

    order = np.argsort([state['total energy'] for state in states]).tolist()
    for i, tm in enumerate([state['transition moment'] for state in states]):
        if np.linalg.norm(tm) > epsilon:
            index = order.pop(i)
            order = [index] + order

    return order


def correct_order2(list, order):

    alist = np.array(list)
    orddered_list = alist[order]

    return orddered_list.tolist()


def interpolate_data(points, data , y_range, z_range):
    from scipy.interpolate import griddata

    Y, Z = np.meshgrid(y_range, z_range)

    grid_z2 = griddata(points, data, (Y, Z), method='cubic')

    return grid_z2.T.flatten()


def triplot(data1, data2, label1, label2, y_range, z_range, wireframe=False, pdf=None,
            zlabel='Energy [eV]', zlevels=np.arange(-0.2, 0.2 + 0.025, 0.025), show_plot=True):

    # from matplotlib.colors import LinearSegmentedColormap
    # from matplotlib.colors import BoundaryNorm
    # from matplotlib.ticker import MaxNLocator

    cmap = plt.get_cmap('PiYG')

    Y, Z = np.meshgrid(z_range, y_range)

    plt.figure(1)

    plt.title(label1)
    CS = plt.contourf(Y, Z, np.array(data1).reshape(len(y_range), len(z_range)), levels=zlevels, cmap=cmap)
    CS2 = plt.contour(CS, levels=CS.levels[::1], colors='black')
    plt.clabel(CS2, inline=1, fontsize=10)
    plt.xlabel('distance Z [Å]')
    plt.ylabel('distance Y [Å]')
    plt.xlim([-4, 4])
    plt.ylim([-4, 4])

    cbar = plt.figure(1).colorbar(CS)

    cbar.ax.set_ylabel(zlabel)

    if pdf is not None:
        pdf.savefig()

    plt.figure(2)

    plt.title(label2)
    CS = plt.contourf(Y, Z, np.array(data2).reshape(len(y_range), len(z_range)), levels=zlevels, cmap=cmap)
    CS2 = plt.contour(CS, levels=CS.levels[::1], colors='black')
    plt.clabel(CS2, inline=1, fontsize=10)
    plt.xlabel('distance Z [Å]')
    plt.ylabel('distance Y [Å]')
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
    ax.set_xlabel('distance Z [Å]')
    ax.set_ylabel('distance Y [Å]')
    ax.set_zlabel(zlabel)

    if pdf is not None:
        pdf.savefig(fig)

    if show_plot:
        plt.show()


def biplot2(data1, data2, label1, label2, y_range, z_range, pdf=None,
            zlabel='Energy [eV]', zrange=(-0.2, 0.2), zp=0.025, show_plot=True):

    from scipy import interpolate

    Y, Z = np.meshgrid(z_range, y_range)

    print('----------------')
    print(Y.shape)
    print(len(Y), len(data1))

    #f1 = interpolate.interp2d(Y, Z, data1, kind='cubic')
    f1 = interpolate.interp2d(Y, Z, data1, kind='linear')
    f2 = interpolate.interp2d(Y, Z, data2, kind='linear')

    x = np.arange(-4, 4, zp)

    plt.figure(3)

    plt.xlim([-4, 4])
    plt.ylim([zrange[0], zrange[1]])
    plt.xlabel('distance Z [Å]')
    plt.ylabel(zlabel)

    plt.plot(x, f1(x, 0), label=label1)
    plt.plot(x, f2(x, 0), label=label2)
    plt.legend()

    if pdf is not None:
        pdf.savefig()

    if show_plot:
        plt.show()


def biplot(data1, data2, label1, label2, y_range, z_range, pdf=None,
            zlabel='Energy [eV]', zrange=(-0.2, 0.2), title=None, show_plot=True):

    data1 = np.array(data1).reshape([len(y_range), len(z_range)])[len(z_range)//2]
    data2 = np.array(data2).reshape([len(y_range), len(z_range)])[len(z_range)//2]

    plt.title(title)
    plt.xlim([-4, 4])
    plt.ylim([zrange[0], zrange[1]])
    plt.xlabel('distance Z [Å]')
    plt.ylabel(zlabel)

    plt.plot(z_range, data1, label=label1)
    plt.plot(z_range, data2, label=label2)
    plt.legend()

    if pdf is not None:
        pdf.savefig()

    if show_plot:
        plt.show()


#############################
folder = args.output_folder
#############################

with open(args.filename, 'rb') as input:
    calculation_data = pickle.load(input)
    print('Loaded data from: {}'.format(args.filename))

do_full = args.full
interpolate = args.interpolate

for slide_y in calculation_data['range_y']:
    for slide_z in calculation_data['range_z']:
        print(slide_y, slide_z)
        if '{}_{}'.format(slide_y, slide_z) in calculation_data:
            print(calculation_data['{}_{}'.format(slide_y, slide_z)])
            data_i = calculation_data['{}_{}'.format(slide_y, slide_z)]
            print(data_i['states_info'])
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

states_orders = [get_order_states(data['states_info'])[:2] for data in total_data]

#exit()
###################

data_1 = []
data_2 = []
for diab, coeff in [[data['diabatic_energies'], data['coefficients']] for data in total_data]:

    factor = coeff['S_01'][1] * coeff['S_10'][1]
    data_1.append(diab['V_DC'] * np.sign(factor))
    # data_1.append(diab['E_LE'][0])

    factor = coeff['S_01'][1] * coeff['S_10'][1]
    data_2.append(diab['V_DC'] * np.sign(factor))
    # data_2.append(diab['E_LE'][0])

print(data_1)
print(states_orders)
data_1, data_2 = correct_order(data_1, data_2, states_orders)

#if interpolate:
#    data_1 = interpolate_data(points, data_1, y_range, z_range)
#    data_2 = interpolate_data(points, data_2, y_range, z_range)

wdc1 = np.array(data_1)
wdc2 = np.array(data_2)


with PdfPages(folder + 'W_DC.pdf') as pdf:
    triplot(data_1, data_2, 'W_DC_1', 'W_DC_2', y_range, z_range, pdf=pdf, wireframe=True, show_plot=args.show_plots)
    biplot(data_1, data_2, 'W_DC_1', 'W_DC_2', y_range, z_range, show_plot=args.show_plots)

exit()
# ##################### W_DC
if True:
    data_1 = []
    for i, e in enumerate([data['diabatic_contributions']['W_DC'][0] for data in total_data]):
        data_1.append(e)

    data_2 = []
    for i, e in enumerate([data['diabatic_contributions']['W_DC'][0] for data in total_data]):
        data_2.append(e)

    data_1, data_2 = correct_order(data_1, data_2, states_orders)

    if interpolate:
        data_1 = interpolate_data(points, data_1, y_range, z_range)
        data_2 = interpolate_data(points, data_2, y_range, z_range)


    wdc1 = np.array(data_1)
    wdc2 = np.array(data_2)

    #pdf = PdfPages(folder + 'W_DC.pdf')
    triplot(data_1, data_2, 'W_DC_1', 'W_DC_2', y_range, z_range, wireframe=True, show_plot=args.show_plots)
    #pdf.close()


    #factor = coefficients['S{}_C10'.format(i + 1)] * coefficients['S{}_C01'.format(i + 1)]
    #w_dc = factor / np.abs(factor) * diabatic_energies['V_DC']


####################

data_1 = []
data_2 = []
for diab, coeff in [[data['diabatic_energies'], data['coefficients']] for data in total_data]:

    factor = coeff['S1_CCA'] * coeff['S1_CAC']
    data_1.append(diab['V_CT'] * np.sign(factor))

    factor = coeff['S2_CCA'] * coeff['S2_CAC']
    data_2.append(diab['V_CT'] * np.sign(factor))



data_1, data_2 = correct_order(data_1, data_2, states_orders)

if interpolate:
    data_1 = interpolate_data(points, data_1, y_range, z_range)
    data_2 = interpolate_data(points, data_2, y_range, z_range)

wct1 = np.array(data_1)
wct2 = np.array(data_2)

triplot(data_1, data_2, 'W_CT_1', 'W_CT_2', y_range, z_range, wireframe=False, show_plot=args.show_plots)


#################### W_CT
if True:
    data_1 = []
    for i, e in enumerate([data['diabatic_contributions']['W_CT_1'] for data in total_data]):
        data_1.append(e)

    data_2 = []
    for i, e in enumerate([data['diabatic_contributions']['W_CT_2'] for data in total_data]):
        data_2.append(e)

    data_1, data_2 = correct_order(data_1, data_2, states_orders)

    if interpolate:
        data_1 = interpolate_data(points, data_1, y_range, z_range)
        data_2 = interpolate_data(points, data_2, y_range, z_range)

    with PdfPages(folder + 'W_CT.pdf') as pdf:
        triplot(data_1, data_2, 'W_CT_1', 'W_CT_2', y_range, z_range, pdf=pdf, show_plot=args.show_plots)
        biplot(data_1, data_2, 'W_CT_1', 'W_CT_2', y_range, z_range, show_plot=args.show_plots)

    wct1 = np.array(data_1)
    wct2 = np.array(data_2)


##################### W_e
if False:
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


    with PdfPages(folder + 'W_e.pdf') as pdf:
        triplot(data_1, data_2, 'W_e_1', 'W_e_2', y_range, z_range, pdf=pdf, show_plot=args.show_plots)
        biplot(data_1, data_2, 'W_e_1', 'W_e_2', y_range, z_range, pdf=pdf, show_plot=args.show_plots)

    we1 = data_1
    we2 = data_2

    #exit()

############
data_1 = []
data_2 = []
for diab, coeff in [[data['diabatic_energies'], data['coefficients']] for data in total_data]:

    #factor = coeff['S1_C10'] * coeff['S1_CAC'] + coeff['S1_C01'] * coeff['S1_CCA']

    factor = coeff['S1_C01'] * coeff['S1_CAC']
    factor2 = coeff['S1_C10'] * coeff['S1_CCA']

    #data_1.append(diab['V_e'] * np.sign(factor))
    data_1.append(diab['V_e'] * np.sign(factor))

    #factor = coeff['S2_C10'] * coeff['S2_CCA'] + coeff['S2_C01'] * coeff['S2_CAC']

    factor = coeff['S2_C01'] * coeff['S2_CAC']
    factor2 = coeff['S2_C10'] * coeff['S2_CCA']

    # data_2.append(diab['V_e'] * np.sign(factor))
    data_2.append(diab['V_e'] * np.sign(factor))


#data_1, data_2 = correct_order(data_1, data_2, states_orders)

if interpolate:
    data_1 = interpolate_data(points, data_1, y_range, z_range)
    data_2 = interpolate_data(points, data_2, y_range, z_range)


with PdfPages(folder + 'W_e.pdf') as pdf:
    triplot(data_1, data_2, 'W_e_1', 'W_e_2', y_range, z_range, pdf=pdf, show_plot=args.show_plots)
    biplot(data_1, data_2, 'W_e_1', 'W_e_2', y_range, z_range, pdf=pdf, show_plot=args.show_plots)

we1 = data_1
we2 = data_2

#exit()
######################## W_h
if False:
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
        triplot(data_1, data_2, 'W_h_1', 'W_h_2', y_range, z_range, pdf=pdf, show_plot=args.show_plots)
        biplot(data_1, data_2, 'W_h_1', 'W_h_2', y_range, z_range, pdf=pdf, show_plot=args.show_plots)

    wh1 = data_1
    wh2 = data_2

##
data_1 = []
data_2 = []
for diab, coeff in [[data['diabatic_energies'], data['coefficients']] for data in total_data]:

    factor1 = coeff['S1_C10'] * coeff['S1_CAC']
    factor2 = coeff['S1_C01'] * coeff['S1_CCA']
    norm = np.abs(factor1) + np.abs(factor2)
    data_1.append(diab['V_h'] * np.sign(factor1))
    #data_1.append((diab['V_h'] * factor1 + diab['V_h_2'] * factor2))
    #data_1.append((diab['V_h'] * np.sign(factor1) + diab['V_h_2'] * np.sign(factor2)))

    #data_1.append(norm)


    factor1 = coeff['S2_C10'] * coeff['S2_CAC']
    factor2 = coeff['S2_C01'] * coeff['S2_CCA']
    norm = np.abs(factor1) + np.abs(factor2)

    data_2.append(diab['V_h'] * np.sign(factor1))
    #data_2.append((diab['V_h'] * np.sign(factor1) + diab['V_h_2'] * np.sign(factor2)))
    #data_2.append(norm)

    #factor = coeff['S2_C10'] * coeff['S2_CAC'] + coeff['S2_C01'] * coeff['S2_CCA']
    #data_2.append(diab['V_h'] * factor)


data_1, data_2 = correct_order(data_1, data_2, states_orders)

if interpolate:
    data_1 = interpolate_data(points, data_1, y_range, z_range)
    data_2 = interpolate_data(points, data_2, y_range, z_range)

with PdfPages(folder + 'W_h.pdf') as pdf:
    triplot(data_1, data_2, 'W_h_1', 'W_h_2', y_range, z_range, pdf=pdf, show_plot=args.show_plots)
    biplot(data_1, data_2, 'W_h_1', 'W_h_2', y_range, z_range, pdf=pdf, show_plot=args.show_plots)

wh1 = data_1
wh2 = data_2

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
    triplot(data_1, data_2, 'lambda 1', 'lambda 2', y_range, z_range, wireframe=True, pdf=pdf, zlabel='', show_plot=args.show_plots)
    biplot(data_1, data_2, 'lambda 1', 'lambda 2', y_range, z_range, pdf=pdf, show_plot=args.show_plots)

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
    triplot(data_1, data_2, 'E1-1', 'E1-2', y_range, z_range, wireframe=True, pdf=pdf, show_plot=args.show_plots)
    biplot(data_1, data_2, 'E1-1', 'E1-2', y_range, z_range, pdf=pdf, show_plot=args.show_plots)

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
    triplot(data_1, data_2, 'E_LE', 'E_CT', y_range, z_range, pdf=pdf, zlevels=None, wireframe=False, show_plot=args.show_plots)
    biplot(data_1, data_2, 'E_LE', 'E_CT', y_range, z_range, show_plot=args.show_plots)

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
    triplot(data_1, data_2, 'E2-1', 'E2-2', y_range, z_range, wireframe=True, pdf=pdf, show_plot=args.show_plots)
    biplot(data_1, data_2, 'E2-1', 'E2-2', y_range, z_range, show_plot=args.show_plots)


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
    triplot(data_1, data_2, 'E_1', 'E_2', y_range, z_range, pdf=pdf, zlevels=None, show_plot=args.show_plots)
    biplot(data_1, data_2, 'E_1', 'E_2', y_range, z_range, pdf=pdf, show_plot=args.show_plots)





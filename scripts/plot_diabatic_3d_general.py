import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import axes3d  # necessary !!
from pyqchem.order_states import get_order_states_list, correct_order_list

mpl.rcParams.update({'font.size': 13})

import argparse

# Argument parser
parser = argparse.ArgumentParser(description='Plot 3D data')
parser.add_argument('filename', metavar='filename', type=str,
                    help='filename for input')

parser.add_argument('--output_folder', metavar='distance', type=str, default='3d_plot/',
                    help='folder to store PDF plots')

parser.add_argument('--show_plots', action='store_true',
                   help='show plots while running')

parser.add_argument('--full', action='store_true',
                   help='apply periodicity (centered at 0,0)')

parser.add_argument('--interpolate', action='store_true',
                   help='interpolate missing points')


args = parser.parse_args()


def interpolate_data(points, data , y_range, z_range):
    from scipy.interpolate import griddata

    Y, Z = np.meshgrid(y_range, z_range)

    grid_z2 = griddata(points, data, (Y, Z), method='cubic')

    return grid_z2.T.flatten()



default_range = 0.4

def triplot2(datas, labels, y_range, z_range, wireframe=False, pdf=None, clabels=False,
            zlabel='Energy [eV]', zlevels=np.arange(-default_range, default_range + 0.025, 0.025),
             show_plot=True, colors=('b', 'r', 'y', 'm')):

    # from matplotlib.colors import LinearSegmentedColormap
    # from matplotlib.colors import BoundaryNorm
    # from matplotlib.ticker import MaxNLocator

    cmap = plt.get_cmap('PiYG')

    Y, Z = np.meshgrid(z_range, y_range)

    for iplot, data in enumerate(datas):
        plt.figure(iplot+1)

        plt.title(labels[iplot])
        #zlevels = np.arange(-default_range, default_range + 0.025, 0.025)
        CS = plt.contourf(Y, Z, np.array(data).reshape(len(y_range), len(z_range)), levels=zlevels, cmap=cmap)
        # print(CS.levels[::1])

        # force symmetry
        #dist = CS.levels[1] - CS.levels[0]
        #if np.abs(CS.levels[0]) < (CS.levels[-1]):
        #    zlevels = np.arange(-CS.levels[-1], CS.levels[-1], dist)
        #else:
        #    zlevels = np.arange(CS.levels[0], -CS.levels[0], dist)
        #CS = plt.contourf(Y, Z, np.array(data).reshape(len(y_range), len(z_range)), levels=zlevels, cmap=cmap)

        CS2 = plt.contour(CS, levels=CS.levels[::1], colors='black', linestyles='solid')
        if clabels:
            plt.clabel(CS2, inline=1, fontsize=10)

        plt.xlabel('X [Å]')
        plt.ylabel('Y [Å]')
        plt.xlim([-4, 4])
        plt.ylim([-4, 4])

        cbar = plt.figure(iplot+1).colorbar(CS)

        cbar.ax.set_ylabel(zlabel)

        if pdf is not None:
            pdf.savefig()

        if not show_plot:
            plt.close()

    fig = plt.figure()

    ax = fig.add_subplot(111, projection='3d')

    fake2Dline = []
    for iplot, data in enumerate(datas):
        if wireframe:
            ax.plot_wireframe(Y, Z, np.array(data).reshape(len(y_range), len(z_range)), color=colors[iplot])
        else:
            ax.plot_surface(Y, Z, np.array(data).reshape(len(y_range), len(z_range)), color=colors[iplot])
        fake2Dline.append(mpl.lines.Line2D([0], [0], linestyle="none", c=colors[iplot], marker='o'))
    ax.legend(fake2Dline, labels, numpoints=1)
    ax.set_xlabel('X [Å]')
    ax.set_ylabel('Y [Å]')
    ax.set_zlabel(zlabel)

    if pdf is not None:
        pdf.savefig(fig)

    if show_plot:
        plt.show()
    else:
        plt.close()



def triplot3(datas, labels, y_range, z_range, wireframe=False, pdf=None, clabels=False,
            zlabel='Energy [eV]', zlevels=np.arange(-default_range, default_range + 0.025, 0.025),
             show_plot=True, colors=('b', 'r', 'y', 'm')):

    # from matplotlib.colors import LinearSegmentedColormap
    # from matplotlib.colors import BoundaryNorm
    # from matplotlib.ticker import MaxNLocator

    cmap = plt.get_cmap('PiYG')

    Y, Z = np.meshgrid(z_range, y_range)

    for iplot, data in enumerate(datas):
        plt.figure(iplot+1)

        plt.title(labels[iplot])
        #zlevels = np.arange(-default_range, default_range + 0.025, 0.025)
        CS = plt.contourf(Y, Z, np.array(data).reshape(len(y_range), len(z_range)), levels=zlevels, cmap=cmap)
        # print(CS.levels[::1])

        # force symmetry
        #dist = CS.levels[1] - CS.levels[0]
        #if np.abs(CS.levels[0]) < (CS.levels[-1]):
        #    zlevels = np.arange(-CS.levels[-1], CS.levels[-1], dist)
        #else:
        #    zlevels = np.arange(CS.levels[0], -CS.levels[0], dist)
        #CS = plt.contourf(Y, Z, np.array(data).reshape(len(y_range), len(z_range)), levels=zlevels, cmap=cmap)

        CS2 = plt.contour(CS, levels=CS.levels[::1], colors='black', linestyles='solid')
        if clabels:
            plt.clabel(CS2, inline=1, fontsize=10)

        plt.xlabel('X [Å]')
        plt.ylabel('Y [Å]')
        plt.xlim([-4, 4])
        plt.ylim([-4, 4])

        cbar = plt.figure(iplot+1).colorbar(CS)

        cbar.ax.set_ylabel(zlabel)

        if pdf is not None:
            pdf.savefig()

        if not show_plot:
            plt.close()

    fig = plt.figure()

    ax = fig.add_subplot(111, projection='3d')

    fake2Dline = []
    if wireframe:
        ax.plot_wireframe(Y, Z, np.array(datas[0]).reshape(len(y_range), len(z_range)), color=colors[0])
    else:
        ax.plot_surface(Y, Z, np.array(datas[0]).reshape(len(y_range), len(z_range)), color=colors[0])
    fake2Dline.append(mpl.lines.Line2D([0], [0], linestyle="none", c=colors[0], marker='o'))
    ax.set_xlabel('X [Å]')
    ax.set_ylabel('Y [Å]')
    ax.set_zlabel(zlabel)


    fig2 = plt.figure()

    ax2 = fig2.add_subplot(111, projection='3d')

    fake2Dline = []
    if wireframe:
        ax2.plot_wireframe(Y, Z, np.array(datas[1]).reshape(len(y_range), len(z_range)), color=colors[1])
    else:
        ax2.plot_surface(Y, Z, np.array(datas[1]).reshape(len(y_range), len(z_range)), color=colors[1])
    fake2Dline.append(mpl.lines.Line2D([0], [0], linestyle="none", c=colors[1], marker='o'))
    ax2.set_xlabel('X [Å]')
    ax2.set_ylabel('Y [Å]')
    ax2.set_zlabel(zlabel)



    if pdf is not None:
        pdf.savefig(fig)
        pdf.savefig(fig2)

    if show_plot:
        plt.show()
    else:
        plt.close()

def biplot_interpolated(data1, data2, label1, label2, y_range, z_range, pdf=None,
                        zlabel='Energy [eV]', zrange=(-default_range, default_range), zp=0.025, show_plot=True):

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
    plt.xlabel('X [Å]')
    plt.ylabel(zlabel)

    plt.plot(x, f1(x, 0), label=label1)
    plt.plot(x, f2(x, 0), label=label2)
    plt.legend()

    if pdf is not None:
        pdf.savefig()

    if show_plot:
        plt.show()
    else:
        plt.close()


def biplot(data1, data2, label1, label2, y_range, z_range, pdf=None,
           zlabel='Energy [eV]', zrange=(-default_range, default_range), title=None, show_plot=True,
           direction=0, colors=[None, None]):

    plt.figure(3)
    plt.gcf().subplots_adjust(left=0.15)

    if direction == 0:
        data1 = np.array(data1).reshape([len(y_range), len(z_range)])[len(z_range)//2]
        data2 = np.array(data2).reshape([len(y_range), len(z_range)])[len(z_range)//2]
        plt.xlabel('X [Å]')

    if direction == 1:
        data1 = np.array(data1).reshape([len(y_range), len(z_range)])[:,len(y_range)//2]
        data2 = np.array(data2).reshape([len(y_range), len(z_range)])[:,len(y_range)//2]
        plt.xlabel('Y [Å]')

    plt.title(title)
    plt.xlim([-4, 4])
    if zrange is not None:
        plt.ylim([zrange[0], zrange[1]])
    plt.ylabel(zlabel)

    plt.plot(z_range, data1, label=label1, color=colors[0])
    plt.plot(z_range, data2, label=label2, color=colors[1])
    plt.legend()

    if pdf is not None:
        pdf.savefig()

    if show_plot:
        plt.show()
    else:
        plt.close()


def multibiplot(data, labels, y_range, z_range, pdf=None,
                zlabel='Energy [eV]', zrange=(-default_range, default_range), title=None, show_plot=True,
                direction=0):

    plt.figure(3)

    if direction == 0:
        for i, dat in enumerate(data):
            data[i] = np.array(dat).reshape([len(y_range), len(z_range)])[len(z_range)//2]
            #data[i] = np.array(dat).reshape([len(y_range), len(z_range)])[len(z_range)//2]
            plt.xlabel('X [Å]')

    if direction == 1:
        for i, dat in enumerate(data):

            data[i] = np.array(dat).reshape([len(y_range), len(z_range)])[:,len(y_range)//2]
            # data2 = np.array(data2).reshape([len(y_range), len(z_range)])[:,len(y_range)//2]
            plt.xlabel('Y [Å]')

    plt.title(title)
    plt.xlim([-4, 4])
    if zrange is not None:
        plt.ylim([zrange[0], zrange[1]])
    plt.ylabel(zlabel)
    for dat, l in zip(data, labels):
        plt.plot(z_range, dat, label=l)
        # plt.plot(z_range, data2, label=label2)
    plt.legend()

    if pdf is not None:
        pdf.savefig()

    if show_plot:
        plt.show()
    else:
        plt.close()


def multibiplot_2axis(data, labels, data2, labels2, y_range, z_range, pdf=None,
                      zlabel='Energy [eV]', zlabel2=' ', zrange=(-default_range, default_range),
                      zrange2=(-default_range, default_range), title=None, show_plot=True, direction=0,
                      baseline=None):

    #plt.figure(3)
    f, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    if direction == 0:
        ax1.set_xlabel('X [Å]')
        for i, dat in enumerate(data):
            data[i] = np.array(dat).reshape([len(y_range), len(z_range)])[len(z_range)//2]
        for i, dat in enumerate(data2):
            data2[i] = np.array(dat).reshape([len(y_range), len(z_range)])[len(z_range)//2]

    if direction == 1:
        ax1.set_xlabel('Y [Å]')
        for i, dat in enumerate(data):
            data[i] = np.array(dat).reshape([len(y_range), len(z_range)])[:,len(y_range)//2]
            # data2 = np.array(data2).reshape([len(y_range), len(z_range)])[:,len(y_range)//2]
        for i, dat in enumerate(data2):
            data2[i] = np.array(dat).reshape([len(y_range), len(z_range)])[:,len(y_range)//2]

    plt.title(title)
    ax1.set_xlim([-4, 4])
    if zrange is not None:
        ax1.set_ylim([zrange[0], zrange[1]])
        ax1.set_ylabel(zlabel)
    if zrange2 is not None:
        ax2.set_ylim([zrange2[0], zrange2[1]])
        ax2.set_ylabel(zlabel2)

    lns1 = []
    for dat, l in zip(data, labels):
        lns1.append(ax1.plot(z_range, dat, label=l))
        # plt.plot(z_range, data2, label=label2)

    lns2 = []
    for dat, l in zip(data2, labels2):
        lns2.append(ax2.plot(z_range, dat, label=l, color='grey'))

    lns = [j for i in lns1 for j in i] + [j for i in lns2 for j in i]
    labs = [l.get_label() for l in lns]
    plt.legend(lns, labs, loc=0)

    if baseline is not None:
        ax1.axhline(baseline, color='k', linestyle=':')

    if pdf is not None:
        pdf.savefig()

    if show_plot:
        plt.show()
    else:
        plt.close()

def get_data(total_data, y_range, z_range, points,
             type=None, property=None, interpolate=False, average=False, states_orders=None):

    try:
        nval = len(total_data[0][type][property])
    except TypeError:
        nval = 1

    data_list = []
    for i in range(nval):
        if average:
            data_list.append([np.average(data[type][property]) for data in total_data])
        else:
            data_list.append(np.array([data[type][property][i] for data in total_data]))

    if states_orders is not None:
        data_list = correct_order_list(data_list, states_orders)

    if interpolate:
        for i in range(nval):
            data_list[i] = interpolate_data(points, data_list[i], y_range, z_range)

    return data_list


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

# define order states
states_orders = [get_order_states_list(data['states_info']) for data in total_data]
# states_orders = [get_order_states_list(data['states_info'], eps_energy=0) for data in total_data]


########################  D_i  ######################

data_1 = [np.average(data['diabatic_energies']['E_LE']) for data in total_data]
data_2 = [np.average(data['diabatic_energies']['E_CT']) for data in total_data]

if interpolate:
    data_1 = interpolate_data(points, data_1, y_range, z_range)
    data_2 = interpolate_data(points, data_2, y_range, z_range)

e_le = np.array(data_1)
e_ct = np.array(data_2)

with PdfPages(folder + 'LE_CT.pdf') as pdf:
    triplot2([e_le, e_ct], ['$E_{LE}$', '$E_{CT}$'], y_range, z_range, pdf=pdf, wireframe=True,
             show_plot=args.show_plots, zlevels=None)
    biplot(e_le, e_ct, '$E_{LE}$', '$E_{CT}$', y_range, z_range, show_plot=args.show_plots,
           pdf=pdf, direction=0, zrange=[8, 14])
    biplot(e_le, e_ct, '$E_{LE}$', '$E_{CT}$', y_range, z_range, show_plot=args.show_plots,
           pdf=pdf, direction=1, zrange=[8, 14])


de_le, = get_data(total_data, y_range, z_range, points,
                  type='diabatic_energies', property='dE_LE',
                  interpolate=interpolate, average=True)

de_ct, = get_data(total_data, y_range, z_range, points,
                  type='diabatic_energies', property='dE_CT',
                  interpolate=interpolate, average=True)

d2le = get_data(total_data, y_range, z_range, points,
                type='asymmetries', property='delta2_le',
                interpolate=interpolate, average=False,
                states_orders=states_orders)

d2ct = get_data(total_data, y_range, z_range, points,
                type='asymmetries', property='delta2_ct',
                interpolate=interpolate, average=False,
                states_orders=states_orders)

lmb2 = get_data(total_data, y_range, z_range, points,
                type='asymmetries', property='lambda2',
                interpolate=interpolate, average=False,
                states_orders=states_orders)


di_1 = lmb2[0]*(e_ct - e_le) + d2le[0]*de_le + d2ct[0]*de_ct
di_2 = lmb2[1]*(e_ct - e_le) + d2le[1]*de_le + d2ct[1]*de_ct


with PdfPages(folder + 'delta.pdf') as pdf:
    triplot2([di_1, di_2], ['$\Delta_1$', '$\Delta_2$'], y_range, z_range, pdf=pdf, wireframe=True,
             show_plot=args.show_plots, zlevels=np.arange(-0.4, 0.4 + 0.02, 0.02))
    biplot(di_1, di_2, '$\Delta_1$', '$\Delta_2$', y_range, z_range, show_plot=args.show_plots,
           pdf=pdf, direction=0, zrange=[0.0, 0.35])
    biplot(di_1, di_2, '$\Delta_1$', '$\Delta_2$', y_range, z_range, show_plot=args.show_plots,
           pdf=pdf, direction=1, zrange=[0.0, 0.35])


########################  omega_DC  ######################

w_dc = get_data(total_data, y_range, z_range, points,
                type='diabatic_contributions', property='W_DC',
                interpolate=interpolate, average=False,
                states_orders=states_orders)

o_dc_1 = np.sqrt((1 - lmb2[0])**2 - d2le[0]**2) * w_dc[0]
o_dc_2 = np.sqrt((1 - lmb2[1])**2 - d2le[1]**2) * w_dc[1]

with PdfPages(folder + 'omega_dc.pdf') as pdf:
    triplot2([o_dc_1, o_dc_2], ['$\Omega_{DC}^{(1)}$', '$\Omega_{DC}^{(2)}$'], y_range, z_range, pdf=pdf, wireframe=True,
             show_plot=args.show_plots, zlevels=np.arange(-0.4, 0.4 + 0.05, 0.05))
    biplot(o_dc_1, o_dc_2, '$\Omega_{DC}^{(1)}$', '$\Omega_{DC}^{(2)}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=0)
    biplot(o_dc_1, o_dc_2, '$\Omega_{DC}^{(1)}$', '$\Omega_{DC}^{(2)}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=1)


########################  omega_CT  ######################

w_ct = get_data(total_data, y_range, z_range, points,
                type='diabatic_contributions', property='W_CT',
                interpolate=interpolate, average=False,
                states_orders=states_orders)

o_ct_1 = np.sqrt(lmb2[0]**2 - d2ct[0]**2) * w_ct[0]
o_ct_2 = np.sqrt(lmb2[1]**2 - d2ct[1]**2) * w_ct[1]

with PdfPages(folder + 'omega_ct.pdf') as pdf:
    triplot2([o_ct_1, o_ct_2], ['$\Omega_{CT}^{(1)}$', '$\Omega_{CT}^{(2)}$'], y_range, z_range, pdf=pdf, wireframe=True,
             show_plot=args.show_plots, zlevels=np.arange(-0.02, 0.02 + 0.005, 0.005))
    biplot(o_ct_1, o_ct_2, '$\Omega_{CT}^{(1)}$', '$\Omega_{CT}^{(2)}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=0)
    biplot(o_ct_1, o_ct_2, '$\Omega_{CT}^{(1)}$', '$\Omega_{CT}^{(2)}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=1)


########################  omega_h  ######################

o_h = get_data(total_data, y_range, z_range, points,
               type='diabatic_contributions', property='Omega_h',
               interpolate=interpolate, average=False,
               states_orders=states_orders)

w_h = get_data(total_data, y_range, z_range, points,
               type='diabatic_contributions', property='W_h',
               interpolate=interpolate, average=False,
               states_orders=states_orders)

oh_1 = 2*np.sqrt(lmb2[0] * (1-lmb2[0]) + lmb2[0]*d2le[0] + (1-lmb2[0])*d2ct[0] + d2le[0]*d2ct[0]) * w_h[0]
oh_2 = 2*np.sqrt(lmb2[1] * (1-lmb2[1]) + lmb2[1]*d2le[1] + (1-lmb2[1])*d2ct[1] + d2le[1]*d2ct[1]) * w_h[1]


with PdfPages(folder + 'omega_h.pdf') as pdf:
    triplot2([o_h[0], o_h[1]], ['$\Omega_{h}^{(1)}$', '$\Omega_{h}^{(2)}$'], y_range, z_range, pdf=pdf, wireframe=True,
             show_plot=args.show_plots, zlevels=None)
    biplot(o_h[0], o_h[1], '$\Omega_{h}^{(1)}$', '$\Omega_{h}^{(2)}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=0)
    biplot(o_h[0], o_h[1], '$\Omega_{h}^{(1)}$', '$\Omega_{h}^{(2)}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=1)


with PdfPages(folder + 'omega_hb.pdf') as pdf:
    triplot2([oh_1, oh_2], ['$\Omega_{h}^{(1)}$', '$\Omega_{h}^{(2)}$'], y_range, z_range, pdf=pdf, wireframe=True,
             show_plot=args.show_plots, zlevels=None)
    biplot(oh_1, oh_2, '$\Omega_{h}^{(1)}$', '$\Omega_{h}^{(2)}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=0)
    biplot(oh_1, oh_2, '$\Omega_{h}^{(1)}$', '$\Omega_{h}^{(2)}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=1)


########################  omega_e  ######################

o_e = get_data(total_data, y_range, z_range, points,
               type='diabatic_contributions', property='Omega_e',
               interpolate=interpolate, average=False,
               states_orders=states_orders)

w_e = get_data(total_data, y_range, z_range, points,
               type='diabatic_contributions', property='W_e',
               interpolate=interpolate, average=False,
               states_orders=states_orders)

oe_1 = 2 * np.sqrt(lmb2[0] * (1-lmb2[0]) + lmb2[0]*d2le[0] - (1-lmb2[0])*d2ct[0] - d2le[0]*d2ct[0]) * w_e[0]
oe_2 = 2 * np.sqrt(lmb2[1] * (1-lmb2[1]) + lmb2[1]*d2le[1] - (1-lmb2[1])*d2ct[1] - d2le[1]*d2ct[1]) * w_e[1]


with PdfPages(folder + 'omega_e.pdf') as pdf:
    triplot2([o_e[0], o_e[1]], ['$\Omega_{e}^{(1)}$', '$\Omega_{e}^{(2)}$'], y_range, z_range, pdf=pdf, wireframe=True,
             show_plot=args.show_plots, zlevels=None)
    biplot(o_e[0], o_e[1], '$\Omega_{e}^{(1)}$', '$\Omega_{e}^{(2)}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=0)
    biplot(o_e[0], o_e[1], '$\Omega_{e}^{(1)}$', '$\Omega_{e}^{(2)}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=1)


with PdfPages(folder + 'omega_eb.pdf') as pdf:
    triplot2([oe_1, oe_2], ['$\Omega_{e}^{(1)}$', '$\Omega_{e}^{(2)}$'], y_range, z_range, pdf=pdf, wireframe=True,
             show_plot=args.show_plots, zlevels=None)
    biplot(oe_1, oe_2, '$\Omega_{e}^{(1)}$', '$\Omega_{e}^{(2)}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=0)
    biplot(oe_1, oe_2, '$\Omega_{e}^{(1)}$', '$\Omega_{e}^{(2)}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=1)


######################### adiabatic energies ############################

osx_1 = o_h[0] + o_e[0]
osx_2 = o_h[1] + o_e[1]

with PdfPages(folder + 'omega_sx.pdf') as pdf:
    triplot2([osx_1, osx_2], ['$\Omega_{SX}^{(1)}$', '$\Omega_{SX}^{(2)}$'], y_range, z_range, pdf=pdf, wireframe=True,
             show_plot=args.show_plots, zlevels=np.arange(-0.12, 0.12 + 0.005, 0.005))
    biplot(osx_1, osx_2, '$\Omega_{SX}^{(1)}$', '$\Omega_{e}^{(2)}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=0, zrange=[-1.5, 1.5])
    biplot(osx_1, osx_2, '$\Omega_{SX}^{(1)}$', '$\Omega_{e}^{(2)}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=1, zrange=[-1.5, 1.5])


with PdfPages(folder + 'omega_sx2.pdf') as pdf:
    triplot2([osx_1, osx_2], ['$\Omega_{SX}^{(1)}$', '$\Omega_{SX}^{(2)}$'], y_range, z_range, pdf=pdf, wireframe=True,
             show_plot=args.show_plots, zlevels=np.arange(-0.01, 0.01 + 0.001, 0.001))
    biplot(osx_1, osx_2, '$\Omega_{SX}^{(1)}$', '$\Omega_{e}^{(2)}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=0, zrange=[-0.2, 0.2])
    biplot(osx_1, osx_2, '$\Omega_{SX}^{(1)}$', '$\Omega_{e}^{(2)}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=1, zrange=[-0.2, 0.2])


######################### figure 6 ############################


with PdfPages(folder + 'figure6.pdf') as pdf:
    multibiplot_2axis([o_dc_1, osx_1, o_ct_1, di_1], ['$\Omega^{(1)}_{DC}$', '$\Omega_{SX}^{(1)}$', '$\Omega_{CT}^{(1)}$', '$\Delta_1$'],
                      [lmb2[0]], ['$\lambda^2_1$'], y_range, z_range, pdf=pdf, zlabel2='$\lambda^2$',
                      zlabel='Energy [eV]', zrange=(-1.5, 1.5), zrange2=(0.0, 0.2), title='State 1',
                      show_plot=args.show_plots, direction=0, baseline=0.0)

    multibiplot_2axis([o_dc_1, osx_1, o_ct_1, di_1], ['$\Omega^{(1)}_{DC}$', '$\Omega_{SX}^{(1)}$', '$\Omega_{CT}^{(1)}$', '$\Delta_1$'],
                      [lmb2[0]], ['$\lambda^2_1$'], y_range, z_range, pdf=pdf, zlabel2='$\lambda^2$',
                      zlabel='Energy [eV]', zrange=(-1.5, 1.5), zrange2=(0.0, 0.2), title='State 1',
                      show_plot=args.show_plots, direction=1, baseline=0.0)

    multibiplot_2axis([o_dc_2, osx_2, o_ct_2, di_2], ['$\Omega^{(2)}_{DC}$', '$\Omega_{SX}^{(2)}$', '$\Omega_{CT}^{(2)}$', '$\Delta_2$'],
                      [lmb2[1]], ['$\lambda^2_2$'], y_range, z_range, pdf=pdf, zlabel2='$\lambda^2$',
                      zlabel='Energy [eV]', zrange=(-1.5, 1.5), zrange2=(0.0, 0.2), title='State 2',
                      show_plot=args.show_plots, direction=0, baseline=0.0)

    multibiplot_2axis([o_dc_2, osx_2, o_ct_2, di_2], ['$\Omega^{(2)}_{DC}$', '$\Omega_{SX}^{(2)}$', '$\Omega_{CT}^{(2)}$', '$\Delta_2$'],
                      [lmb2[1]], ['$\lambda^2_2$'], y_range, z_range, pdf=pdf, zlabel2='$\lambda^2$',
                      zlabel='Energy [eV]', zrange=(-1.5, 1.5), zrange2=(0.0, 0.2), title='State 2',
                      show_plot=args.show_plots, direction=1, baseline=0.0)


######################### adiabatic energies ############################

data_1 = [data['adiabatic_energies']['E_1'] for data in total_data]
data_2 = [data['adiabatic_energies']['E_2'] for data in total_data]
data_3 = [data['adiabatic_energies']['E_3'] for data in total_data]
data_4 = [data['adiabatic_energies']['E_4'] for data in total_data]

data_1, data_2, data_3, data_4 = correct_order_list([data_1, data_2, data_3, data_4], states_orders)

if interpolate:
    data_1 = interpolate_data(points, data_1, y_range, z_range)
    data_2 = interpolate_data(points, data_2, y_range, z_range)
    data_3 = interpolate_data(points, data_3, y_range, z_range)
    data_4 = interpolate_data(points, data_4, y_range, z_range)

e_1 = np.array(data_1)
e_2 = np.array(data_2)
e_3 = np.array(data_3)
e_4 = np.array(data_4)


with PdfPages(folder + 'figure4.pdf') as pdf:
    triplot3([e_1, e_2], ['$^1A_g$', '$^1B_{3u}$'], y_range, z_range, pdf=pdf, wireframe=True,
             show_plot=args.show_plots, zlevels=None, colors=['black', 'red'])
    biplot(e_1, e_2, '$^1A_g$', '$^1B_{3u}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=0, zrange=[7.6, 9.2], colors=['black', 'red'])
    biplot(e_1, e_2, '$^1A_g$', '$^1B_{3u}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=1, zrange=[7.6, 9.2], colors=['black', 'red'])

e_1b = e_le + di_1 + o_dc_1 + o_ct_1 + o_h[0] + o_e[0]
e_2b = e_le + di_2 + o_dc_2 + o_ct_2 + o_h[1] + o_e[1]

with PdfPages(folder + 'test.pdf') as pdf:
    #triplot2([e_1-e_1b, e_2-e_2b], ['$E^{(1)}$', '$E^{(2)}$'], y_range, z_range, pdf=pdf, wireframe=True,
    #         show_plot=args.show_plots, zlevels=None)
    biplot(e_1-e_1b, e_2-e_2b, '$E^{(1)}$', '$E^{(2)}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=0)
    biplot(e_1-e_1b, e_2-e_2b, '$E^{(1)}$', '$E^{(2)}$', y_range, z_range,
           show_plot=args.show_plots, pdf=pdf, direction=1)


exit()

########################  omega 1,J,K  ######################

o_one_dc = get_data(total_data, y_range, z_range, points,
                    type='diabatic_contributions', property='Omega_One_dc',
                    interpolate=interpolate, average=False,
                    states_orders=states_orders)

o_j_dc = get_data(total_data, y_range, z_range, points,
                  type='diabatic_contributions', property='Omega_J_dc',
                  interpolate=interpolate, average=False,
                  states_orders=states_orders)
o_k_dc = get_data(total_data, y_range, z_range, points,
                  type='diabatic_contributions', property='Omega_K_dc',
                  interpolate=interpolate, average=False,
                  states_orders=states_orders)



with PdfPages(folder + 'omega_other.pdf') as pdf:
    #triplot2([o_dc_one_1, o_dc_one_2], ['$\Omega_{e}^{(1)}$', '$\Omega_{e}^{(2)}$'], y_range, z_range, pdf=pdf, wireframe=True,
    #          show_plot=args.show_plots, zlevels=None)
    multibiplot([o_one_dc[0], o_j_dc[0], o_k_dc[0]],
                ['$\Omega_{DC,1}^{(1)}$', '$\Omega_{DC,J}^{(1)}$', '$\Omega_{DC,K}^{(1)}$'], y_range, z_range,
                 show_plot=args.show_plots, pdf=pdf, direction=0)
    multibiplot([o_one_dc[0], o_j_dc[0], o_k_dc[0]],
                ['$\Omega_{DC,1}^{(1)}$', '$\Omega_{DC,J}^{(1)}$', '$\Omega_{DC,K}^{(1)}$'], y_range, z_range,
                 show_plot=args.show_plots, pdf=pdf, direction=1)


exit()

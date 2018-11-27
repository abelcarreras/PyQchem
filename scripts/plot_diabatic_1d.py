import pickle
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyqchem.order_states import get_order_states_list, correct_order_list



# Argument parser
parser = argparse.ArgumentParser(description='Plot 3D data')
parser.add_argument('filename', metavar='filename', type=str,
                    help='filename for input')

parser.add_argument('--output_folder', metavar='distance', type=str, default='1d_plot/',
                    help='folder to store PDF plots')

parser.add_argument('--show_plots', action='store_true',
                   help='show plots while running')


args = parser.parse_args()


def multiplot(data , labels, x_range, title=None, factor=None, range_y=(-1, 1),
              ylabel='Energy [eV]', show_plots=False, colors=None, baseline=None):

    if factor is not None:
        for i, dat in enumerate(data):
            data[i] = factor * np.array(dat)

    f = plt.figure(1)

    plt.title(title)
    for i, dat in enumerate(data):
        if colors is not None:
            plt.plot(x_range, dat, label=labels[i], color=colors[i])
        else:
            plt.plot(x_range, dat, label=labels[i])

    plt.xlim(3.5, 6.0)
    if range_y is not None:
        plt.ylim(*range_y)
    plt.legend()
    plt.xlabel('Z [Å]')
    plt.ylabel(ylabel)

    if baseline is not None:
        plt.axhline(baseline, color='k', linestyle=':')

    if show_plots:
        plt.show()
    else:
        plt.close()

    return f


def multiplot_2axis(data , labels, data2, labels2, x_range, title=None, factor=None, range_y=(-1, 1),
                    ylabel='Energy [eV]', show_plots=False, colors=None, style1=None, style2=None,
                    colors2=None, ylabel2=None, range_y2=(-1, 1), baseline=None):

    if style1 is None:
        style1 = [None for i in range(len(data))]

    if style2 is None:
        style2 = [None for i in range(len(data2))]


    if factor is not None:
        for i, dat in enumerate(data):
            data[i] = factor[0] * np.array(dat)
        for i, dat in enumerate(data2):
            data2[i] = factor[1] * np.array(dat)

    f, ax1 = plt.subplots()

    plt.title(title)
    plt.xlim(3.5, 6.0)

    lns1 = []
    for i, dat in enumerate(data):
        if colors is not None:
            lns1.append(ax1.plot(x_range, dat, label=labels[i], color=colors[i], linestyle=style1[i]))
        else:
            lns1.append(ax1.plot(x_range, dat, label=labels[i], linestyle=style1[i]))

    if range_y is not None:
        ax1.set_ylim(*range_y)
    ax1.set_ylabel(ylabel)

    ax2 = ax1.twinx()

    lns2 = []
    for i, dat in enumerate(data2):
        if colors2 is not None:
            lns2.append(ax2.plot(x_range, dat, label=labels2[i], color=colors2[i]))
        else:
            lns2.append(ax2.plot(x_range, dat, label=labels2[i]))

    if range_y2 is not None:
        ax2.set_ylim(*range_y2)

    ax2.set_ylabel(ylabel2)

    lns = [j for i in lns1 for j in i] + [j for i in lns2 for j in i]
    labs = [l.get_label() for l in lns]
    plt.legend(lns, labs, loc=0)

    if baseline is not None:
        plt.axhline(baseline, color='k', linestyle=':')

    ax1.set_xlabel('Z [Å]')

    if show_plots:
        plt.show()
    else:
        plt.close()

    return f


#############################
folder = args.output_folder
#############################


with open(args.filename, 'rb') as input:
    calculation_data = pickle.load(input)
    print('Loaded data from calculation_data.pkl')

total_data = []
d_coordinate = []
for distance in calculation_data['distance']:
        print('{}'.format(distance))
        if '{}'.format(distance) in calculation_data:
            data = calculation_data['{}'.format(distance)]
            total_data.append(data)
            d_coordinate.append(distance)

states_orders = [get_order_states_list(data['states_info']) for data in total_data]

########################### W_DC ###########################

data_1 = [data['diabatic_contributions']['W_DC'][0] for data in total_data]
data_2 = [data['diabatic_contributions']['W_DC'][1] for data in total_data]
data_3 = [data['diabatic_contributions']['W_DC'][2] for data in total_data]
data_4 = [data['diabatic_contributions']['W_DC'][3] for data in total_data]

data_1, data_2, data_3, data_4 = correct_order_list([data_1, data_2, data_3, data_4], states_orders)

f = multiplot([data_1, data_2], ['W_DC_1', 'W_DC_2'], d_coordinate, title='W_DC', show_plots=args.show_plots)
f.savefig(folder + "W_DC.pdf", bbox_inches='tight')

wdc1 = np.array(data_1)
wdc2 = np.array(data_2)
wdc3 = np.array(data_3)
wdc4 = np.array(data_4)


########################### W_CT ###########################

data_1 = [data['diabatic_contributions']['W_CT'][0] for data in total_data]
data_2 = [data['diabatic_contributions']['W_CT'][1] for data in total_data]
data_3 = [data['diabatic_contributions']['W_CT'][2] for data in total_data]
data_4 = [data['diabatic_contributions']['W_CT'][3] for data in total_data]

data_1, data_2, data_3, data_4 = correct_order_list([data_1, data_2, data_3, data_4], states_orders)

f = multiplot([data_1, data_2], ['W_CT_1', 'W_CT_2'], d_coordinate, title='W_CT', show_plots=args.show_plots)
f.savefig(folder + "W_CT.pdf", bbox_inches='tight')

wct1 = np.array(data_1)
wct2 = np.array(data_2)
wct3 = np.array(data_3)
wct4 = np.array(data_4)

########################### W_e ############################

data_1 = [data['diabatic_contributions']['W_e'][0] for data in total_data]
data_2 = [data['diabatic_contributions']['W_e'][1] for data in total_data]
data_3 = [data['diabatic_contributions']['W_e'][2] for data in total_data]
data_4 = [data['diabatic_contributions']['W_e'][3] for data in total_data]

data_1, data_2, data_3, data_4 = correct_order_list([data_1, data_2, data_3, data_4], states_orders)

f = multiplot([data_1, data_2], ['W_e_1', 'W_e_2'], d_coordinate, title='W_e', show_plots=args.show_plots)
f.savefig(folder + "W_e.pdf", bbox_inches='tight')

we1 = np.array(data_1)
we2 = np.array(data_2)
we3 = np.array(data_3)
we4 = np.array(data_4)

########################### W_h ############################

data_1 = [data['diabatic_contributions']['W_h'][0] for data in total_data]
data_2 = [data['diabatic_contributions']['W_h'][1] for data in total_data]
data_3 = [data['diabatic_contributions']['W_h'][2] for data in total_data]
data_4 = [data['diabatic_contributions']['W_h'][3] for data in total_data]

data_1, data_2, data_3, data_4 = correct_order_list([data_1, data_2, data_3, data_4], states_orders)

f = multiplot([data_1, data_2], ['W_h_1', 'W_h_2'], d_coordinate, title='W_h')
f.savefig(folder + "W_h.pdf", bbox_inches='tight')

wh1 = np.array(data_1)
wh2 = np.array(data_2)
wh3 = np.array(data_3)
wh4 = np.array(data_4)

###################### W per state ###########################

f = multiplot([wh1, we1, wct1, wdc1], ['W_h', 'W_e', 'W_CT', 'W_DC'], d_coordinate,
              title='State 1', show_plots=args.show_plots)
f.savefig(folder + "W_1.pdf", bbox_inches='tight')

f = multiplot([wh2, we2, wct2, wdc2], ['W_h', 'W_e', 'W_CT', 'W_DC'], d_coordinate,
              title='State 2', show_plots=args.show_plots)
f.savefig(folder + "W_2.pdf", bbox_inches='tight')


###################################### lambda ###########################

data_1 = [data['lambda'][0] for data in total_data]
data_2 = [data['lambda'][1] for data in total_data]
data_3 = [data['lambda'][2] for data in total_data]
data_4 = [data['lambda'][3] for data in total_data]

#data_1 = []
#for i, e in enumerate(np.array([data['lambda'] for data in total_data]).T[0]):
#    data_1.append(e)

#data_2 = []
#for i, e in enumerate(np.array([data['lambda'] for data in total_data]).T[1]):
#    data_2.append(e)

data_1, data_2, data_3, data_4 = correct_order_list([data_1, data_2, data_3, data_4], states_orders)

f = multiplot([data_1, data_2], ['lambda 1', 'lambda 2'], d_coordinate, range_y=[0, 0.6],
              ylabel='', title='lambda', show_plots=args.show_plots)
f.savefig(folder + "lambda.pdf", bbox_inches='tight')

l1 = np.array(data_1)
l2 = np.array(data_2)
l3 = np.array(data_3)
l4 = np.array(data_4)

f = multiplot([np.square(data_1), np.square(data_2)], ['lambda^2 1', 'lambda^2 2'], d_coordinate,
              range_y=[0, 0.6], ylabel='', title='lambda^2', show_plots=args.show_plots)
f.savefig(folder + "lambda2.pdf", bbox_inches='tight')


##########################  OMEGA  #######################


data_1 = 2 * l1 * np.sqrt(1 - l1**2) * we1
data_2 = 2 * l2 * np.sqrt(1 - l2**2) * we2

data_3 = 2 * l1 * np.sqrt(1 - l1**2) * wh1
data_4 = 2 * l2 * np.sqrt(1 - l2**2) * wh2

f = multiplot([data_1, data_2, data_3, data_4], ['Omega_e 1', 'Omega_e 2', 'Omega_h 1', 'Omega_h 2'],
              d_coordinate, ylabel='', title='Omega', show_plots=args.show_plots)
f.savefig(folder + "Omega.pdf", bbox_inches='tight')

oe_1 = np.array(data_1)
oe_2 = np.array(data_2)
oh_1 = np.array(data_3)
oh_2 = np.array(data_4)


########################## E1 (Superexchange) ###########################

data_1 = 2 * l1 * np.sqrt(1 - l1**2) * (we1 + wh1)
data_2 = 2 * l2 * np.sqrt(1 - l2**2) * (we2 + wh2)

f = multiplot([data_1, data_2], ['E1-1', 'E1-2'], d_coordinate, show_plots=args.show_plots)
f.savefig(folder + "E1(superexchange).pdf", bbox_inches='tight')

e_11 = np.array(data_1)
e_12 = np.array(data_2)


#######################  diabatic_energies  ######################

data_1 = [np.average(data['diabatic_energies']['E_LE']) for data in total_data]
data_2 = [np.average(data['diabatic_energies']['E_CT']) for data in total_data]

e_le = np.array(data_1)
e_ct = np.array(data_2)

f = multiplot([data_1, data_2], ['$E_{LE}$', '$E_{CT}$'], d_coordinate, title='Diabatic energies',
              range_y=[6, 13], show_plots=args.show_plots)
f.savefig(folder + "figure1b.pdf", bbox_inches='tight')

#######################  Second order term  (E2)  #####################

data_1 = l1**2 * (e_ct - e_le + wct1 - wdc1)
data_2 = l2**2 * (e_ct - e_le + wct2 - wdc2)

f = multiplot([data_1, data_2], ['E2-1', 'E2-2'], d_coordinate)
f.savefig(folder + "E2.pdf", bbox_inches='tight')

e_21 = np.array(data_1)
e_22 = np.array(data_2)


###################### Diabatic energies ##################

data_1 = [data['diabatic_energies']['V_DC'] for data in total_data]
data_2 = [data['diabatic_energies']['V_CT'] for data in total_data]
data_3 = [data['diabatic_energies']['V_e'][0] for data in total_data]
data_4 = [data['diabatic_energies']['V_h'][0] for data in total_data]
data_5 = [data['diabatic_energies']['V_e'][1] for data in total_data]
data_6 = [data['diabatic_energies']['V_h'][1] for data in total_data]

f = multiplot([data_1, data_2,  data_3,  data_4,  data_5, data_6],
              ['V_DC', "V_CT", "V_e", "V_h", "V_e'", "V_h'"],
              d_coordinate, title='Diabatic energies', show_plots=args.show_plots)

f.savefig(folder + "diabatic_energies.pdf", bbox_inches='tight')

###################### Adiabatic energies ##################

data_1 = [data['adiabatic_energies']['E_1'] for data in total_data]
data_2 = [data['adiabatic_energies']['E_2'] for data in total_data]
data_3 = [data['adiabatic_energies']['E_3'] for data in total_data]
data_4 = [data['adiabatic_energies']['E_4'] for data in total_data]

data_1, data_2, data_3, data_4 = correct_order_list([data_1, data_2, data_3, data_4], states_orders)

f = multiplot([data_1, data_2,  data_3,  data_4,],
              #['$1 - {}^1A_g$', '$1 - {}^1B_{3u}$', '$2- {}^1B_{3u}$', '$2 - {}^1A_g$'],
              ['${}^1A_g$', '${}^1B_{3u}$',None , None],
              d_coordinate, title='Adiabatic energies', range_y=[6, 13], show_plots=args.show_plots,
              colors=['black', 'red', 'black', 'red']
              )

f.savefig(folder + "figure1a.pdf", bbox_inches='tight')

e1 = np.array(data_1)
e2 = np.array(data_2)
e3 = np.array(data_3)
e4 = np.array(data_4)


#######################  J Coloumb ######################


mu = 1.7528
angs2au = 1.889725989
har2ev = 27.2113962
j_coul = np.array([(mu**2)/(z * angs2au)**3 * har2ev for z in d_coordinate])

#######################  adiabatic_energies (calculated) ######################

e1c = e_le + wdc1 + e_11 + e_21
e2c = e_le + wdc2 + e_12 + e_22

f = multiplot([e1-e1c, e2-e2c],
              ['diff 1', 'diff 2'],
              d_coordinate, title='Adiabatic energies difference\n (original-calculated)',
              range_y=None, show_plots=args.show_plots)

f.savefig(folder + "test.pdf", bbox_inches='tight')

#######################  Figure 2  #####################

f = multiplot_2axis([wdc1, e_11, e_21, -j_coul], ['$W^{(1)}_{DC}$', '$\Omega^{(1)}$', '$E(\lambda^2_1)$', '$J_{Coul}$'],
                    [np.square(l1)], ['$\lambda_1^2$'], d_coordinate, ylabel2='$\lambda^2$',
                    colors=[None, None, None, '#1f77b4'], style1=['-', '-', '-', '--'],
                    range_y2=[-0.3, 0.3], title='State 1', colors2=['grey'], baseline=0)
f.savefig(folder + "figure2a.pdf", bbox_inches='tight')

f = multiplot_2axis([wdc2, e_12, e_22, j_coul], ['$W^{(2)}_{DC}$', '$\Omega^{(2)}$', '$E(\lambda^2_2)$', '$J_{Coul}$'],
                    [np.square(l2)], ['$\lambda_2^2$'], d_coordinate, ylabel2='$\lambda^2$',
                    colors=[None, None, None, '#1f77b4'], style1=['-', '-', '-', '--'],
                    range_y2=[-0.3, 0.3], title='State 2', colors2=['grey'], baseline=0)

f.savefig(folder + "figure2b.pdf", bbox_inches='tight')

f = multiplot_2axis([oe_1, oh_1, e_11], ['$\Omega_e^{(1)}$', '$\Omega_h^{(1)}$', '$\Omega^{(1)}$'],
                    [np.square(l1)], ['$\lambda_1^2$'], d_coordinate, ylabel2='$\lambda^2$',
                    range_y2=[-0.3, 0.3], title='State 1', colors2=['grey'], baseline=0)

f.savefig(folder + "figure2c.pdf", bbox_inches='tight')

f = multiplot_2axis([oe_2, oh_2, e_12], ['$\Omega_e^{(2)}$', '$\Omega_h^{(2)}$', '$\Omega^{(2)}$'],
                    [np.square(l2)], ['$\lambda_2^2$'], d_coordinate, ylabel2='$\lambda^2$',
                    range_y2=[-0.3, 0.3], title='State 2', colors2=['grey'], baseline=0)

f.savefig(folder + "figure2d.pdf", bbox_inches='tight')


f = multiplot_2axis([oe_1, oh_1], ['$\Omega_e^{(1)}$', '$\Omega_h^{(1)}$'],
                    [np.square(l1)], ['$\lambda_1^2$'], d_coordinate, ylabel2='$\lambda^2$',
                    style1=['-', '-.'], colors=['black', 'black'],
                    range_y2=[-0.3, 0.3], title='State 1', colors2=['grey'], baseline=0)

f.savefig(folder + "figure2c_x.pdf", bbox_inches='tight')

f = multiplot_2axis([oe_2, oh_2], ['$\Omega_e^{(2)}$', '$\Omega_h^{(2)}$'],
                    [np.square(l2)], ['$\lambda_2^2$'], d_coordinate, ylabel2='$\lambda^2$',
                    style1=['-', '-.'], colors=['red', 'red'],
                    range_y2=[-0.3, 0.3], title='State 2', colors2=['grey'], baseline=0)

f.savefig(folder + "figure2d_x.pdf", bbox_inches='tight')

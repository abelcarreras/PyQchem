import pickle
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pyqchem.order_states import get_order_states_list, correct_order_list
import matplotlib

matplotlib.rcParams.update({'font.size': 13})

# Argument parser
parser = argparse.ArgumentParser(description='Plot 3D data')
parser.add_argument('filename', metavar='filename', type=str,
                    help='filename for input')

parser.add_argument('--output_folder', metavar='distance', type=str, default='1d_plot/',
                    help='folder to store PDF plots')

parser.add_argument('--show_plots', action='store_true',
                   help='show plots while running')


args = parser.parse_args()


def multiplot(data , labels, x_range, title=None, factor=None, range_y=(-1.0, 1.0),
              ylabel='Energy [eV]', show_plots=False, colors=None, baseline=None, style=None):

    if style is None:
        style = [None for i in range(len(data))]

    if factor is not None:
        for i, dat in enumerate(data):
            data[i] = factor * np.array(dat)

    f = plt.figure(1)

    plt.title(title)
    for i, dat in enumerate(data):
        if colors is not None:
            plt.plot(x_range, dat, label=labels[i], color=colors[i], linestyle=style[i])
        else:
            plt.plot(x_range, dat, label=labels[i], linestyle=style[i])

    plt.xlim(3.5, 6.0)
    if range_y is not None:
        plt.ylim(*range_y)
    plt.legend()
    #leg = plt.legend(loc='center right')
    #leg.get_frame().set_linewidth(0.0)

    plt.xlabel('Z [Å]')
    plt.ylabel(ylabel)

    if baseline is not None:
        plt.axhline(baseline, color='k', linestyle=':')

    if show_plots:
        plt.show()
    else:
        plt.close()

    return f


def multiplot_2axis(data , labels, data2, labels2, x_range, title=None, factor=None, range_y=(-1.0, 1.0),
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
        ax1.axhline(baseline, color='k', linestyle=':')

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

states_orders = [get_order_states_list(data['states_info'], eps_energy=0) for data in total_data]

########################### W_DC ###########################

data_1 = [data['diabatic_contributions']['W_DC'][0] for data in total_data]
data_2 = [data['diabatic_contributions']['W_DC'][1] for data in total_data]
data_3 = [data['diabatic_contributions']['W_DC'][2] for data in total_data]
data_4 = [data['diabatic_contributions']['W_DC'][3] for data in total_data]

data_1, data_2, data_3, data_4 = correct_order_list([data_1, data_2, data_3, data_4], states_orders)

f = multiplot([data_1, data_2], ['$W_{DC}^{(1)}$', '$W_{DC}^{(2)}$'], d_coordinate, title='W_DC',
              show_plots=args.show_plots, range_y=[None, None])
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

f = multiplot([data_1, data_2], ['$W_{CT}^{(1)}$', '$W_{CT}^{(2)}$'], d_coordinate, title='W_CT',
              show_plots=args.show_plots, range_y=[None, None])
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

f = multiplot([data_1, data_2], ['$W_e^{(1)}$', '$W_e^{(2)}$'], d_coordinate, title='W_e',
              show_plots=args.show_plots, range_y=[None, None])
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

f = multiplot([data_1, data_2], ['$W_h^{(1)}$', '$W_h^{(2)}$'], d_coordinate, title='W_h',
              show_plots=args.show_plots, range_y=[None, None])
f.savefig(folder + "W_h.pdf", bbox_inches='tight')

wh1 = np.array(data_1)
wh2 = np.array(data_2)
wh3 = np.array(data_3)
wh4 = np.array(data_4)

###################### W per state ###########################

f = multiplot([wh1, we1, wct1, wdc1], ['$W_h$', '$W_e$', '$W_{CT}$', '$W_{DC}$'], d_coordinate,
              title='State 1', show_plots=args.show_plots, range_y=[None, None])
f.savefig(folder + "W_1.pdf", bbox_inches='tight')

f = multiplot([wh2, we2, wct2, wdc2], ['$W_h$', '$W_e$', '$W_{CT}$', '$W_{DC}$'], d_coordinate,
              title='State 2', show_plots=args.show_plots, range_y=[None, None])
f.savefig(folder + "W_2.pdf", bbox_inches='tight')


###################################### lambda ###########################

data_1 = [data['asymmetries']['lambda2'][0] for data in total_data]
data_2 = [data['asymmetries']['lambda2'][1] for data in total_data]
data_3 = [data['asymmetries']['lambda2'][2] for data in total_data]
data_4 = [data['asymmetries']['lambda2'][3] for data in total_data]

#data_1 = []
#for i, e in enumerate(np.array([data['lambda'] for data in total_data]).T[0]):
#    data_1.append(e)

#data_2 = []
#for i, e in enumerate(np.array([data['lambda'] for data in total_data]).T[1]):
#    data_2.append(e)

data_1, data_2, data_3, data_4 = correct_order_list([data_1, data_2, data_3, data_4], states_orders)

f = multiplot([data_1, data_2], ['$\lambda^{2(1)}$', '$\lambda^{2(2)}$'], d_coordinate, range_y=[None, None],
              ylabel='', title='lambda^2', show_plots=args.show_plots)
f.savefig(folder + "lambda2.pdf", bbox_inches='tight')

l1 = np.array(data_1)
l2 = np.array(data_2)
l3 = np.array(data_3)
l4 = np.array(data_4)

##########################  OMEGA  #######################

data_1 = [data['diabatic_contributions']['Omega_e'][0] for data in total_data]
data_2 = [data['diabatic_contributions']['Omega_e'][1] for data in total_data]

data_3 = [data['diabatic_contributions']['Omega_h'][0] for data in total_data]
data_4 = [data['diabatic_contributions']['Omega_h'][1] for data in total_data]

#data_1 = 2 * l1 * np.sqrt(1 - l1**2) * we1
#data_2 = 2 * l2 * np.sqrt(1 - l2**2) * we2

#data_3 = 2 * l1 * np.sqrt(1 - l1**2) * wh1
#data_4 = 2 * l2 * np.sqrt(1 - l2**2) * wh2

f = multiplot([data_1, data_2, data_3, data_4], ['$\Omega_e^{(1)}$', '$\Omega_e^{(2)}$', '$\Omega_h^{(1)}$', '$\Omega_h^{(2)}$'],
              d_coordinate, ylabel='', title='Omega', show_plots=args.show_plots, range_y=[None, None])
f.savefig(folder + "Omega.pdf", bbox_inches='tight')

oe_1 = np.array(data_1)
oe_2 = np.array(data_2)
oh_1 = np.array(data_3)
oh_2 = np.array(data_4)


########################## delta ###########################

e_le = np.array([data['diabatic_energies']['E_LE'] for data in total_data])
e_ct = np.array([data['diabatic_energies']['E_CT'] for data in total_data])
de_le = np.array([data['diabatic_energies']['dE_LE'] for data in total_data])
de_ct = np.array([data['diabatic_energies']['dE_CT'] for data in total_data])

d2le1 = np.array([data['asymmetries']['delta2_le'][0] for data in total_data])
d2le2 = np.array([data['asymmetries']['delta2_le'][1] for data in total_data])

d2ct1 = np.array([data['asymmetries']['delta2_ct'][0] for data in total_data])
d2ct2 = np.array([data['asymmetries']['delta2_ct'][1] for data in total_data])


data_1 = d2le1 * de_le + d2ct1 * de_ct + l1 * (e_ct - e_le)
data_2 = d2le2 * de_le + d2ct2 * de_ct + l2 * (e_ct - e_le)


f = multiplot([data_1, data_2], ['$\Delta_1$', '$\Delta_2$'], d_coordinate, show_plots=args.show_plots, #range_y=[None, None],
              title='$\Delta$')
f.savefig(folder + "delta.pdf", bbox_inches='tight')

f = multiplot([e_le, e_ct, de_le, de_ct], ['$E_{LE}$', '$E_{CT}$', '$\Delta E_{LE}$', '$\Delta E_{CT}$'],
              d_coordinate, show_plots=args.show_plots, range_y=[None, None], title='term1 pieces')
f.savefig(folder + "diabatic_1.pdf", bbox_inches='tight')

f = multiplot([d2le1, d2ct1, d2le2, d2ct2], ['$\delta^{2(1)}_{LE}$', '$\delta^{2(1)}_{CT}$',
                                             '$\delta^{2(2)}_{LE}$', '$\delta^{2(2)}_{CT}$'],
              d_coordinate, show_plots=args.show_plots, title='term1 pieces')
f.savefig(folder + "diabatic_2.pdf", bbox_inches='tight')


delta_1 = np.array(data_1)
delta_2 = np.array(data_2)

#######################  J Coloumb ######################

mu = 1.7528
angs2au = 1.889725989
har2ev = 27.2113962
j_coul = np.array([(mu**2)/(z * angs2au)**3 * har2ev for z in d_coordinate])

#######################  Omega DC  #####################

data_1 = np.sqrt((1-l1)**2 - d2le1**2)*wdc1
data_2 = np.sqrt((1-l2)**2 - d2le2**2)*wdc2

add_coul = True
if add_coul:
    f = multiplot([data_1, data_2, j_coul, -j_coul],
                  ['$\Omega_{DC}^{(1)}$', '$\Omega_{DC}^{(2)}$', '$J_{Coul}$', '$J_{Coul}$'],
                  d_coordinate, title='$\Omega_{DC}$', colors=['#1f77b4', '#ff7f0e', 'grey', 'grey'],
                  style=['-', '-', '--', '--'])
else:
    f = multiplot([data_1, data_2],
                  ['$\Omega_{DC}^{(1)}$', '$\Omega_{DC}^{(2)}$'],
                  d_coordinate, title='$\Omega_{DC}$')

f.savefig(folder + "omega_DC.pdf", bbox_inches='tight', show_plots=args.show_plots, range_y=[None, None])

o_dc1 = data_1
o_dc2 = data_2

#######################  Omega CT  #####################

data_1 = np.sqrt(l1**2 - d2ct1**2)*wct1
data_2 = np.sqrt(l1**2 - d2ct2**2)*wct2

f = multiplot([data_1, data_2], ['$\Omega_{CT}^{(1)}$', '$\Omega_{CT}^{(2)}$'], d_coordinate, title='$\Omega_{CT}$')
f.savefig(folder + "omega_CT.pdf", bbox_inches='tight', show_plots=args.show_plots, range_y=[None, None])

o_ct1 = data_1
o_ct2 = data_2

###################### Omega h ###########################

o_h1 = np.array([data['diabatic_contributions']['Omega_h'][0] for data in total_data])
o_h2 = np.array([data['diabatic_contributions']['Omega_h'][1] for data in total_data])

f = multiplot([o_h1, o_h2], ['$\Omega_{h}^{(1)}$', '$\Omega_{h}^{(2)}$'], d_coordinate, title='$\Omega_{h}$')
f.savefig(folder + "omega_h.pdf", bbox_inches='tight', show_plots=args.show_plots, range_y=[None, None])

###################### Omega e ###########################

o_e1 = np.array([data['diabatic_contributions']['Omega_e'][0] for data in total_data])
o_e2 = np.array([data['diabatic_contributions']['Omega_e'][1] for data in total_data])

f = multiplot([o_e1, o_e2], ['$\Omega_{e}^{(1)}$', '$\Omega_{e}^{(2)}$'], d_coordinate, title='$\Omega_{e}$')
f.savefig(folder + "omega_e.pdf", bbox_inches='tight', show_plots=args.show_plots, range_y=[None, None])

#######################  diabatic_energies  ######################

data_1 = [np.average(data['diabatic_energies']['E_LE']) for data in total_data]
data_2 = [np.average(data['diabatic_energies']['E_CT']) for data in total_data]

e_le = np.array(data_1)
e_ct = np.array(data_2)

f = multiplot([data_1, data_2], ['$E_{LE}$', '$E_{CT}$'], d_coordinate, #title='Diabatic energies',
              range_y=[6, 13], show_plots=args.show_plots)
f.savefig(folder + "figure2b.pdf", bbox_inches='tight')


###################### Diabatic energies ##################

data_1 = [data['diabatic_energies']['V_DC'] for data in total_data]
data_2 = [data['diabatic_energies']['V_CT'] for data in total_data]
data_3 = [data['diabatic_energies']['V_e'][0] for data in total_data]
data_4 = [data['diabatic_energies']['V_h'][0] for data in total_data]
data_5 = [data['diabatic_energies']['V_e'][1] for data in total_data]
data_6 = [data['diabatic_energies']['V_h'][1] for data in total_data]

f = multiplot([data_1, data_2,  data_3,  data_4,  data_5, data_6],
              ['$V_{DC}$', "$V_{CT}$", "$V_e$", "$V_h$", "$V_e'$", "$V_h'$"],
              d_coordinate, #title='Diabatic energies',
              show_plots=args.show_plots)

f.savefig(folder + "diabatic_energies.pdf", bbox_inches='tight')

###################### Adiabatic energies ##################

data_1 = [data['adiabatic_energies']['E_1'] for data in total_data]
data_2 = [data['adiabatic_energies']['E_2'] for data in total_data]
data_3 = [data['adiabatic_energies']['E_3'] for data in total_data]
data_4 = [data['adiabatic_energies']['E_4'] for data in total_data]

data_1, data_2, data_3, data_4 = correct_order_list([data_1, data_2, data_3, data_4], states_orders)

f = multiplot([data_1, data_2,  data_3,  data_4,],
              #['$1 - {}^1A_g$', '$1 - {}^1B_{3u}$', '$2- {}^1B_{3u}$', '$2 - {}^1A_g$'],
              ['${}^1A_g$', '${}^1B_{3u}$', None, None],
              #['$E_1$', '$E_2$', '$E_3$' , '$E_4$'],
              d_coordinate, # title='Adiabatic energies',
              range_y=[6, 13], show_plots=args.show_plots,
              colors=['black', 'red', 'black', 'red'],
              )

f.savefig(folder + "figure2a.pdf", bbox_inches='tight')

e1 = np.array(data_1)
e2 = np.array(data_2)
e3 = np.array(data_3)
e4 = np.array(data_4)


f = multiplot([data_1, data_2,  data_3,  data_4, e_le, e_ct],
              #['$1 - {}^1A_g$', '$1 - {}^1B_{3u}$', '$2- {}^1B_{3u}$', '$2 - {}^1A_g$'],
              ['${}^1A_g$', '${}^1B_{3u}$', None, None, '$E_{LE}$', '$E_{CT}$'],
              #['$E_1$', '$E_2$', '$E_3$' , '$E_4$'],
              d_coordinate, # title='Adiabatic energies',
              range_y=[6, 13], show_plots=args.show_plots,
              colors=['black', 'red', 'black', 'red', None, None],
              )

f.savefig(folder + "figure2_mix.pdf", bbox_inches='tight')



#######################  Figure 3  #####################

o_sx_1 = o_e1 + o_h1
o_sx_2 = o_e2 + o_h2

di_1 = l1*(e_ct - e_le) + d2le1*de_le + d2ct1*de_ct
di_2 = l2*(e_ct - e_le) + d2le2*de_le + d2ct2*de_ct


f = multiplot_2axis([o_dc1, o_sx_1, o_ct1, di_1, -j_coul], ['$\Omega^{(1)}_{DC}$', '$\Omega_{SX}^{(1)}$',
                                                            '$\Omega_{CT}^{(1)}$', '$\Delta_1$', '$J_{Coul}$'],
                    [l1], ['$\lambda_1^2$'], d_coordinate, ylabel2='$\lambda^2$',
                    colors=[None, None, None, None, '#1f77b4'], style1=['-', '-', '-', '-', '--'],
                    range_y2=[0, 0.2], title='State 1', colors2=['grey'], baseline=0)
f.savefig(folder + "figure3a.pdf", bbox_inches='tight')


f = multiplot_2axis([o_dc2, o_sx_2, o_ct2, di_2, j_coul], ['$\Omega^{(2)}_{DC}$', '$\Omega_{SX}^{(2)}$',
                                                           '$\Omega_{CT}^{(2)}$', '$\Delta_2$', '$J_{Coul}$'],
                    [l2], ['$\lambda_2^2$'], d_coordinate, ylabel2='$\lambda^2$',
                    colors=[None, None, None, None, '#1f77b4'], style1=['-', '-', '-', '-','--'],
                    range_y2=[0, 0.2], title='State 2', colors2=['grey'], baseline=0)

f.savefig(folder + "figure3b.pdf", bbox_inches='tight')


f = multiplot_2axis([oe_1, oh_1, o_sx_1], ['$\Omega_e^{(1)}$', '$\Omega_h^{(1)}$', '$\Omega_{SX}^{(1)}$'],
                    [l1], ['$\lambda_1^2$'], d_coordinate, ylabel2='$\lambda^2$',
                    range_y2=[0.0, 0.2], title='State 1', colors2=['grey'], baseline=0)

f.savefig(folder + "figure3c.pdf", bbox_inches='tight')

f = multiplot_2axis([oe_2, oh_2, o_sx_1], ['$\Omega_e^{(2)}$', '$\Omega_h^{(2)}$', '$\Omega_{SX}^{(2)}$'],
                    [l2], ['$\lambda_2^2$'], d_coordinate, ylabel2='$\lambda^2$',
                    range_y2=[0.0, 0.2], title='State 2', colors2=['grey'], baseline=0)

f.savefig(folder + "figure3d.pdf", bbox_inches='tight')



# Figure omega

data_1 = [data['diabatic_contributions']['Omega_One_dc'][0] for data in total_data]
data_2 = [data['diabatic_contributions']['Omega_One_dc'][0] for data in total_data]
data_3 = [data['diabatic_contributions']['Omega_One_dc'][0] for data in total_data]
data_4 = [data['diabatic_contributions']['Omega_One_dc'][0] for data in total_data]

data_1, data_2, data_3, data_4 = correct_order_list([data_1, data_2, data_3, data_4], states_orders)

o_one_dc = data_1

data_1 = [data['diabatic_contributions']['Omega_J_dc'][0] for data in total_data]
data_2 = [data['diabatic_contributions']['Omega_J_dc'][0] for data in total_data]
data_3 = [data['diabatic_contributions']['Omega_J_dc'][0] for data in total_data]
data_4 = [data['diabatic_contributions']['Omega_J_dc'][0] for data in total_data]

data_1, data_2, data_3, data_4 = correct_order_list([data_1, data_2, data_3, data_4], states_orders)

o_j_dc = data_1

data_1 = [data['diabatic_contributions']['Omega_K_dc'][0] for data in total_data]
data_2 = [data['diabatic_contributions']['Omega_K_dc'][0] for data in total_data]
data_3 = [data['diabatic_contributions']['Omega_K_dc'][0] for data in total_data]
data_4 = [data['diabatic_contributions']['Omega_K_dc'][0] for data in total_data]

data_1, data_2, data_3, data_4 = correct_order_list([data_1, data_2, data_3, data_4], states_orders)

o_k_dc = data_1

f = multiplot_2axis([o_dc1, o_one_dc, o_j_dc,  o_k_dc, -j_coul],
                    ['$\Omega_{DC}^{(1)}$', '$\Omega_{DC,1}^{(1)}$', '$\Omega_{DC,J}^{(1)}$', '$\Omega_{DC,K}^{(1)}$', '$J_{Coul}$'],
                    [l1], ['$\lambda_1^2$'], d_coordinate, ylabel2='$\lambda^2$',
                    colors=[None, None, None, None, '#2ca02c'], style1=['-', '-', '-', '-', '--'],
                    range_y2=[0, 0.2], colors2=['grey'], baseline=0)

f.savefig(folder + "omega_other.pdf", bbox_inches='tight')

exit()

f = multiplot_2axis([oe_1, oh_1], ['$\Omega_e^{(1)}$', '$\Omega_h^{(1)}$'],
                    [np.square(l1)], ['$\lambda_1^2$'], d_coordinate, ylabel2='$\lambda^2$',
                    style1=['-', '-.'], colors=['black', 'black'],
                    range_y2=[-0.3, 0.3], title='State 1', colors2=['grey'], baseline=0)

f.savefig(folder + "figure3c_x.pdf", bbox_inches='tight')

f = multiplot_2axis([oe_2, oh_2], ['$\Omega_e^{(2)}$', '$\Omega_h^{(2)}$'],
                    [np.square(l2)], ['$\lambda_2^2$'], d_coordinate, ylabel2='$\lambda^2$',
                    style1=['-', '-.'], colors=['red', 'red'],
                    range_y2=[-0.3, 0.3], title='State 2', colors2=['grey'], baseline=0)

f.savefig(folder + "figure2d_x.pdf", bbox_inches='tight')



exit()

#######################  adiabatic_energies (calculated) ######################

e1c = e_le + wdc1 + e_11 + e_21
e2c = e_le + wdc2 + e_12 + e_22

f = multiplot([e1-e1c, e2-e2c],
              ['diff 1', 'diff 2'],
              d_coordinate, title='Adiabatic energies difference\n (original-calculated)',
              range_y=None, show_plots=args.show_plots)

f.savefig(folder + "test.pdf", bbox_inches='tight')

#######################  Figure 3  #####################

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


import numpy as np
import matplotlib.pyplot as plt


def plot_obtention(input, x_data, y_data):
    # "matplotlib" help: https://aprendeconalf.es/docencia/python/manual/matplotlib/
    # https://chartio.com/resources/tutorials/how-to-save-a-plot-to-a-file-using-matplotlib/
    # https://pythonspot.com/matplotlib-bar-chart/

    fuente = "Sans"
    size = 12

    y_pos = np.arange(len(x_data))
    plt.bar(y_pos, y_data, align='center', alpha=0.5, color='red')
    plt.xticks(y_pos, x_data)

    plt.title(input, fontsize=size, fontname=fuente)
    plt.ylabel('Energies (eV)', fontsize=size, fontname=fuente)
    plt.xlabel('number of state', fontsize=size, fontname=fuente)
    plt.axis([min(x_data) - 2, max(x_data), min(y_data), max(y_data) + 0.1 * max(y_data)])
    # plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    # plt.grid(True)

    plt.plot()
    figure_name = input + '.png'
    plt.savefig(figure_name)

    # plt.show()
    plt.close()


def get_bar_chart(input, x_list, y_list, x_title, y_title, main_title):
    y_pos = (list(x_list))
    plt.bar(x_list, y_list, align='center', width=0.5, color='r', edgecolor="black")
    plt.xticks(y_pos)

    fuente = 'serif'

    plt.xlabel(x_title, fontsize=14, fontfamily=fuente)
    plt.ylabel(y_title, fontsize=14, fontfamily=fuente)
    plt.title(main_title, fontsize=18, fontfamily=fuente)

    plt.plot()
    figure_name = input + '_' + main_title + '.png'
    plt.savefig(figure_name)
    plt.show()
    plt.close()


def plot_g_tensor_vs_states(presentation_matrix, x_title, y_title, main_title, save_picture):
    plt.plot(presentation_matrix[:, 0], presentation_matrix[:, 1], 'r', label='gxx')
    plt.plot(presentation_matrix[:, 0], presentation_matrix[:, 2], 'b', label='gyy')
    plt.plot(presentation_matrix[:, 0], presentation_matrix[:, 3], 'k', label='gzz')

    # MARKER TYPES: https://matplotlib.org/2.1.1/api/_as_gen/matplotlib.pyplot.plot.html
    plt.plot(presentation_matrix[:, 0], presentation_matrix[:, 1], 'ro') #, label='gxx')
    plt.plot(presentation_matrix[:, 0], presentation_matrix[:, 2], 'bo') #, label='gyy')
    plt.plot(presentation_matrix[:, 0], presentation_matrix[:, 3], 'ks') #, label='gzz')

    fuente = 'serif'

    plt.xlabel(x_title, fontsize=20, fontfamily=fuente)
    plt.ylabel(y_title, fontsize=20, fontfamily=fuente)
    plt.title(main_title, fontsize=26, fontfamily=fuente)

    # plt.locator_params(nbins=10)

    plt.grid()
    plt.legend()

    if (save_picture == 0):
        plt.show()
        plt.close()
    else:
        figure_name = main_title + '_sos_analysis.png'
        plt.savefig(figure_name)

def sos_analysis_and_plot(input, save_picture):
    """"
    PROGRAM TO CALCULATE THE G-TENSOR OVER FROM
    AN INITIAL TO A FINAL NUMBER OF STATES IN THE
    SUM-OVER-STATES EXPANSION
    """
    from g_read import get_number_of_states, get_eigenenergies, get_selected_states,get_spin_orbit_couplings,\
        get_spin_matrices, get_orbital_matrices, get_symmetry_states
    from g_operations import from_energies_SOC_to_g_values

    totalstates = get_number_of_states(input)

    presentation_list = []
    # presentation_list.append(['State', 'gxx', 'gyy', 'gzz'])

    for i in range(1, totalstates + 1):
        states_ras = list(range(1, i + 1))

        eigenenergies_ras, excitation_energies_ras = get_eigenenergies(input, totalstates, states_ras)

        SOC_ras = get_spin_orbit_couplings(input, totalstates, states_ras, selected_SOC=0)

        ras_G_matrix, ras_g_values, eigenvalues, eigenvector = from_energies_SOC_to_g_values(input, states_ras, totalstates,
                                                                   excitation_energies_ras, SOC_ras)

        state_symmetries, ordered_state_symmetries = get_symmetry_states(input, totalstates)

        # presentation_list.append([ordered_state_symmetries[i-1], np.round(ras_g_values.real[0], 3), np.round(ras_g_values.real[1], 3),
        #                           np.round(ras_g_values.real[2], 3)])
        presentation_list.append([i, np.round(ras_g_values.real[0], 3), np.round(ras_g_values.real[1], 3),
                                  np.round(ras_g_values.real[2], 3)])

    presentation_matrix = np.array(presentation_list, dtype=object)

    # To presents deviation from previous g-values instead of the total g-values:
    presentation_matrix_deviation = np.array(presentation_list, dtype=object)
    for ndim in [1,2,3]:
        for i in range(1, len(presentation_matrix)):
            presentation_matrix_deviation[i,ndim] = (presentation_matrix[i,ndim] - presentation_matrix[i-1,ndim])

    print("--------------------------------")
    print(" SUM-OVER-STATE ANALYSIS")
    print("--------------------------------")
    print('\n'.join([''.join(['{:^20}'.format(item) for item in row]) for row in (presentation_matrix[:, :])]))

    plot_g_tensor_vs_states(presentation_matrix_deviation, x_title='State', y_title='g-values deviations (ppt)',
                            main_title='g-tensor sum-over-states analysis', save_picture=0)

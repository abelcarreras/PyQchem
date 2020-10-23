import numpy as np
import matplotlib.pyplot as plt


def draw_electron(pos, spin=1, color='black'):
    a = 0.5
    s = 0.05
    if spin == 1:
        d = 0.1
        plt.plot([pos[0]+d, pos[0]+0.1+d], [pos[1]+a, pos[1]+a*0.6], 'k-', lw=2, color=color)
    else:
        d = -0.1
        plt.plot([pos[0]+d, pos[0]-0.1+d], [pos[1]+s, pos[1]+a*0.4+s], 'k-', lw=2, color=color)

    plt.plot([pos[0]+d, pos[0]+d], [pos[1]+s, pos[1] + a], 'k-', lw=2, color=color)


def plot_configuration(alpha, beta, index=1):

    levels_alpha = range(len(alpha))

    plt.grid(axis='y')

    x = [index + 1] * len(levels_alpha)

    plt.annotate(index+1, xy=(index+1, 0), xytext=(0, -4), size=10,
                 ha="center", va="top", textcoords="offset points")

    for xi, yi, oa, ob in zip(x, levels_alpha, alpha, beta):

        plt.plot([xi-0.3, xi+0.3], [yi, yi], 'k-', lw=3, color='blue')
        if oa == '1':
            draw_electron([xi, yi], spin=1, color='r')
        if ob == '1':
            draw_electron([xi, yi], spin=-1, color='r')

    plt.margins(0.2)


def plot_state(state, with_amplitude=False, orbital_range=None):

    plt.figure(figsize=(len(state['configurations']), 5))

    n_orbitals = 1
    amplitude_list = []
    for j, conf in enumerate(state['configurations']):
        # print(conf)
        alpha = ''.join(np.array(conf['occupations']['alpha'], dtype=str))
        beta = ''.join(np.array(conf['occupations']['beta'], dtype=str))

        if orbital_range is not None:
            alpha = alpha[orbital_range[0]: orbital_range[1]]
            beta = beta[orbital_range[0]: orbital_range[1]]

        plot_configuration(alpha, beta, index=j)
        amplitude_list.append(conf['amplitude'])

        n_orbitals = len(alpha)

    #plt.axis('off')
    frame1 = plt.gca()
    frame1.axes.get_yaxis().set_visible(False)

    frame1.axes.xaxis.set_ticklabels([])
    #frame1.axes.yaxis.set_ticklabels([])

    for tick in frame1.axes.get_xticklines():
        tick.set_visible(False)
    for tick in frame1.axes.get_yticklines():
        tick.set_visible(False)

    plt.xlabel('Configurations')


    if with_amplitude:
        plt.plot(range(1, len(amplitude_list) + 1), np.square(amplitude_list) * n_orbitals,
                 label='amplitudes', color='green')
        plt.ylabel('Amplitude')
        plt.legend()


def plot_diabatization(diabatic_states, atoms_ranges=(1, 2)):

    for i, state in enumerate(diabatic_states):

        bars = range(1, atoms_ranges[-1]+1)

        for pos in atoms_ranges[:-1]:
            plt.axvline((pos+0.5), color='black')

        plt.figure(i+1)
        plt.suptitle('Mulliken analysis')
        plt.title('Diabatic state {}'.format(i+1))
        plt.bar(bars, state['mulliken']['attach'], label='Attachment')
        plt.bar(bars, state['mulliken']['detach'], label='Detachment')
        plt.plot(bars, state['mulliken']['total'], label='Total', color='r')
        plt.xlabel('Atoms')
        plt.ylabel('Charge [e-]')
        plt.legend()

    plt.show()

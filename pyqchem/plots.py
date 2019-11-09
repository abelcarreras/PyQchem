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


def plot_state(alpha, beta, index=1):

    levels_alpha = range(len(alpha))

    plt.grid(axis='y')

    x = [index + 1] * len(levels_alpha)

    plt.annotate(index+1, xy=(index+1, 0), xytext=(0, -4), size=10,
                 ha="center", va="top", textcoords="offset points")

    for xi, yi, oa, ob in zip(x,levels_alpha, alpha, beta):

        plt.plot([xi-0.3, xi+0.3], [yi, yi], 'k-', lw=3, color='blue')
        if oa == '1':
            draw_electron([xi, yi], spin=1, color='r')
        if ob == '1':
            draw_electron([xi, yi], spin=-1, color='r')

    plt.margins(0.2)



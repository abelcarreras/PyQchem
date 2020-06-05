import numpy as np
from pyqchem.units import DEBYE_TO_AU
import matplotlib.pyplot as plt
from pyqchem.plots import plot_state


def print_excited_states(parsed_data, include_conf_rasci=False, include_mulliken_rasci=False):
    for i, state in enumerate(parsed_data):
        print('\nState {}'.format(i+1))
        if 'multiplicity' in state:
            print('Multiplicity', state['multiplicity'])

        if state['transition_moment'] is not None:
            print('Osc. strength: {:6.4f}'.format(state['oscillator_strength']))
            print('Transition DM: ', '{:6.4f} {:6.4f} {:6.4f}'.format(*state['transition_moment']),
                  ' ', state['dipole_moment_units'])

        print('Energy: ', state['excitation_energy'], ' ', state['excitation_energy_units'])

        if include_conf_rasci:
            print('    Configurations')
            print('  Hole   Alpha  Beta   Part   Amplitude')
            for j, conf in enumerate(state['configurations']):
                print(' {:^6} {:^6} {:^6} {:^6} {:8.3f}'.format(conf['hole'], conf['alpha'], conf['beta'], conf['part'], conf['amplitude']))

        if include_mulliken_rasci:
            print('Mulliken analysis')
            print('         Attach    Detach    Total ')
            for i_atom, (at, det) in enumerate(zip(state['mulliken']['attach'],
                                                   state['mulliken']['detach'])):
                print('{:5}  {:8.4f}  {:8.4f}  {:8.4f}'.format(i_atom + 1, at, det, at + det))


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


def plot_rasci_state_configurations(states):
    for i, state in enumerate(states):
        plt.figure(figsize=(len(state['configurations']), 5))
        plt.title('State {}'.format(i+1))
        amplitude_list = []
        for j, conf in enumerate(state['configurations']):
            plot_state(conf['alpha'], conf['beta'], index=j)
            amplitude_list.append(conf['amplitude'])

        plt.plot(range(1, len(amplitude_list)+1), np.square(amplitude_list)*len(state['configurations'][0]['alpha']), label='amplitudes')
        plt.xlabel('Configurations')
        plt.ylabel('Amplitude')
        plt.axis('off')
        plt.legend()

    plt.show()

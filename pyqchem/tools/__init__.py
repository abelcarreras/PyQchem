import matplotlib.pyplot as plt
from pyqchem.plots import plot_configuration
from pyqchem.structure import Structure
from urllib.request import urlopen
import requests as req
import numpy as np
import json
import warnings
from requests.exceptions import ConnectionError


def print_excited_states(parsed_data, include_conf_rasci=False, include_mulliken_rasci=False):
    for i, state in enumerate(parsed_data):
        print('\nState {}'.format(i+1))
        if 'multiplicity' in state:
            print('Multiplicity', state['multiplicity'])

        if state['transition_moment'] is not None:
            #print('Osc. strength: {:6.4f}'.format(state['oscillator_strength']))
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


def plot_rasci_state_configurations(states):
    for i, state in enumerate(states):
        plt.figure(figsize=(len(state['configurations']), 5))
        plt.title('State {}'.format(i+1))
        amplitude_list = []
        for j, conf in enumerate(state['configurations']):
            plot_configuration(conf['alpha'], conf['beta'], index=j)
            amplitude_list.append(conf['amplitude'])

        plt.plot(range(1, len(amplitude_list)+1), np.square(amplitude_list)*len(state['configurations'][0]['alpha']), label='amplitudes')
        plt.xlabel('Configurations')
        plt.ylabel('Amplitude')
        plt.axis('off')
        plt.legend()

    plt.show()


def submit_notice(message, service='pushbullet', pb_token=None, sp_url=None):
    """
    Submit a notification using webhooks

    :param message: The message to send
    :param service: pushbullet or samepage
    :param pb_token: pushbullet token
    :param sp_url: samepage url
    :return:
    """

    if service.lower() == 'pushbullet':
        url = 'https://api.pushbullet.com/v2/pushes'
        bot_message = {'body': message, 'type': 'note'}
        message_headers = {'Content-Type': 'application/json; charset=UTF-8',
                           'Access-Token': pb_token}
    elif service.lower() == 'samepage':
        url = sp_url
        bot_message = {'text': message}
        message_headers = {'Content-Type': 'application/json'}
    else:
        print('client not found!')
        return

    try:
        r = req.post(url=url,
                     headers=message_headers,
                     data=json.dumps(bot_message))
        r.close()

    except ConnectionError:
        warnings.warn('Connection error: Message was not delivered')


def rotate_coordinates(coordinates, angle, axis, atoms_list=None, center=(0, 0, 0)):
    """
    Rotate the coordinates (or range of coordinates) with respect a given axis

    :param coordinates: coordinates to rotate
    :param angle: rotation angle
    :param axis: rotation axis
    :param atoms_list: list of atoms to rotate (if None then rotate all)
    :return: rotated coordinates
    """
    coordinates = np.array(coordinates) - np.array(center)

    cos_term = 1 - np.cos(angle)
    rot_matrix = [[axis[0]**2*cos_term + np.cos(angle), axis[0]*axis[1]*cos_term - axis[2]*np.sin(angle), axis[0]*axis[2]*cos_term + axis[1]*np.sin(angle)],
                  [axis[1]*axis[0]*cos_term + axis[2]*np.sin(angle), axis[1]**2*cos_term + np.cos(angle), axis[1]*axis[2]*cos_term + axis[0]*np.sin(angle)],
                  [axis[2]*axis[0]*cos_term + axis[1]*np.sin(angle), axis[1]*axis[2]*cos_term + axis[0]*np.sin(angle), axis[2]**2*cos_term + np.cos(angle)]]

    if atoms_list is not None:
        coordinates[atoms_list] = np.dot(coordinates[atoms_list], rot_matrix)
    else:
        coordinates = np.dot(coordinates, rot_matrix) + np.array(center)

    return coordinates.tolist()


def get_geometry_from_pubchem(entry, type='name'):
    """
    Get structure form PubChem database

    :param entry: entry data
    :param type: data type: 'name', 'cid'
    :return: Structure
    """

    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
    input_1 = 'compound/{}/'.format(type)
    output_format = "JSON"
    additional = "record_type=3d"

    apiurl = base + input_1 + output_format + '?' + additional
    postdata = '{}={}'.format(type, entry).encode()

    from urllib.error import HTTPError

    try:
        response = urlopen(apiurl, postdata)
    except HTTPError as e:
        string = e.read().decode("utf-8")
        json_data = json.loads(string)
        fault = json_data['Fault']
        if 'Details' in fault:
            raise Exception(fault['Details'][0])
        else:
            raise Exception(fault['Message'])

    string = response.read().decode("utf-8")
    json_data = json.loads(string)

    conformers = json_data['PC_Compounds'][0]['coords'][0]['conformers'][0]
    atoms = json_data['PC_Compounds'][0]['atoms']

    positions = np.array([conformers['x'], conformers['y'], conformers['z']]).T
    atomic_numbers = atoms['element']

    if 'charge' in atoms:
        charge = np.add.reduce([c_atom['value'] for c_atom in atoms['charge']])
    else:
        charge = 0

    return Structure(coordinates=positions,
                     atomic_numbers=atomic_numbers,
                     charge=charge)


if __name__ == '__main__':

    mol = get_geometry_from_pubchem('acetone', type='name')
    print(mol)

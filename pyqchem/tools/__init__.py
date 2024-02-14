import matplotlib.pyplot as plt
from pyqchem.plots import plot_configuration
from pyqchem.structure import Structure
from urllib.request import urlopen
import requests as req
import numpy as np
import json
import warnings
from requests.exceptions import ConnectionError
from pyqchem.tools.geometry import rotate_coordinates


def print_excited_states(parsed_data, include_conf_rasci=False, include_mulliken_rasci=False):
    """
    Prints excited states in nice format. It works for CIS/TDDFT/RASCI methods

    :param parsed_data: parsed data (excited states) dictionary entry from CIS/RASCI/TDDFT calculation
    :param include_conf_rasci: print also configuration data (only RASCI method)
    :param include_mulliken_rasci: print also mulliken analysis (only RASCI method)
    :return: None
    """
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
    """
    Prints
    :param states: parsed data (excited states) dictionary entry from RASCI calculation
    :return: None
    """
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


def submit_notice(message,
                  service='pushbullet',
                  pb_token=None,
                  sp_url=None,
                  gc_key=None, gc_token=None, gc_thread=None,
                  slack_token=None, slack_channel=None, ):
    """
    Submit a notification using webhooks

    :param message: The message to send
    :param service: pushbullet, samepage, google_chat
    :param pb_token: pushbullet token
    :param sp_url: samepage url
    :param gc_key: google chat key
    :param gc_token: google chat token
    :param gc_thread: google chat thread
    :param slack_token: slack bot token (xoxb-xxx.xxx.xxx),
    :param slack_channel: slack channel
    :return: server response
    """

    if service.lower() == 'pushbullet':
        if pb_token is None:
            warnings.warn('Message: you have to specify {}'.format(['pb_token']))
        url = 'https://api.pushbullet.com/v2/pushes'
        bot_message = json.dumps({'body': message, 'type': 'note'})
        message_headers = {'Content-Type': 'application/json; charset=UTF-8',
                           'Access-Token': pb_token}
    elif service.lower() == 'samepage':
        if sp_url is None:
            warnings.warn('You have to specify {}'.format(['sp_url']))

        url = sp_url
        bot_message = json.dumps({'text': message})
        message_headers = {'Content-Type': 'application/json'}

    elif service.lower() == 'google_chat':
        if gc_key is None or gc_token is None or gc_thread is None:
            warnings.warn('Message: you have to specify {}'.format(['gc_key', 'gc_token', 'gc_token']))

        url = 'https://chat.googleapis.com/v1/spaces/<space>/messages?key={}\&token={}\&threadKey={}'.format(gc_key, gc_token, gc_thread)
        bot_message = {'text': message}
        message_headers = {'Content-Type': 'application/json'}

    elif service.lower() == 'slack':
        if slack_token is None or slack_channel is None:
            warnings.warn('Message: you have to specify {}'.format(['slack_token', 'slack_channel']))

        url = 'https://slack.com/api/chat.postMessage'
        blocks = None
        bot_message = {'token': slack_token,
                       'channel': slack_channel,
                       'text': message,
                       'blocks': json.dumps(blocks) if blocks else None}
        message_headers = {}

    else:
        warnings.warn('Message: client not found!')
        return

    try:
        with req.post(url=url, headers=message_headers, data=bot_message) as r:
            response = r.json()

        return response

    except ConnectionError:
        warnings.warn('Connection error: Message was not delivered')


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
                     charge=charge,
                     name=str(entry))


if __name__ == '__main__':

    mol = get_geometry_from_pubchem('acetone', type='name')
    print(mol)

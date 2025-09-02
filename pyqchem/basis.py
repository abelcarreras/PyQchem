import requests as req
import unicodedata
import re
import numpy as np
import six
from pyqchem.structure import atom_data
from copy import deepcopy


def _txt_to_basis_dict(basis_txt):
    """
    read basis in gaussian/qchem format

    :param basis_txt: string containing the basis set in gaussian/qchem format
    :return: basis set dicitonary
    """

    symbol = basis_txt[0].split()[0]

    def is_number(s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    basis_pure = basis_txt[1:]

    section_marks = []
    for i, line in enumerate(basis_pure):
        if not is_number(line.split()[0]):
            section_marks.append(i)

    shells = []
    for i in section_marks[:-1]:
        type, n_func, _  = basis_pure[i].split()
        n_func = int(n_func)

        if type.upper() in ['SP']:
            p_exponent, con_coefficients, p_con_coefficients = np.array([line.split()
                                                                         for line in basis_pure[i + 1:i + n_func + 1]],
                                                                        dtype=float).T
        else:
            p_exponent, con_coefficients = np.array([line.split()
                                                     for line in basis_pure[i + 1:i + n_func + 1]],
                                                    dtype=float).T
            p_con_coefficients = np.zeros_like(p_exponent)


        shells.append({'shell_type': type,
                       'p_exponents': list(p_exponent),
                       'con_coefficients': list(con_coefficients),
                       'p_con_coefficients': list(p_con_coefficients)})

    return {'symbol': symbol,
            'shells': shells}


def get_basis_element_from_ccRepo(element,
                                  program='Gaussian',
                                  basis='cc-pVDZ'):
    """
    Downloads basis set for a particular atomic element from ccRepo database

    :param element: atomic element
    :param program: format of basis set
    :param basis: basis name
    :return: citation, description, basis set dictionary
    """
    from lxml import html

    # Check main page element list
    with req.get("http://www.grant-hill.group.shef.ac.uk/ccrepo/") as resp:
        element_list = []
        for m in re.finditer('class="sym" href=', resp.text):
            # print(resp.text[m.end():m.end()+50].split('"')[1] )
            element_list.append(resp.text[m.end():m.end()+50].split('"')[1])

    element_dict = {atom[1]: atom[2].lower() for atom in atom_data[1:]}

    try:
        resp = req.get("http://www.grant-hill.group.shef.ac.uk/ccrepo/{}".format(element_dict[element]))
    except KeyError as e:
        raise Exception('Atom label {} not recognized'.format(e.args[0]))

    n_ini = resp.text.find('form-inline')

    # determinte the php web address for particular element
    php_site = resp.text[n_ini:n_ini+100].split('"')[2]

    # print('http://www.grant-hill.group.shef.ac.uk/ccrepo/{}/{}'.format(element, php_site))
    # r = req.post(url = 'http://www.grant-hill.group.shef.ac.uk/ccrepo/titanium/tibasis.php', data = data)
    r = req.post(url='http://www.grant-hill.group.shef.ac.uk/ccrepo/{}/{}'.format(element_dict[element], php_site),
                 data={'basis': basis, 'program': program})

    # parse the data
    tree = html.fromstring(r.content)
    r.close()

    header = tree.xpath('//div[@class="container"]/text()')
    citation = unicodedata.normalize('NFC', six.text_type(header[1])).strip()
    description = unicodedata.normalize('NFC', header[2]).strip()


    basis_data = tree.xpath('/html/body/div/nobr/text()')
    basis_clean = [unicodedata.normalize('NFKC', six.text_type(line)).strip() for line in basis_data]
    # basis_name = basis_clean[1].split('"')[1]

    return citation, description, basis_clean[2:]


def get_basis_from_ccRepo(structure, basis, full=False, if_missing=None):
    """
    Get basis from ccRepo

    :param structure: Structure
    :param basis: basis set label (string)
    :param full: if False only list basis for unique atoms
    :param if_missing: backup basis to use if basis is missing for a particular atom
    :return: basis set dictionary
    """

    symbols = structure.get_symbols()
    if not full:
        symbols = np.unique(symbols)

    atoms = []
    for symbol in symbols:
        citation, description, basis_data = get_basis_element_from_ccRepo(symbol,
                                                                          program='Gaussian',
                                                                          basis=basis)

        if len(basis_data) == 0:
            if if_missing is not None:
                citation, description, basis_data = get_basis_element_from_ccRepo(symbol,
                                                                                  program='Gaussian',
                                                                                  basis=if_missing)
                if len(basis_data) == 0:
                    raise Exception('Basis {}, {} not found for atom {} in ccRepo'.format(basis, if_missing, symbol))
            else:
                raise Exception('Basis {} not found for atom {} in ccRepo'.format(basis, symbol))

        # print(description)
        # print(citation)
        atoms.append(_txt_to_basis_dict(basis_data))

    basis_set = {'name': basis,
                 'primitive_type': 'gaussian',
                 'atoms': atoms}

    return basis_set


def get_basis_from_BSE(structure, basis, full=False, if_missing=None):
    """
    get basis from Basis Set Exchange

    :param structure: Structure object (the molecule)
    :param basis: basis set label (string)
    :param full: if False only list basis for unique atoms
    :return: basis set dictionary
    """

    base_url = "http://basissetexchange.org"
    basis = basis.lower()

    headers = {
        'User-Agent': 'BSE Example Python Script',
        'From': 'bse@molssi.org'
    }

    r = req.get(base_url + '/api/basis/{}/format/gaussian94'.format(basis),
                     headers=headers
                     )

    description = '\n'.join([line for line in r.text.split('\n') if '!'  in line])

    r = req.get(base_url + '/api/references/{}/format/bib'.format(basis),
                     headers=headers
                     )

    citation = r.text

    symbols = structure.get_symbols()
    if not full:
        symbols = np.unique(symbols)

    atoms = []
    for symbol in symbols:

        params = {'elements': [symbol]}
        r = req.get(base_url + '/api/basis/{}/format/gaussian94'.format(basis),
                    params=params,
                    headers=headers
                    )
        # https://www.basissetexchange.org/basis/coemd-2/format/gaussian94/?version=0&elements=6,7,8
        if r.status_code != 200:
            if if_missing is not None:
                r = req.get(base_url + '/api/basis/{}/format/gaussian94'.format(if_missing),
                            params=params,
                            headers=headers
                            )
                if r.status_code != 200:
                    raise Exception('Basis {}, {} not found for atom {} in BSE'.format(basis, if_missing, symbol))
            else:
                #raise RuntimeError("Could not obtain data from the BSE. Check the error information above")
                raise Exception('Basis {} not found for atom {} in BSE'.format(basis, symbol))

        basis_data = []
        for line in r.text.split('\n'):
            if len(line) !=0 and line[0] not in ['!']:
                basis_data.append(line.replace('D+', 'E+').replace('D-', 'E-'))


        atoms.append(_txt_to_basis_dict(basis_data))

    basis_set = {'name': basis,
                 'primitive_type': 'gaussian',
                 'atoms': atoms}

    return basis_set


def basis_to_txt(basis):
    """
    convert from basis set dictionary to string in Gaussian/Q-Chem format

    :param basis: basis set dictionary
    :return: basis set string in gaussian/Q-Chem format
    """

    basis_txt = ''
    for atom in basis['atoms']:
        basis_txt += atom['symbol'] + '  0\n'
        for shell in atom['shells']:
            # check if shell is pure or cartesian
            if shell['shell_type'].endswith('_'):
                shell_type = shell['shell_type'][:-1]
            else:
                shell_type = shell['shell_type']

            basis_txt += '{} {} {}\n'.format(shell_type.upper(), len(shell['p_exponents']), 1.00)
            for p, c, pc in zip(shell['p_exponents'], shell['con_coefficients'], shell['p_con_coefficients']):
                if shell['shell_type'].upper() in ['SP']:
                    basis_txt += '{:15.10e} {:15.10e} {:15.10e} \n'.format(p, c, pc)
                else:
                    basis_txt += '{:15.10e} {:15.10e} \n'.format(p, c)

        basis_txt += '****\n'
    return basis_txt


def get_purecard(basis):
    """
    returns a string indicating for each shell wether is pure or Cartesian.
    To be used in Purecard keyword in Q-Chem input

    :param basis: basis set dictionary
    :return: string for purecard keyword
    """

    # default 2111
    keyword = {'d': 1, 'f': 1, 'g': 1, 'h': 2}

    # check if basis is pure
    for atom in basis['atoms']:
        for shell in atom['shells']:
            if shell['shell_type'].endswith('_'):
                keyword[shell['shell_type'][:-1].lower()] = 1
            else:
                keyword[shell['shell_type'].lower()] = 2

    return '{}{}{}{}'.format(keyword['h'], keyword['g'], keyword['f'], keyword['d'])


def trucate_basis(basis, shells=()):
    """
    delete shells from the basis to reduce its size

    :param basis: basis dictionary
    :param shells: shell type list to be removed ex: ('S', 'P', 'D', 'H', 'I', ..)
    :return: truncated basis
    """

    basis_trunc = deepcopy(basis)
    for check_shell in shells:
        for atom in basis_trunc['atoms']:
            delete_list = []
            for i, shell in enumerate(atom['shells']):
                if shell['shell_type'].upper() == check_shell.upper():
                    del shell
                    delete_list.append(i)
            for i in sorted(delete_list, reverse=True):
                del atom['shells'][i]

    return basis_trunc

if __name__ == '__main__':

    citation, description, basis = get_basis_element_from_ccRepo('C',
                                                                 program='Gaussian',
                                                                 basis='cc-pVTZ')

    print(citation)
    print(description)
    print('-----------------------')
    for line in basis:
        print(line)

    basis = _txt_to_basis_dict(basis)
    print(basis)


    from pyqchem.structure import Structure


    # create molecule
    molecule = Structure(coordinates=[[0.0, 0.0, 0.0000],
                                      [0.0, 0.0, 1.5811]],
                         symbols=['Se', 'H'],
                         charge=-1,
                         multiplicity=1)

    basis = get_basis_from_ccRepo(molecule,
                                  basis='cc-pVTZ',
                                  full=False)

    print(basis_to_txt(basis))

    print('----------------------')
    exit()
    basis = get_basis_from_BSE(molecule,
                               basis='cc-pVTZ',
                               full=False)

    print(basis_to_txt(basis))


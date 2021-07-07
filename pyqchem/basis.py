import requests as req
from lxml import html
import unicodedata
import re
import numpy as np
import six
from pyqchem.structure import atom_data
from copy import deepcopy


def _txt_to_basis_dict(basis_txt):
    # read basis in gaussian/qchem format

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

    # Check main page element list
    with req.get("http://www.grant-hill.group.shef.ac.uk/ccrepo/") as resp:
        element_list = []
        for m in re.finditer('class="sym" href=', resp.text):
            # print(resp.text[m.end():m.end()+50].split('"')[1] )
            element_list.append(resp.text[m.end():m.end()+50].split('"')[1])

    element_dict = {atom[1]: atom[2].lower() for atom in atom_data[1:]}

    resp = req.get("http://www.grant-hill.group.shef.ac.uk/ccrepo/{}".format(element_dict[element]))
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
                    raise Exception('Basis {}, {} not found for atom {}'.format(basis, if_missing, symbol))
            else:
                raise Exception('Basis {} not found for atom {}'.format(basis, symbol))

        # print(description)
        # print(citation)
        atoms.append(_txt_to_basis_dict(basis_data))

    basis_set = {'name': basis,
                 'primitive_type': 'gaussian',
                 'atoms': atoms}

    return basis_set


def basis_to_txt(basis):
    # write basis in qchem/gaussian format
    basis_txt = ''

    for atom in basis['atoms']:
        basis_txt += atom['symbol'] + '\n'
        for shell in atom['shells']:
            basis_txt += '{} {} {}\n'.format(shell['shell_type'].upper(), len(shell['p_exponents']), 1.00)
            for p, c, pc in zip(shell['p_exponents'], shell['con_coefficients'], shell['p_con_coefficients']):
                if shell['shell_type'].upper() in ['SP']:
                    basis_txt += '{:15.10e} {:15.10e} {:15.10e} \n'.format(p, c, pc)
                else:
                    basis_txt += '{:15.10e} {:15.10e} \n'.format(p, c)

        basis_txt += '****\n'
    return basis_txt


def trucate_basis(basis, shells=()):

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

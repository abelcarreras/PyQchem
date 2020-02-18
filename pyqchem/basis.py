import requests as req
from lxml import html
import unicodedata
import re


def get_basis_from_ccRepo(element,
                          program='Gaussian',
                          basis='cc-pVDZ'):

    # Check main page element list
    with req.get("http://www.grant-hill.group.shef.ac.uk/ccrepo/") as resp:
        element_list = []
        for m in re.finditer('class="sym" href=', resp.text):
            # print(resp.text[m.end():m.end()+50].split('"')[1] )
            element_list.append(resp.text[m.end():m.end()+50].split('"')[1])

    # define symbol to element dictionary (more should be added)
    element_dict = {'H': 'hydrogen',
                    'C': 'carbon',
                    'N': 'nitrogen',
                    'O': 'oxygen'}

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
    citation = unicodedata.normalize('NFC', header[1]).strip()
    description = unicodedata.normalize('NFC', header[2]).strip()

    basis_data = tree.xpath('/html/body/div/nobr/text()')
    basis_clean = [unicodedata.normalize('NFKC', line).strip() for line in basis_data]

    return citation, description, basis_clean


if __name__ == '__main__':

    citation, description, basis = get_basis_from_ccRepo('O',
                                                         program='Gaussian',
                                                         basis='cc-pVTZ')

    print(citation)
    print(description)
    print('-----------------------')
    for line in basis:
        print(line)
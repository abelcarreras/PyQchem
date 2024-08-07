[![Build, Test & Upload](https://github.com/abelcarreras/PyQchem/actions/workflows/test-publish.yaml/badge.svg)](https://github.com/abelcarreras/PyQchem/actions/workflows/test-publish.yaml)
[![Coverage Status](https://coveralls.io/repos/github/abelcarreras/PyQchem/badge.svg?branch=master)](https://coveralls.io/github/abelcarreras/PyQchem?branch=master)
[![PyPI version](https://badge.fury.io/py/pyqchem.svg)](https://badge.fury.io/py/pyqchem)
[![Documentation Status](https://readthedocs.org/projects/pyqchem/badge/?version=master)](https://pyqchem.readthedocs.io/en/master/?badge=master)

PyQchem
=======
Python wrapper for Q-Chem (https://www.q-chem.com)  
Online manual: https://pyqchem.readthedocs.io/

Main features
-------------
- Easy to use clean python interface for Q-Chem
- No special q-chem installation needed (reads Q-Chem environment)
- Output parsers for different type of calculations 
- Custom basis set and guess support using high level interface
- Calculation Cache system powered by SQLite database
- python 2.7.x/3.5+ compatibility

Installation instructions
-------------------------
1. Requirements

- numpy
- scipy
- matplotlib
- requests
- lxml
- wfnsympy (optional)
- paramiko (optional)
- pymatgen (optional)

2a. From pypi repository (recommended)
```shell
pip install pyqchem --user
```

2b. Installation from source

```shell
python setup.py install --user
```

Examples 
--------
**Simple pythonic API to define your input**

```python
from pyqchem import Structure, QchemInput, get_output_from_qchem
from pyqchem.parsers.basic import basic_parser_qchem

molecule = Structure(coordinates=[[0.0, 0.0, 0.0],
                                  [0.0, 0.0, 0.9]],
                     symbols=['H', 'H'],
                     charge=0,
                     multiplicity=1)

qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      basis='6-31G')

data = get_output_from_qchem(qc_input,
                             processors=4,
                             parser=basic_parser_qchem)

# obtain a python dictionary
print('Energy: ', data['scf_energy'])
```

**Link calculations in powerful workflows**

```python
from pyqchem import QchemInput, get_output_from_qchem
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.parsers.parser_frequencies import basic_frequencies
from pyqchem.tools import get_geometry_from_pubchem

molecule = get_geometry_from_pubchem('aspirin')

qc_input = QchemInput(molecule,
                      jobtype='opt',
                      exchange='hf',
                      basis='sto-3g',
                      geom_opt_tol_gradient=300,
                      geom_opt_tol_energy=100,
                      geom_opt_tol_displacement=1200)

parsed_data, electronic_structure = get_output_from_qchem(qc_input,
                                                          processors=4,
                                                          parser=basic_optimization,
                                                          read_fchk=True)

qc_input = QchemInput(parsed_data['optimized_molecule'],
                      jobtype='freq',
                      exchange='hf',
                      basis='sto-3g',
                      scf_guess=electronic_structure['coefficients'])

parsed_data = get_output_from_qchem(qc_input,
                                    processors=4,
                                    parser=basic_frequencies)


for i, mode in enumerate(parsed_data['modes']):
    print('mode:                      {}'.format(i+1))
    print('frequency (cm-1):          {:10.2f}'.format(mode['frequency']))
    print('force constant (mdyne/A):  {:10.5f}\n'.format(mode['force_constant']))

```
**Custom basis support**

```python
from pyqchem import QchemInput, Structure
from pyqchem.basis import get_basis_from_BSE


molecule = Structure(coordinates=[[0.0, 0.0, 0.0000],
                                  [0.0, 0.0, 1.5811]],
                     symbols=['Se', 'H'],
                     charge=-1,
                     multiplicity=1)

basis_custom = get_basis_from_BSE(molecule, 'cc-pVTZ')

qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      basis=basis_custom)

```

**Handle qchem errors**


```python
from pyqchem import get_output_from_qchem
from pyqchem.errors import OutputError, ParserError
from pyqchem.parsers.parser_rasci import parser_rasci

try:
    parsed_data = get_output_from_qchem(qc_input,
                                        processors=4,
                                        parser=parser_rasci,
                                        )

except OutputError as e:
    print('Calculation ended with errors. Error lines:')
    print(e.error_lines)
    
    # Try to parse your output anyway
    try: 
        parsed_data = parser_rasci(e.full_output)
    except ParserError:
        print('Failed parsing')
        exit()


print('Energy: ', parsed_data['scf_energy'])
```
Additional example scripts are found at the [examples folder](/examples)

Contact info
------------
Abel Carreras  
abelcarreras83@gmail.com

Donostia International Physics Center (DIPC)  
Donostia-San Sebastian (Spain)
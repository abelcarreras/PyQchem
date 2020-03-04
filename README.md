[![Build Status](https://travis-ci.com/abelcarreras/PyQchem.svg?branch=master)](https://travis-ci.com/abelcarreras/PyQchem)
[![Coverage Status](https://coveralls.io/repos/github/abelcarreras/PyQchem/badge.svg?branch=development)](https://coveralls.io/github/abelcarreras/PyQchem?branch=development)
[![PyPI version](https://badge.fury.io/py/pyqchem.svg)](https://badge.fury.io/py/pyqchem)

PyQchem
=======
Python wrapper for Q-Chem (https://www.q-chem.com)

Main features
-------------
- Easy to use clean python interface for Q-Chem
- No special q-chem compilation needed
- Output parser support 

Requirements
------------
- Python 2.7.x/3.5+
- numpy
- scipy
- matplolib
- wfnsympy
- requests
- lxml
------------
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
                     atomic_elements=['H', 'H'],
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
print('Energy: ', data['scf energy'])
```

**Link calculations in powerful workflows**

```python
from pyqchem import QchemInput, get_output_from_qchem
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.parsers.parser_frequencies import basic_frequencies


parsed_data, electronic_structure = get_output_from_qchem(qc_input,
                                                          processors=4,
                                                          parser=basic_optimization,
                                                          read_fchk=True)

opt_molecule = parsed_data['optimized_molecule']


qc_input = QchemInput(opt_molecule,
                      jobtype='freq',
                      exchange='hf',
                      basis='sto-3g',
                      scf_guess=electronic_structure['coefficients'])

parsed_data = get_output_from_qchem(qc_input,
                                    processors=4,
                                    parser=basic_frequencies)


for mode, freq in enumerate(parsed_data['frequencies']):

    force_constants = parsed_data['force_constants'][mode]

    print('mode:                      {}'.format(mode+1))
    print('frequency (cm-1):          {:10.2f}'.format(freq))
    print('force constant (mdyne/A):  {:10.5f}\n'.format(force_constants))

```

**Handle qchem errors like a pro!**

```python
from pyqchem import get_output_from_qchem
from pyqchem.errors import OutputError, ParserError
from pyqchem.parsers.parser_rasci import rasci as parser_rasci

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


print('Energy: ', parsed_data['scf energy'])
```

Contact info
------------
Abel Carreras  
abelcarreras83@gmail.com

Donostia International Physics Center (DIPC)  
Donostia-San Sebastian (Spain)
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
```
# python setup.py install --user
```


Interface 
---------
**Simple pythonic API to define your input**

```
molecule = Structure(coordinates=[[0.0, 0.0, 0.0],
                                  [0.0, 0.0, 0.9]],
                     atomic_elements=['H', 'H'],
                     charge=0,
                     multiplicity=1)

qc_input = QchemInput(molecule,
                      jobtype='sp',
                      exchange='hf',
                      basis='6-31G(d,p)')

data = get_output_from_qchem(qc_input,
                             processors=4,
                             parser=basic_parser_qchem)

# obtain a python dictionary
print('Energy: ', data['scf energy'])
```

**Link calculations in powerful workflows**

```
parsed_data, electronic_structure = get_output_from_qchem(qc_input,
                                                          processors=4,
                                                          parser=basic_optimization,
                                                          read_fchk=True)

opt_molecule = parsed_data['optimized_molecule']


qc_input = create_qchem_input(opt_molecule,
                              jobtype='freq',
                              exchange='hf',
                              basis='sto-3g',
                              scf_guess=electronic_structure['coefficients'])

parsed_data = get_output_from_qchem(qc_input,
                                    processors=4,
                                    parser=basic_frequencies)

```

**Handle qchem errors like a pro!**

```
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
        parsed_data = parser_rasci(e.output)
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
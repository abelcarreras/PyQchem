Extra information
=================

Electronic structure
--------------------
The electronic structure dictionary is designed to contain the data
from sources other than the output. Most of its contents are from the
fchk file, but also contains data from scratch files such as the hessian
and the fock matrix. Other data will be included in this dictionary in
the future.

The basic structure of the electronic structure dictionary is the following:

..  code-block:: bash

    root
     ├── basis
     │   ├── name
     │   ├── primitive_type
     │   └── atoms(list)
     │         ├── shells (list)
     │         ├── symbol
     │         └── atomic_number
     ├── coefficients
     │   ├── alpha
     │   └── beta (optional)
     ├── mo_energies
     │   ├── alpha
     │   └── beta (optional)
     ├── number_of_electrons
     │   ├── alpha
     │   └── beta
     ├── nato_coefficients (optional)
     │   ├── alpha
     │   └── beta
     ├── nato_occupancies (optional)
     │   ├── alpha
     │   └── beta
     ├── structure
     └── overlap


Using the information of this dictionary a Fchk file can be generated. This may be used to visualize the molecular
orbitals, electronic density and other properties using an (external) visualization program.

..  code-block:: python

    from pyqchem.file_io import build_fchk
    with open('file.fchk', 'w') as f:
        f.write(build_fchk(electronic_structure))


While electronic structure is a simple dictionary, its elements are designed to be interoperable along the
pyqchem code such as guess and basis. Some examples of this interoperability can be found in the examples folder.

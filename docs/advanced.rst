.. highlight:: rst

Advanced input
==============

The structure of the electronic structure dictionary is as follows:

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
orbitals using an external program.

..  code-block:: python

    from pyqchem.file_io import build_fchk
    with open('file.fchk', 'w') as f:
        f.write(build_fchk(electronic_structure))


While electronic structure is a simple dictionary, its elements are designed to be interoperable along the
pyqchem code. Some examples of this interoperability will be shown in the following sections.

Custom guess
------------

*structure* entry contains a Structure type object that can be used, for instance, in *QchemInput* to create a new input.
*scf_guess* requires a dictionary {'alpha': [], 'beta':[]} containing the molecular orbitals coefficients matrix. This
is exactly what *coefficients* entry contains. Using this two objects it becomes easy to use the electronic structure
of a previous calculation as a guess of a new calculation (see **frequencies_simple.py** example) :

..  code-block:: python

    qc_input = QchemInput(electronic_structure['structure],
                          scf_guess=electronic_structure['coefficients']
                          jobtype='sp',
                          exchange='hf',
                          basis='6-31G')


Custom basis set
----------------

The same can be done for basis set. *QchemInput* basis argument accepts predefined basis sets included in Q-Chem as
a label string (e.g. 'sto-3g', '6-31g(d,p)',..) but also accepts custom basis sets. These basis sets should be written
as a python dictionary following the same structure as the one output in *electronic_structure*. These basis can be used
directly in *QchemInput*:

..  code-block:: python

    qc_input = QchemInput(electronic_structure['structure],
                          scf_guess=electronic_structure['coefficients']
                          jobtype='sp',
                          exchange='hf',
                          basis=electronic_structure['basis])


However this is may not very useful if the basis in *electronic_structure* is one of the predefined basis in Q-Chem.
PyQchem include a helper function to retrieve a basis set from *ccRepo* (http://www.grant-hill.group.shef.ac.uk/ccrepo/)
repository. This function require as argument Structure object and the name of the basis set (see: **custom_basis.py** example):

..  code-block:: python

    from pyqchem.basis import get_basis_from_ccRepo

    basis_custom_repo = get_basis_from_ccRepo(molecule, 'cc-pVTZ')
    qc_input = QchemInput(molecule,
                          jobtype='sp',
                          exchange='hf',
                          basis=basis_custom_repo)


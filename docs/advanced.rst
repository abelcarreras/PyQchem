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


Dual basis set
--------------
The use of dual basis set can improve the performance of Q-Chem calculations. This can be used, for example, to use as
a guess a previous calculations tha uses a smaller basis set. The keyword to use this is *basis2* and works in the same way
as *basis*. Usual *basis* keyword defines the new basis and *basis2* keyword defines the previous (and smaller) basis.

..  code-block:: python

    # Initial calculation using sto-3g basis set
    qc_input = QchemInput(molecule,
                          jobtype='sp',
                          exchange='hf',
                          basis='sto-3g',
                          )

    _, ee = get_output_from_qchem(qc_input, return_electronic_structure=True)

    # Precise calculation with larger 6-31G basis using previous MO as guess
    qc_input = QchemInput(molecule,
                          jobtype='sp',
                          exchange='hf',
                          basis='6-31g',
                          basis2=ee['basis'],  # previous basis can also be read from electronic structure dictionary
                          scf_guess=ee['coefficients'] # previous MO coefficients to be used as a guess
                          )

Usage of Solvent
----------------
Usage of solvent is implemented in pyQchem by the use of *solvent_method* and *solvent_params*. *solvent_method*
is a strightforward of the keyword with the same name in Q-Chem while *solvent_params* is a dictionary that
contains the keywords in the section **$solvent** in Q-Chem input. For PCM that requiere additional
parameters *pcm_params* keyword is used which implements the keywords of **$pcm** section in Q-Chem input.

..  code-block:: python

    qc_input = create_qchem_input(molecule,
                                  jobtype='sp',
                                  exchange='hf',
                                  basis='sto-3g',
                                  unrestricted=True,
                                  solvent_method='pcm',
                                  solvent_params={'Dielectric': 8.93},  # Cl2CH2
                                  pcm_params={'Theory': 'CPCM',
                                              'Method': 'SWIG',
                                              'Solver': 'Inversion',
                                              'Radii': 'Bondi'}
                                  )

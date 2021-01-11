.. highlight:: rst

Tutorial
========

PyQchem is python interface for Q-Chem. It allows to create Q-Chem inputs, execute Q-Chem from python, parse its
outputs and store the results in convenient python dictionaries. This is especially useful to create complex workflows
to automate frequent tasks using python programming language.

As a preparation for the incoming talk I prepared a series of exercises to introduce the very basics of PyQchem. These
exercises are preceded by a detailed explanation of the PyQchem functionality that you may need to complete them. For further
information you can check the PyQchem manual available online at: https://pyqchem.readthedocs.io/

Once you complete these exercises, you will be able to:

- Submit a simple Q-Chem calculation from python and obtain the desired results
- Use a simple loops to automate the creation of q-Chem inputs.
- Combine two Q-Chem calculation by using the outputs of the former in the input of the later.


Execution
---------
In order to run Q-Chem calculations using PyQchem a installation of Q-Chem is necessary. Q-Chem is currently installed
in ATLAS cluster along with PyQchem, hence if you use ATLAS it is not necessary any further installation to complete these exercises.
PyQchem will be loaded along Q-Chem when loading qchem_group/qchem_trunk modules as usual ::

    export MODULEPATH=/scratch/user/SOFTWARE/privatemodules:$MODULEPATH
    export QCSCRATCH=/scratch/user/QCHEM_SCRATCH
    module load qchem_group

    python script.py > output.txt

..  note::
    Since the exercises contained in this document are quite short, you can run them directly in ATLAS
    (so you don't have to wait for the queue system) just connect to ATLAS and export MODULEPATH (as shown above)
    and load qchem_group module.

However, it may be useful to have PyQchem installed locally in your computer to prepare scripts, specially when using some
advanced python editors like PyCharm, Clion, VCode, etc.
PyQchem can be downloaded and installed in a MAC/Linux machine from the official python repository (https://pypi.org)
by using the command ::

    pip install pyqchem

In some python installations this command may need sudo permissions. If this is the case, then you can specify user's home
installation path by ::

    pip install pyqchem --user

this way the installation will be done in your home and will not require superuser permissions.

Basic concepts
--------------
In order to perform a Q-Chem calculation, two main pieces of information are needed as input:

1) the molecular structure
2) the parameters of the calculation (method/basis set/etc..)


To define the molecular structure, PyQchem uses a *Structure* class that can be imported as:

..  code-block:: python

    from pyqchem import Structure

Using this class you can create an instance of this class. An instance can be understood as a variable with type *Structure*.
To create an instance it is necessary to initialize it with the necessary parameters. In the case of *Structure* class
these parameters are *coordinates*, *symbols*, *charge* and *multiplicity*:

..  code-block:: python

    from pyqchem import Structure

    hydroxide = Structure(coordinates=[[0.0, 0.0, 0.0],
                                       [0.0, 0.0, 0.9]],

                          symbols=['O', 'H'],
                          charge=-1,
                          multiplicity=1)

As can be seen in the example above, both *coordinates* and *symbols* are simple python lists separated by comma and
limited by **[** and **]** characters. In the case of *symbols* list, its elements are character strings which are surrounded
by quotes. *charge* and *multiplicity* are simple integer numbers and are optional parameters. In case of not being defined,
charge will be set to 0 and multiplicity to 1.

This instance has all the methods defined for the *Structure* class. Methods can be understood as functions associated
to a particular class instance. To access to these methods the operator '.' is used followed by the method's name.
The following example shows some of the methods defined in *Structure* class:

..  code-block:: python

    # variables
    ne = hydroxide.number_of_electrons
    alpha = hydroxide.alpha_electrons
    beta = hydroxide.beta_electrons

    print('Number of electrons:', ne, '(', alpha, beta, ')')

    # functions
    xyz_file_txt = hydroxide.get_xyz(title='hydroxide anion')
    print(xyz_file_txt)

..  note::
    In this example xyz_file_txt is a string that is printed on screen using *print()* function. You can write
    this string into a file using python language, for example: ::

        open('file_name.xyz').write(xyz_file_txt)



The definition of the parameters is done by the *QchemInput* class.
Using this class we define an instance of this class as:

..  code-block:: python

    from pyqchem import QchemInput

    oh_input = QchemInput(molecule,
                          jobtype='sp',
                          exchange='hf',
                          basis='6-31G',
                          unrestricted=True)

In a similar way as the *Structure* class, to initialize a *QchemInput* instance several parameters are necessary.
In this case the first parameter is *molecule*. *molecule* is an instance of the *Structure* class, like the one that we
defined before (hydroxide). All the other parameters are optional and have default parameters in case of not being defined.
The name of these parameters is designed to be equal or similar to the respective Q-Chem keywords. The list of available
parameters is being updated continuously and can be found in: https://github.com/abelcarreras/PyQchem/blob/master/pyqchem/qc_input.py

As in the case of *Structure* class, several methods are defined for *QchemInput*. The main one is **get_txt()**.
This method returns a string containing the input in Q-Chem format. This can be used to check the exact input that
will be submitted to Q-Chem to do the calculation.

..  code-block:: python

    input_txt = oh_input.get_txt()
    print(input_txt)

Other useful methods are **get_copy()** and **update_input()**. These methods are useful to modify already created inputs. For example,
in case you want to prepare multiple different inputs with few differences you can create a general input, make
multiple copies of it and modify them:

..  code-block:: python

    input_txt = oh_input.get_txt()
    print(input_txt)


    general_input = QchemInput(molecule,
                               jobtype='sp',
                               exchange='hf')

    input_basis_1 = general_input.get_copy()
    input_basis_2 = general_input.get_copy()

    input_basis_1.update_input({'basis': 'sto-3g', 'mem_total' : 2000})
    input_basis_2.update_input({'basis': '6-31G', 'mem_total' : 1000})


Finally, to run the calculation **get_output_from_qchem** function is used. The first argument of this function is
a *QchemInput* instance. There are several optional parameters for this function mainly related to computer stuff
(which do not affect the results of the calculation). A good representative is *processors*, that indicate the number of
processor cores to use in the calculation (in openMP compilation) or the number of MPI processes (in MPI compilation).

..  note::
    In ATLAS cluster Q-Chem is compiled using openMP.

..  code-block:: python

    from pyqchem import get_output_from_qchem

    output = get_output_from_qchem(oh_input,
                                   processors=4)

    print(output)


The output of this function is a string containing the full Q-Chem output. In this example the output
is printed in the screen. Combining all these functions together we obtain a simple script that runs
a single Q-Chem calculation and prints its output:


..  code-block:: python

    from pyqchem import Structure, QchemInput, get_output_from_qchem

    hydroxide = Structure(coordinates=[[0.0, 0.0, 0.0],
                                       [0.0, 0.0, 0.9]],

                          symbols=['O', 'H'],
                          charge=-1,
                          multiplicity=1)

    oh_input = QchemInput(hydroxide,
                          jobtype='sp',
                          exchange='hf',
                          basis='6-31G',
                          unrestricted=True)


    output = get_output_from_qchem(oh_input,
                                   processors=4)

    print(output)



Practical exercises
"""""""""""""""""""

a) Use PyQchem to write a script that generates a set of Q-Chem inputs to do a HF calculation of the methane molecule
using the following basis sets: STO-3G, 6-31G, DZ, cc-pVDZ and aug-cc-pVDZ. Use python's **open()** function to store
these inputs in different files.

b) Modify the previous script to run the generated inputs using **get_output_from_qchem()** function to obtain the corresponding
Q-Chem outputs. Store these outputs in different files.

c) [ADVANCED] Make use of python language tools such as *list comprehension* and *for/while* loops to make this
exercise, obtaining a cleaner and more extendable code.

Parsing data
------------
Being able to automatically generate Q-Chem inputs and outputs can be pretty useful. However the key feature of
PyQchem is the use of parsers to extract the output information. A parser is just a function that takes a text string,
finds the important data and places it in an organized structure. In the case of PyQchem this structure is a python dictionary.

Here an example of such a function:

..  code-block:: python

    def parser_example(output):
        data_dict = {}
        enum = output.find('Total energy in the final basis set')
        data_dict['scf_energy'] = float(output[enum: enum+100].split()[8])
        return data_dict

    print(output)

This function does 3 main things:

* Define a python dictionary using {} syntax.
* Get the location of the data that we are interested in the output, in this case the SCF energy.
* Convert the interesting data from text format to number format using *float()* function and store them in the dictionary.

..  note::
    The use of *[ini: fin]* in strings to get a substring is called slicing. This is very useful in parser functions
    since you can divide a long output string in small strings that contain the data. On the other hand *split()*
    method divides a text in words and generates a list that can be accessed by indices. ::

        text = 'this may be a long text with lots of words'
        subtext = text[0: 11]   # Result: 'this may be'
        words = subtext.split() # Result: ['this', 'may', 'be']
        word = words[1]         # Result: 'may'

..  note::
    In contrast to other languages like Fortran Python indices start from 0 (not 1!).

Parser functions can be explicitly written in the python script just after getting the Q-Chem output:

..  code-block:: python

    def parser_example(output):
        data_dict = {}
        enum = output.find('Total energy in the final basis set')
        data_dict['scf_energy'] = float(output[enum: enum+100].split()[8])
        return data_dict

    (...)

    output = get_output_from_qchem(oh_input)

    parsed_data = parser_example(output)
    print(parsed_data)  # Result: {'scf_energy': 1.234567}

The above example will print a dictionary with a single item with key **'scf_energy'** and the energy as a value. A python
dictionary works in a similar way as list/vectors, but instead of accessing the elements with an integer index we use
a key string.

..  code-block:: python

    energy = parser_data['scf_energy']
    print ('The energy is ', energy)

As may be expected, a dictionary can contain multiple items so accessing them via keys is a basic functionality.
The values of a dictionary can be almost anything: strings, numbers, lists ... and even other dictionaries. This
generates a very common structure of dictionaries inside dictionaries used to organize the data in a tree-like structure.

..  code-block:: python

    sub_dict = {}
    dict = {}

    sub_dict['inside'] = [4, 3, 5]
    dict['outside'] = sub_dict

    print(dict['outside']['inside'])  # Result: [4, 3, 5]
    print(dict['outside']['inside'][2])  # Result: 5

..  note::
    Technically a dictionary key can be other objects aside from strings but to make it simple we will use strings.

The use of parsers in PyQchem is kind of a basic feature, for this reason **get_output_from_qchem()** function has an
optional argument that requires a parser function:

..  code-block:: python

    def parser_example(output):
        data_dict = {}
        enum = output.find('Total energy in the final basis set')
        data_dict['scf_energy'] = float(output[enum: enum+100].split()[8])
        return data_dict

    (...)

    parsed_data = get_output_from_qchem(oh_input, parser=parser_example)

    print(parsed_data)  # Result: {'scf_energy': 1.234567}


As can be observed in the example above, using parser argument transforms the output of **get_output_from_qchem** into
a dictionary with the parsed output. This makes the script shorter and cleaner.
PyQchem package includes parsers written for the most common types of calculations. These can be found in the parsers
folder: (https://github.com/abelcarreras/PyQchem/tree/master/pyqchem/parsers). To use them, you just need to use *import*
statement:

..  code-block:: python

    from pyqchem.parsers.basic import basic_parser_qchem

    (...)

    parsed_data = get_output_from_qchem(oh_input,
                                        parser=basic_parser_qchem)


Practical exercises
"""""""""""""""""""

a) Create a parser function to get the following properties from a HF calculation: *Sum of atomic charges* & *Sum of spin   charges*.
You can use the same system as in the fist example (methane with STO-3G basis set) to test it.

b) Use the basic parser included in PyQchem (*basic_parser_qchem*) to write a script
that calculates and the orbital energies of methane molecule.

c) [ADVANCED] Use a loop (*for/while*) to calculate the scf energy of the hydrogen molecule at different geometries (bond length)
to study the dissociation of hydrogen molecule. Print the results as two columns (bond length and scf energy)

..  note::
    During the execution a *calculation_data.pkl* file is generated. This stores data of previous calculations
    (see manual for more information). Modifying the parser may make this data obsolete, if something unexpected happens
    modifying the parser try removing this file. Also, see *force_recalculation=True* argument of **get_output_from_qchem()** function.


Linking calculations
--------------------
One of the strongest reasons to use a library like PyQchem is the ability to link different calculations together.
This means prepare inputs from output data of previous calculations. A typical example is the calculation of the normal
modes frequencies of a previously optimized structure. This can be done in PyQchem in the following way:

..  code-block:: python

    from pyqchem.parsers.parser_frequencies import basic_frequencies
    from pyqchem.parsers.parser_optimization import basic_optimization

    (...)

    opt_input = create_qchem_input(molecule,
                                   jobtype='opt',
                                   exchange='hf',
                                   basis='sto-3g')

    parsed_opt_data = get_output_from_qchem(opt_input, parser=basic_optimization)

    opt_molecule = parsed_opt_data['optimized_molecule']

    freq_input = create_qchem_input(opt_molecule,
                                    jobtype='freq',
                                    exchange='hf',
                                    basis='sto-3g')

    parsed_data = get_output_from_qchem(freq_input, parser=basic_frequencies)

    print(parsed_data)


In this example, the optimized structure is obtained from the parsed output of the optimization calculation.
In this parser the value of **'optimize_molecule'** entry is already an instance of the *Structure* class so it can be
used directly in the frequencies calculation input.

This script is pretty convenient but it can be done even better. In order to take maximum profit of a previous
calculation, the already optimized electronic structure (molecular orbitals) can be used as a initial guess in the frequencies
calculation. To do this, it is necessary to get the orbitals coefficients, which are not present in the usual output.
PyQchem obtains the electronic structure data from the *FChk* file. The request of the *FChk* generation is done directly
in the **get_output_from_qchem** function by using the argument *read_fchk=True*. This modifies the output of this function
returning two pieces of data (a list of two elements): the parsed output & the parsed *FChK* data:

..  code-block:: python

    from pyqchem.parsers.parser_frequencies import basic_frequencies
    from pyqchem.parsers.parser_optimization import basic_optimization

    (...)

    parsed_opt_data, electronic_structure = get_output_from_qchem(opt_input,
                                                                  parser=basic_optimization,
                                                                  read_fchk=True)

    print(electronic_structure)

if you print *electronic_structure* you will notice that it is a dictionary containing the entries of a usual FChk file.
Due to the standard format of this file all data is parsed so it is not necessary to indicate a parser. The format
of this dictionary is designed for inter-operation with the different functions of PyQchem. A simple example is the use
of the molecular orbitals coefficients as an initial guess:

..  code-block:: python

    mo_coefficients = electronic_structure['coefficients']

    freq_input = create_qchem_input(opt_molecule,
                                    jobtype='freq',
                                    exchange='hf',
                                    basis='sto-3g',
                                    scf_guess=mo_coefficients)

In this case, **electronic_structure['coefficients']** contains a NxN square matrix with the coefficients of the
molecular orbitals where each row corresponds to a molecular orbital.

Practical exercises
"""""""""""""""""""

a) Use PyQchem to optimize the SO2 molecule using HF and minimum basis set (STO-3G). From
the optimized structure perform 3 additional optimizations using larger 3 different basis sets: SV, DZ & TZ.
Get the scf_energies from each optimization and store the optimized structures in XYZ files.

b) (ADVANCED) Perform a frequencies calculation of the SO2 molecule using PyQchem (with HF/STO-3G) and create
a movie in a XYZ file that shows the vibration of each normal mode. Print the results of the **basic_frequencies** parser
and investigate its contents to find the necessary information (*displacements*).
(https://github.com/abelcarreras/PyQchem/blob/master/pyqchem/parsers/parser_frequencies.py)


..  note::
    To create a movie in XYZ just put all the geometries (one under the other) in the same file. This will be interpreted
    in most molecular visualization software (Ex. VMD) as frames and you will be able to reproduce them as a movie.

    HINT: In python you can combine two strings by the + operator. Ex: ::

        video_xyz = frame1_xyz + frame2_xyz + frame3_xyz


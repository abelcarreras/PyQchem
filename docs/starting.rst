.. highlight:: rst

Get started
===========

A basic pyqhem script is composed of 4 steps: Defining a molecule, preparing an Q-Chem input, runing
the calculation and parsing the information.

Defining molecule
-----------------

The definition of the molecule is done creating an instance of Structure class. Its initialization requires
the coordinates as a list of lists (or Nx3 numpy array), the atomic symbols of the atoms (in the same order
of the coordinates), the charge and the multiplicity. Only coordinates and symbols are mandatory, if
charge/multiplicity are not specified they will be defined as neutral/singlet.


Example of hydroxide anion:

.. code-block:: python

    from pyqchem import Structure

    molecule = Structure(coordinates=[[0.0, 0.0, 0.0],
                                      [0.0, 0.0, 0.9]],
                         symbols=['O', 'H'],
                         charge=-1,
                         multiplicity=1)


Preparing Q-Chem input
----------------------

The Q-Chem input is defined using the QchemInput class. To initialize this class it is necessary an instance of
the Structure class as defined in the previous section (molecule). Afterwards the different Q-Chem keywords are
set as optional arguments to specify the calculation. In general the name of these keywords have the same name as
in Q-Chem. At this moment, only a small set of the keywords available in Q-Chem are implemented in the QchemInput
class. Refer at this class definition to check which ones are implemented. New keywords will be implemented under
request.


Example of single point calculation using unrestricted Hartree-Fock method with 6-31G basis set:

.. code-block:: python

    from pyqchem import QchemInput

    qc_input = QchemInput(molecule,
                          jobtype='sp',
                          exchange='hf',
                          basis='6-31G',
                          unrestricted=True)

.. note::

    In case a particular Q-Chem keyword of *$REM* section in not implemented in PyQchem, *extra_rem_keywords* argument
    can be used to include it. This argument requires a dictionary containing the keywords as dictionary keys and its
    values as dictionary values. Values can be either strings or numbers.

    .. code-block:: python

            qc_input = QchemInput(molecule,
                                  jobtype='sp',
                                  exchange='hf',
                                  basis='6-31G',
                                  extra_rem_keywords={'keyword_1': 'value',
                                                      'keyword_2': 34}


    Also if a non implemented Q-Chem section is required, *extra_sections* argument
    can be used to include it. This argument requires a list of CustomSection objects
    which can be imported from pyqchem.qc_input.


    .. code-block:: python
            from pyqchem.qc_input import CustomSection
            qc_input = QchemInput(molecule,
                                  jobtype='sp',
                                  exchange='hf',
                                  basis='6-31G',
                                  extra_sections=[CustomSection(title='section_title',
                                                                keywords={'sec_keyword_1': 'value',
                                                                          'sec_keyword_2': 34})]


Running calculations
--------------------
Once the QchemInput is prepared you can run it using *get_output_from_qchem* function. This function if the core
of PyChem and makes the connection between PyQchem and Q-Chem. In order for this function to work properly
**$QC** and **$QCSCRATCH** environment variables should be defined. Refer to Q-Chem manual to learn how to set them.
By default *get_output_from_qchem* will assume an openMP compilation of Q-Chem were *processors* indicate the number
of threads to use in the calculation. If your compilation of Q-Chem is MPI then *use_mpi=True* should be set and
*processors* will correspond to the number of processors.

*Get_output_from_qchem* function has two main ways of operation. If no parser is specified, then the output of this
function will be a string containing the full Q-Chem output. This way can be useful to do your own treatment of the
output file or if you are not sure about the information you want to parse.


Example of simple parallel(openMP) calculation using 4 threads:

.. code-block:: python

    from pyqchem import get_output_from_qchem

    output = get_output_from_qchem(qc_input,
                                   processors=4)


The second way is by defining the *parser* optional argument. This indicates that the output will be parsed
using the specified parser function. In the following example *basic_parser_qchem* function is used. This is
imported from the parser collection located at *pyqchem.parsers.*. Using a parser function the output of this
function becomes a python dictionary containing the parsed data.


This is similar to the example shown above using a simple parser (*basic_parser_qchem*) :

..  code-block:: python

    from pyqchem.parsers.basic import basic_parser_qchem
    parsed_data = get_output_from_qchem(qc_input,
                                        processors=4,
                                        parser=basic_parser_qchem,
                                        )


This can be done also in two steps, since the parser (*basic_parser_qchem* in this case) is a just regular python
function that accepts a string as argument.

..  code-block:: python

    output = get_output_from_qchem(qc_input, processors=4)
    parsed_data = basic_parser_qchem(output)


It is simple to create a custom parser by defining a custom function with the following structure:

..  code-block:: python

    def custom_parser_qchem(output):
        """
        output: contains the full Q-Chem output in a string

        return: a dictionary with the parsed data
        """
        ...
        return {'property_1': prop1,
                'property_2': prop2}


Complex parsers may have optional arguments to add more control. This may be used to include parameters such as
precision, max number of cycles/states/etc to read, etc..:

..  code-block:: python

    def custom_parser_qchem(output, custom_option=True, custom_prec=1e-4):
        """
        output: contains the full Q-Chem output in a string
        custom_option: controls option to be used or not
        custom_prec: defines the precision of som data to be read

        return: a dictionary with the parsed data
        """
        ...


        return {'property_1': prop1,
                'property_2': prop2}

to define this optional arguments *get_output_from_qchem* function you should include *parser_parameters* argument
which requires a python dictionary. Each of the entries in this dictionary should be the name of one of the optional
arguments in the parser function whose value is the value of the argument:

..  code-block:: python

    parsed_data = get_output_from_qchem(qc_input,
                                        processors=4,
                                        parser=custom_parser_qchem,
                                        parser_parameters={'custom_option': True, 'custom_prec': 1e-4}
                                        )


Electronic structure
--------------------
Most of the electronic information (molecular orbitals coefficients, electronic density, basis set, etc..) can be found
in fchk file generated by Q-Chem. In usual Q-Chem calculation to generate the *fchk* file it is necessary to include
the keyword *gui=2* in the input file. Using pyqchem this is not necessary, you can request to generate this file
and parser its contents using the argument *read_fchk=True*:

..  code-block:: python

    from pyqchem.parsers.basic import basic_parser_qchem
    parsed_data, electronic_structure = get_output_from_qchem(qc_input,
                                                              processors=4,
                                                              parser=basic_parser_qchem,
                                                              read_fchk=True
                                                              )


as can be observed in the previous example, the return of *get_output_from_qchem* function contains two elements:
*parsed_data* and the *electronic_structure*. *Parsed_data* is a python dictionary that contains the same information
as previously described. *Electronic_structure* is another python dictionary that contains the information parsed from
the FCHK file.

Reusing data efficiently
------------------------
Pyqchem is specially focused in the automation and design of complex Q-Chem workflows. For this reason pyqchem
implements a feature to avoid redundant calculation by storing the parsed data in a pickle file. This works
seamessly, if a calculation is requested with an input *equivalent* to a previous one, the calculation is skip
and stored data is output instead. By default only parsed data is stored, therefore if no parser is provided
the calculation will be recomputed.

The behavior of this feature is controlled by two arguments in *get_output_from_qchem* function:
*force_recalculation* and *store_full_output*. *force_recalculation=True* forces the calculation to be calculated
even if a previous *equivalent* calculation already exists.
If *store_full_output=True* then the raw outputs are also stored. This may produce a significant
increase in size of the storage file, but it can be useful to test new parsers or to use several parsers in
the same output.

..  code-block:: python

    parsed_data = get_output_from_qchem(qc_input,
                                        processors=4,
                                        parser=basic_parser_qchem,
                                        force_recalculation=True,
                                        store_full_output=True
                                        )


It is possible to set a custom storage pickle filename by using *redefine_calculation_data_filename* function.
This may be written at the beginning of the script to define a different storage file for each script if
multiple scripts run in the same directory at the same time.

..  code-block:: python

    from pyqchem.qchem_core import redefine_calculation_data_filename
    redefine_calculation_data_filename('custom_file.pkl')


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


Running calculations
--------------------
Once the QchemInput is prepared you can run it using *get_output_from_qchem* function. This function if the core
of PyChem and makes the connection between PyQchem and Q-Chem. In order for this function to work properly
**$QC** and **$QCSCRATCH** environment variables should be defined. Refer to Q-Chem manual to learn how to set them.
By default *get_output_from_qchem* will assume an openMP compilation of Q-Chem were *processors* indicate the number
of threads to use in the calculation. If your compilation of Q-Chem is MPI then *use_mpi=True* should be set and
*processors* will correspond to the number of processors.

*get_output_from_qchem* function has two main ways of operation. If no parser is specified, then the output of this
function will be a string containing the full Q-Chem output. This way can be useful to do your own treatment of the
output data or if you are not sure about the information you want to parse.


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
    data = get_output_from_qchem(qc_input,
                                 processors=4,
                                 parser=basic_parser_qchem,
                                 )

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


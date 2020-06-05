.. highlight:: rst

Error treatment
===============

In the previous sections it is always assumed that Q-Chem calculations always finish successfully,
but, in general, this is not true. Calculation may fail for several reasons, incorrect input parameters,
slow convergence, problems with the cluster, issues during the parser process, etc..
To handle these errors PyQchem implements two custom error classes: **OutputError** and **ParserError**.

These errors are risen during the *get_output_from_qchem* function execution if something has gone wrong.
If **OutputError** is risen means that Q-Chem calculation did not finish correctly. This may mean that the
input is incorrect or something happened with the computer that lead the calculation to crash. When this
error is risen the traceback will show the last 20 lines of the Q-Chem output where, hopefully, you will find
enough information to determine the cause of the error.

..  code-block:: bash

    Traceback (most recent call last):
      File "/Users/abel/PycharmProjects/qchem_scripts/scripts/test_tddft.py", line 51, in <module>
        store_full_output=True,
      File "/Users/abel/PycharmProjects/qchem_scripts/pyqchem/qchem_core.py", line 330, in get_output_from_qchem
        raise OutputError(output, err)
    pyqchem.errors.OutputError: Error in Q-Chem calculation:
    sh: gdb: command not found
    Unable to call dbx in QTraceback: No such file or directory
     http://arma.sourceforge.net/

     Q-Chem begins on Fri Jun  5 02:31:23 2020

    Host:
    0

         Scratch files written to /Users/abel/Scratch/export/qchem66428//

     Name = B3LP

     Q-Chem fatal error occurred in module libdft/dftcodes.C, line 599:

     Unrecognized exchange functional in XCode


     Please submit a crash report at q-chem.com/reporter


In some cases, specially when working with complex workflows it becomes useful to be able to capture these errors
to handle them automatically without finish all the workflow. This can be done by the usual try/except statements.
By capturing the error it is possible to skip the calculation, recover the error_lines and even recover the full Q-Chem
output:

..  code-block:: python

    try:
        parsed_data = get_output_from_qchem(qc_input,
                                            processors=4,
                                            parser=basic_optimization)
    except OutputError as e:
        print('These are the error lines:\n', e.error_lines)
        print('This is the full output:\n', e.full_output)


Using the full output it is possible to try to parse the information that contains by applying the parser directly:


..  code-block:: python

    try:
        parsed_data = get_output_from_qchem(qc_input,
                                            processors=4,
                                            parser=basic_optimization)
    except OutputError as e:
        print('Something wrong happened!:\n', e.error_lines)
        print('Recovering usable data...')
        parsed_data = basic_optimization(e.full_output)


But if the calculation is really incomplete, or the format of Q-Chem output is incompatible with the parser,
the parsing process may also fail and a **ParserError** will be risen. In this case the output data cannot
be recovered using this parser.

In the same way as **OutputError** a try/except block can be written to capture this error. A sensible to proceed
can be either skip the calculation or try another parser by nesting two try/except blocks :

..  code-block:: python

    try:
        parsed_data = get_output_from_qchem(qc_input,
                                            processors=4,
                                            parser=basic_optimization)
    except OutputError as e:
        print('Something wrong happened!:\n', e.error_lines)
        print('Recovering usable data...')

        try:
            parsed_data = basic_optimization(e.full_output)
        except ParserError:
            print('Trying another parser')
            parsed_data = other_parser(e.full_output)



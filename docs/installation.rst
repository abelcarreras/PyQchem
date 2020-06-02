.. highlight:: rst

Installation
============

PyQchem can be installed directly from the source code (python package) or via PyPI repository.
Q-Chem is not necessary to be installed in the system but, of course, it will be necessary later
on to perform calculations. Still some of pyqchem features can be used without Q-Chem.

Requirements
------------

- Python 2.7.x/3.5+
- numpy
- scipy
- matplolib
- requests
- lxml
- wfnsympy (optional: for symmetry analysis)
- paramiko (optional: for remote calculations)

Install
-------

1) From source code ::

    python setup.py install --user

2) From PyPI repository ::

    pip install pyqchem --user

Q-Chem setup
------------

PyQchem checks two environment variables to locate Q-Chem installation: $QC and $QCSCRATCH.
*$QC* contains the path to the root directory of Q-Chem installation. This directory
should contain a */bin* directory where *qchem* run script is placed. *$QCSCRATCH* contains
the path to scratch directory.

.. note::
    $QCAUX and $QCREF environment variables are not directly used by pyqchem but they are
    used by Q-Chem to properly run. Check Q-Chem manual for further information.


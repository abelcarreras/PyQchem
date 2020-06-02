.. highlight:: rst

Symmetry
========

PyQchem is combined with wfnsympy to analyze the symmetry of the wave function calculated using Q-Chem.
PyQchem implements several symmetry functions that are compatible with *electronic_structure* dictionary.
Due to the limitations of *wfnsympy* at this moment only gaussian type orbitals (GTO) basis can be used.

Orbital classification
----------------------

A simple function included in pyQchem is *get_orbital_classification*, this function classify the molecular
orbitals between PI and SIGMA. For this function to work properly the user needs to define the center ond
orientation of the molecule. Assuming that the molecule is planar *center* should be a point within the plane
that the atoms form and *orientation* is a unitary vector perpendicular to this plane.

The return of this function is a list. Each element of this list correspond to a molecular orbital (in energy order
from lower to higher) and contains two things: a label that indicates the type of orbitals (SIGMA or PI) and a float
number that indicates the degree of accuracy (from 0[None] to 1[Full]).

..  code-block:: python

    from pyqchem.symmetry import get_orbital_classification
    orbital_types = get_orbital_classification(electronic_structure,
                                               center=[0.0, 0.0, 0.0],
                                               orientation=[0.0, 0.0, 1.0])

    for i, ot in enumerate(orbital_types):
        print('{:5}:  {:5} {:5.3f}'.format(i + 1, ot[0], ot[1]))


Sometimes molecules may not be planar but still some notion of pi/sigma symmetry can be extracted,
even if its only local. For this reason *pyqchem* implements several fuctions to manipulate electronic
structures. See **classify_orbitals.py** and **advanced_symmetry.py** as more complete examples.

.. automodule:: pyqchem.utils
    :members:


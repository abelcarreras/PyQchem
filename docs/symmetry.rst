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
even if its only local. For this reason *pyqchem* implements several functions to manipulate electronic
structures (See **classify_orbitals.py** and **advanced_symmetry.py** as more complete examples)

The following example shows the use of two of these functions (*get_plane*, and *crop_electronic_structure*) to
determine the orbital symmetry (PI/SIGMA) of a fragment of a large molecule.
*Get_plane* allows to determine
the plane and orientation of the fragment. This function assumes that all atoms are in the same plane, so if
some atoms are out of the plane it may be more adequate to use this function with the subset of the atoms that
are more or less in a plane.
On the other hand, *crop_electronic_structure* modifies the MO coefficients, by setting all the basis functions
that are not centered in the atoms of the fragment to zero. This allows to do a symmetry measure of the part
of the MO orbitals that is located around the fragment.

..  code-block:: python

    # define the atoms of the fragment
    atoms_list = [0, 1, 2, 3, 4, 5]

    # get coordinates of the fragment
    coord_fragment = electronic_structure['structure'].get_coordinates(fragment=atoms_list)

    # get the plane and orientation of the fragment
    center, normal = get_plane(coord_fragment)

    # Set zero to all coefficients centered in the atoms that are not part of the fragment
    electronic_structure_fragment = crop_electronic_structure(electronic_structure, atoms_list)

    # get classified orbitals
    orbital_types = get_orbital_classification(electronic_structure_fragment,
                                               center=center,
                                               orientation=normal)


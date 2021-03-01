.. highlight:: rst

Introduction
============

PyQchem in a python interface for Q-Chem, a popular general-purpose
quantum chemistry maintained and distributed by Q-Chem, Inc., located
in Pleasanton, California, USA.

PyQchem allows to take advantage of Python's simple and powerful syntax
to automatize Q-Chem calculations. For this purpose PyQchem implements
an input generation class, a calculation submitting function and a set of
flexible parsers that extracts the output information and converts in
a well structured python dictionary. These parsers are intended to be as
homogeneous as possible among the different methods producing a
python dictionary that contains similar entries and can be used with the
same analysis/visualization functions.

The phylosophy of PyQChem is to build a homogenious input and output
interface for the different methods implemented in Q-Chem to make the
life of Q-Chem users easier.

Main features
-------------

- Easy to use & clean python interface
- Easy to install in your personal computer or cluster, no special q-chem compilation needed.
- Output parser support for a variety of calculation.
- Modular implementation that allows to easily extend its functionality by writing new parsers and analysis functions.


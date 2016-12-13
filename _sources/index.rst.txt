.. ns_doc documentation master file, created by
   sphinx-quickstart on Thu Jul 16 23:02:52 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pymatnest documentation
==================================

The ``pymatnest`` package is a software
library written by Noam Bernstein in Fortran 95/python for the purpose of carrying out
nested sampling calculations.

Long term support of the package is ensured by:
 - Noam Bernstein (Naval Research Laboratory)
 - Gabor Csanyi (University of Cambridge)
 - Rob Baldock  (University of Cambridge)
 - Livia Bartok-Partay (University of Cambridge)
 
Copyright 2015-2016.

Most of the publicly available version is in the public
domain.

Features
========

* both constant volume and constant pressure calculations (variable cell shapes)
* single and multicomponent systems
* Monte Carlo and Molecular Dynamics exploration
* can be used with

   - the supplied fortran models

   - ``LAMMPS``

   - ``QUIP``

Contents
========
.. toctree::
   :maxdepth: 3

   intro.rst
   install.rst
   getstart.rst
   ns_run.rst
   ns_rng.rst
   rngstream.rst
   lammpslib.rst
   stacktrace.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


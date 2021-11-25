.. ns_doc documentation master file, created by
   sphinx-quickstart on Thu Jul 16 23:02:52 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

About ``pymatnest``
==================================

The ``pymatnest`` package is a software
library written by Noam Bernstein in Fortran 95/python for the purpose of carrying out
nested sampling calculations.

Long term support of the package is ensured by:
 - Noam Bernstein (Naval Research Laboratory)
 - Gabor Csanyi (University of Cambridge)
 - Rob Baldock  (EPFL)
 - Livia Bartok-Partay (University of Warwick)
 
Copyright 2015-2016.

Most of the publicly available version is in the public
domain.
If you use this code please cite the following articles:

R.J.N. Baldock, N. Bernstein, K. M. Salerno, L.B. Partay, G. Csanyi, *Constant-pressure nested sampling with atomistic dynamics*,
Phys. Rev. E (2017), 96, 043311, http://link.aps.org/doi/10.1103/PhysRevE.96.043311

L.B. Partay, G. Csanyi, N. Bernstein *Nested sampling for materials*
Eur. Phys. J. B (2021) 94, 159, https://link.springer.com/article/10.1140/epjb/s10051-021-00172-1

Features
========

* both constant volume and constant pressure calculations (variable cell shapes)
* single and multicomponent systems
* Monte Carlo and Molecular Dynamics exploration
* 2D systems
* can be used with

   - the supplied fortran models

   - ``LAMMPS``

   (- ``QUIP`` not currently compatible)

Contents
========
.. toctree::
   :maxdepth: 3

   intro.rst
   install.rst
   getstart.rst
   tutorial.rst
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


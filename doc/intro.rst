Nested sampling method 
==================================

:Authors: Noam, Rob, Gabor, Livia
:Date: 2015


Overview
-------------------------------------

The nested sampling method was first introduced by John Skilling[#f1]_ 
in the field of applied probability and inference to efficiently sample 
probability densities in high-dimensional spaces where the regions contributing 
most of the probability mass are exponentially localized.
It is an annealing algorithm that creates a sequence of probability distributions, 
each more concentrated than the previous one near the high-likelihood region of the 
parameter space - i.e. in the context of materials, the low-energy region of configuration space.


Since its original inception, nested sampling has also been applied to atomistic systems[#f2]_ [#f3]_,
and its several advantages mean it became a powerful method to uniformly sample atomic 
configuartion spaces.
 * calculate the partition function, hence all thermodynamic properties become accesible
     - calculate the heat capacity and locate the phase transitions
     - with an order parameter calculate the free energy
 * the sampling process itself is independent of temperature, thus the calculation of thermodynamic quantities is a simple post-processing step
 * no prior knowledge is needed of the phases or phase transitions
 * can be used with both constant volume and constant pressure calculations
 * calculate the entire phase diagram in an automated way
 * considerable computational gain over parallel tempering

.. [#f1] J. Skilling, in Nested Sampling, edited by Rainer Fischer, Roland Preuss, and Udo von Toussaint, AIP Conf. Proc. No. 735 (AIP, New York, 2004), p. 395
.. [#f2] L. B. Partay, A. P. Bartok, and G. Csanyi, J. Phys. Chem. B 114, 10502 (2010)
.. [#f3] L. B. Partay, A. P. Bartok, and G. Csanyi, Phys. Rev. E 89, 022302 (2014)

Algorithm
-------------------------------------

Nested sampling is an iterative method...
Top-down, Size of live set, random walk...etc.


``pymatnest``
-------------------------------------

The ``pymatnest`` package is a software
library written in Fortran 95/python for the purpose of carrying out
nested sampling calculations with a variatey of options suitable for different systems. 
It can be used with the in-built fortran potential models and it also has interfaces with the following packages:

   - ``LAMMPS``

   - ``QUIP``

MC and MD step algorithms
+++++++++++++++++++++++++++++++++++++

.. include:: MC_MD_steps.txt


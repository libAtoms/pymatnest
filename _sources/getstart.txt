Getting started               
==================================

List of input parameters
++++++++++++++++++++++++++++++++
.. automodule:: ns_run
   :members: usage

Example input files
+++++++++++++++++++++++++++++++

Cluster with ``QUIP`` implemented Lennard-Jones potential, using MD
--------------------------------------------------------

This input file will start a calculation of 6 Lennard-Jones atoms 
in a cubic cell with volume 648 Angstrom^3, using the potential implemeted in ``QUIP``
and using MD trajectory for generating a new sample configuration. As the shape 
and volume of the cell will not be changed, a cluster will be formed by the atoms.

.. literalinclude:: ../tests_quip/inputs.test.cluster.MD.quip
    :language: python

The corresponding ``quip_params.cluster.xml`` file is the following:

.. literalinclude:: ../tests_quip/quip_params.cluster.xml

Periodic Lennard-Jones system with ``LAMMPS``, using MD
--------------------------------------------------------

This input file will start a calculation of 64 Lennard-Jones atoms 
in a cell with variable shape and size at the set pressure value, using the potential implemented in ``LAMMPS``
and using MD trajectory for generating a new sample configuration. 

.. literalinclude:: ../tests_lammps/inputs.test.periodic.MD.lammps
    :language: python


Restart a run
++++++++++++++++++++++++++++++++

If one needs to continue a run, i.e. perform more iteration cycles, the input file has to be modified as

::

    #start_species=1 64 1.0
    restart_file=my_output.snapshot.all.extxyz


Comment out the ``start_species`` option and include a keyword ``restart_file``. This defines the name of the file 
where all the walker configurations can be read from. This file should be a produced by
concatenating the last saved snapshot files of all the processors. You might also need to increase the ``n_iter_per_walker``
value if the run already reached the number of iterations set previously.



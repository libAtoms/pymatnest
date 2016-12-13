Getting started               
==================================

List of input parameters
++++++++++++++++++++++++++++++++
.. automodule:: ns_run
   :members: usage

Example input files
+++++++++++++++++++++++++++++++

Cluster with the supplied fortran code Lennard-Jones potential, using MC
-------------------------------------------------------------------------

This input file will start a calculation of 6 Lennard-Jones atoms 
in a cubic cell with volume 648 Angstrom^3, using the potential implemented in the supplied fortran routine,
and using MC trajectory for generating a new sample configuration. As the shape 
and volume of the cell will not be changed, a cluster will be formed by the atoms.
Note that the input parameters given below, such as the number of walkers and length of trajectory, 
are for a short test run. For six atoms this calcualtion should take only a couple of minutes on one processor, and although
the calculation should end up in the global minimum octahedral cluster structure, 
thermodynamic variables will not be necessarily well converged.

.. literalinclude:: ../example_inputs/inputs.test.cluster.MC.fortran
    :language: python

Cluster with ``QUIP`` implemented Lennard-Jones potential, using MD
--------------------------------------------------------------------

This input file will start a calculation of 6 Lennard-Jones atoms 
in a cubic cell with volume 648 Angstrom^3, using the potential implemented in ``QUIP``
and using MD trajectory for generating a new sample configuration. As the shape 
and volume of the cell will not be changed, a cluster will be formed by the atoms.

.. literalinclude:: ../example_inputs/inputs.test.cluster.MD.quip
    :language: python

The corresponding ``quip_params.cluster.xml`` file is the following:

.. literalinclude:: ../example_inputs/quip_params.cluster.xml

Periodic Lennard-Jones system with ``LAMMPS``, using MD
--------------------------------------------------------

This input file will start a calculation of 64 Lennard-Jones atoms 
in a cell with variable shape and size at the set pressure value, using the potential implemented in ``LAMMPS``
and using MD trajectory for generating a new sample configuration. 

.. literalinclude:: ../example_inputs/inputs.test.periodic.MD.lammps
    :language: python

Periodic binary Lennard-Jones system with the supplied fortran code or``LAMMPS``, using MD
---------------------------------------------------------------

This input file will start a calculation of a binary Lennard-Jones system, with 32 A-type atoms and 32 B-type atoms.
Note that in case of multicomponent systems swap moves have to be introduced, when the coordinates of different atomic types are swapped.
The cell is of variable shape and size at the set pressure value.

.. literalinclude:: ../example_inputs/inputs.test.periodic_binary.MD.FORTRAN_LJ
    :language: python

Polymer with ``LAMMPS``, using MD
---------------------------------------------------------------

Start a nested sampling simulation with a polymer. Initially the random placement of atoms has to be turned off, 
and a configuration file has to be read. The initial walk is done with a heuristic choice for E_max.
The ``atom_style full`` used in this example is part of the MOLECULE package, you have to include that when compling
``LAMPS``.

.. literalinclude:: ../example_inputs/inputs.test.cluster.MD.lammps.polymer
    :language: python

The corresponding ``test_start_config_polymer.xyz`` file is the following:

.. literalinclude:: ../example_inputs/test_start_config_polymer.xyz

Some tips on setting the input parameters
+++++++++++++++++++++++++++++++++++++++++

Minimum lattice height: ``MC_cell_min_aspect_ratio``
----------------------------------------------------

The fully flexible cell is introduced to remove the finite size effect whereby 
it may not be possible to arrange a fixed number of particles
in a fixed shape cell into certain crystal structures. In unfortunate cases this can exclude
thermodynamically relevant structures from the results of the calculation.
But for a flexible cell and for a finite number of particles, there exist simulation cells such that parallel
faces of the cell are separated by only a few layers of atoms, and those does not approximate the infinite fluid
in three dimensions. This problem can be solved by the introduction of a “minimum cell height” parameter.

The effect of the chosen minimum cell height is demonstrated in the case of the periodic system of 
64 Lennard-Jonesium particles in the figure below. The legend on the right shows the value of 
``MC_cell_min_aspect_ratio`` used in each simulation. The peak at lower temperature corresponds to melting,
and the peak at higher temperature to evaporation. The location of
the evaporation transition is converged for ``MC_cell_min_aspect_ratio`` ≥ 0.35, but melting requires
a higher ``MC_cell_min_aspect_ratio`` ≥ 0.65.  At low values of ``MC_cell_min_aspect_ratio`` the
system’s behaviour is dominated by the fictitious periodicity imposed by the boundary
conditions.
(source: R. J. N. Baldock, *Classical Statistical Mechanics with Nested Sampling*, Ph.D. thesis, University of Cambridge (2014).)

.. figure:: doc_figure_mlh.jpg
   :align: center
   :width: 500


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



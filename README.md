
The ``pymatnest`` package is a software library for carrying out nested sampling calculations. 
It can be used to explore the energy landscape of different materials, calculate thermodynamic variables at
arbitrary temperatures, locate phase transitions and calculate the phase diagram. It can be used with LAMMPS,
QUIP, and the supplied fortran model, and both with MC and MD.

Publications:

L.B. Partay, A.P. Bartok, G. Csanyi, *Efficient Sampling of Atomic Configurational Spaces*, 
J. Phys. Chem. B (2010), 114, 10502–10512, http://pubs.acs.org/doi/abs/10.1021/jp1012973

L.B. Partay, A.P. Bartok, G. Csanyi, *Nested sampling for materials: The case of hard spheres*, 
Phys. Rev. E (2014), 89, 022302, http://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.022302

R.J.N. Baldock, L.B. Partay, A.P. Bartok, M.C. Payne, G. Csanyi, *Determining pressure-temperature phase diagrams of materials*,
Phys. Rev. B (2016), 93, 174108, http://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.174108 

******

System requirements
------------------------------------------------------------------------------

You have to have ``numpy`` and ``ASE`` installed (check that you have a newer ASE
version otherwise the extended xyz version will not be recognised). 
You have to have ``mpi4py`` too.


Compilation of the streamlined FORTRAN models and MC/MD walkers.
------------------------------------------------------------------------------

1. edit ``Makefile.arch`` (by default Linux appropriate variables are set)
2. ``make``


Using with ``QUIP``
------------------------------------------------------------------------------

Make sure your ``PYTHONPATH`` is set correctly to find the ``quippy`` module.
(Note: Check the ``QUIP_ARCH`` with which you build quippy, as some might result in segmentation fault when running
``ns_run`` (e.g. try gfortan_openmpi if gfortran fails)


Using with ``LAMMPS``
------------------------------------------------------------------------------

These instructions assume the latest (git/svn) version of ``LAMMPS``.  It is not tested how
far back older versions would also work.

Apply the ``MPI`` patch to the ``LAMMPS`` source by doing

    ``cd lammps_top_dir/src``

    ``patch < ns_run_dir/lammps.patch``

where ``ns_run_dir`` is the directory where ``ns_run`` is, and ``lammps_top_dir`` is the ``LAMMPS`` directory.
Create a Makefile for **parallel** lammps in ``lammps_top_dir/src/MAKE``, and define ``-DLIBRARY_MPI_COMM_WORLD=MPI_COMM_SELF`` 
in the ``LMP_INC`` makefile variable. Then

    ``make [machine] mode=shlib``

Copy ``lammps_top_dir/python/lammps.py`` to the directory set in your ``PYTHONPATH``.

Copy ``lammps_top_dir/src/liblammps_[machine].so`` to the same place where you copied ``lammps.py``.

**Important note:** Check the ``lammps.py`` file as the path definition used to have a bug in the line:

``else: self.lib = CDLL(join(modpath,"/liblammps_%s.so" % name),RTLD_GLOBAL)`` 

You HAVE TO delete the ``/`` before ``liblammps`` otherwise it is meant as an absolute path!!!

``LAMMPS_name`` is what you set for ``[machine]`` when compiling ``LAMMPS``.

Note that you **have** to compile a parallel version of ``LAMMPS`` with the source patch provided in ``pymatnest``.  
``LAMMPS`` "serial" compilation still links to fake ``MPI`` routines, which then conflict in unpredictable ways with 
the true mpi routines that ``mpi4py`` includes.

It is possible to use ``OpenMP`` to parallelize each ``LAMMPS`` task.  This has been tested to run, but not for correctness or efficiency.

    ``cd lammps_top_dir/src``

    ``make yes-user-omp``

Add openmp enabling flag (e.g. ``-fopenmp`` for gfortran) to ``CCFLAGS`` in the ``MAKE/MINE/Makefile.[machine]``

    ``make [machine] mode=shlib``

Copy ``liblammps_[machine].so`` as before.

When running:

Set ``OMP_NUM_THREADS`` environment variable for number of ``OpenMP`` threads per task, and
add ``LAMMPS_header_extra='package omp 0'`` input file argument.
Use ``LAMMPS`` pair style that activates omp, e.g. ``pair_style lj/cut/omp 3.00``.
Pass flags to ``mpirun`` to tell it to run fewer MPI tasks than total number of cores assigned to entire job so that cores are 
available for OpenMP parallelization.
Example for OpenMPI, on 8 nodes, with 16 cores each, OpenMP parallelizing each MPI task's ``LAMMPS`` work over all 16 cores:

     ``export OMP_NUM_THREADS=16``

     ``mpirun -np 8 -x OMP_NUM_THREADS --map-by slot:pe=$OMP_NUM_THREADS ns_run < inputs``

Note: the ``-np 8`` may not be needed, depending on your queueing system. 
The ``LAMMPS ASE`` interface (``lammpslib.py``) is a heavily modified version of

<https://svn.fysik.dtu.dk/projects/ase-extra/trunk/ase/calculators/lammpslib.py>

For more information on how the interface works, see the :any:`lammpslib`.

Running 
------------------------------------------------------------------------------

To start a nested sampling run type

   ``ns_run < input``

Example input files can be found in the folder ``./example_inputs``.

For further help see also

   ``ns_run --help``

If you get weird errors about modules and/or ``.so`` files not found, do (in sh syntax)

   ``export PYTHONPATH=ns_run_dir:$PYTHONPATH``

where ``ns_run_dir`` is the directory where ``ns_run`` is.
This appears to be necessary on some HPC machines where mpirun copies the executable,
because ``ns_run`` by default looks for modules in the same directory as the top level 
python script itself. If it is still not sufficient, you might have to copy the entire ``ns_run_dir``
to the directory where the jobs are submitted from.

Running on ARCHER (UK National Supercomputing Service)
------------------------------------------------------------------------------

Install the latest ``ASE`` (3.9 or later) version and add that directory to your ``PYTHONPATH``, as the 
default version on ARCHER is just 3.8.

Copy the whole ``pymatnest`` library to your ``/work`` directory, otherwise the compute nodes will not be
able to read all the relevant python files.

In the job script you have to swap and load appropriate modules.

   ``module load python-compute``

   ``module load pc-numpy``

   ``module load gcc``


Analysis
------------------------------------------------------------------------------

To analyse the results you can use

   ``ns_analyse -M 0.01 -D 0.01 -n 100 file.energies > analysis``

For further help see also

   ``ns_analyse --help``


Temperature averaged analysis workflow
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Merge configurations using
   ``ns_process_traj -t``

Do analysis on output of ``ns_process_traj`` using ``structure_analysis_traj``.

Add T-dependent weights to analyses using ``ns_set_analysis_weights``.  This will write new analysis files, one per temperature per analysis, with do_weights set in the header and each data line prepended by the weight.

Finally, use ``mean_var_correl`` to calculated the weighted mean of each analysis at each temperature.

About the documentation
------------------------------------------------------------------------------

The documentation with example input files and a list of keywords...etc. can be found at
<http://libatoms.github.io/pymatnest/>.

The documentation is generated by Sphinx, using the files within the ``doc`` library.
Modules are autodocumented with ``.. automodule::`` so all the properly formatted comments
in the python code (i.e. within triple quote) appear.
The installation and basic usage guidelines in the documentation are shown as the content of the README.md file
is ``.. included:``-d.
Example inputs are located in the folder ``./example_inputs`` and these files are also included in the documentation together with additional comments. 


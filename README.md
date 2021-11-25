
The ``pymatnest`` package is a software library for carrying out nested sampling calculations. 
It can be used to explore the energy landscape of different materials (both clusters and bulk materials), 
calculate thermodynamic variables at arbitrary temperatures, locate phase transitions and calculate the 
phase diagram. It can be used with LAMMPS, and the supplied fortran models, and both with MC and MD.

If you use pymatnest, please cite the following publications: 
(references are available in bibtex format in the ``NS_publications.bib`` file)

L.B. Partay, A.P. Bartok, G. Csanyi, *Efficient Sampling of Atomic Configurational Spaces*, 
J. Phys. Chem. B (2010), 114, 10502–10512, http://pubs.acs.org/doi/abs/10.1021/jp1012973

R.J.N. Baldock, L.B. Partay, A.P. Bartok, M.C. Payne, G. Csanyi, *Determining pressure-temperature phase diagrams of materials*,
Phys. Rev. B (2016), 93, 174108, http://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.174108 

R.J.N. Baldock, N. Bernstein, K. M. Salerno, L.B. Partay, G. Csanyi, *Constant-pressure nested sampling with atomistic dynamics*,
Phys. Rev. E (2017), 96, 043311, http://link.aps.org/doi/10.1103/PhysRevE.96.043311

Further publications using the pymatnest package:

L.B. Partay *On the performance of interatomic potential models of iron: Comparison of the phase diagrams*,
Comp. Mat. Sci. (2018), 149, 153, https://www.sciencedirect.com/science/article/pii/S0927025618301794

J. Dorrell, L.B. Partay *Thermodynamics and the potential energy landscape: case study of small water clusters*,
Phys. Chem. Chem. Phys. (2019), 21, 7305 https://pubs.rsc.org/en/content/articlehtml/2019/cp/c9cp00474b

******

System requirements
------------------------------------------------------------------------------

You have to have ``numpy`` and ``ASE`` installed (check that you have a newer ASE
version otherwise the extended xyz version will not be recognised). 
You have to have ``mpi4py`` and ``psutil`` too.


Compilation of the streamlined FORTRAN models and MC/MD walkers.
------------------------------------------------------------------------------

1. edit ``Makefile.arch`` (by default Linux appropriate variables are set)
2. ``make``


[//]: # (Using with ``QUIP`` - NO LONGER COMPATIBLE!!!)

[//]: # (Make sure your ``PYTHONPATH`` is set correctly to find the ``quippy`` module.)

[//]: # (Note: Check the ``QUIP_ARCH`` with which you build quippy, as some might result in segmentation fault when running)

[//]: # (``ns_run`` (e.g. try gfortan_openmpi if gfortran fails))


Using with ``LAMMPS``
------------------------------------------------------------------------------

These instructions assume the latest (git/svn) version of ``LAMMPS``.  It is not tested how
far back older versions would also work.

### Nearly mandatory compilation flags (for all ``LAMMPS`` versions)

It is necessary to take advantage of the ``-DLAMMPS_EXCEPTIONS``
flag, which allows lammps crashes to be handled gracefully within python.  Add it to the ``LMP_INC`` variable in the
``LAMMPS`` makefile before compiling.

### Basic instructions for recent versions of ``LAMMPS`` and ``mpi4py`` version 2.0.0 or newer

Create an appropriate ``LAMMPS`` **serial** makefile, and compile with 

- ``make [machine] mode=shlib``

For example ``make serial mode=shlib``.
Install the python files supplied with ``LAMMPS`` :                                    

- ``make install-python``

Copy ``lammps_top_dir/src/liblammps_[machine].so`` to the same place where ``LAMMPS`` installed the python packages. (Likely to be in the PYTHONAPTH, where a lammps directory is created). 

The input file variable ``LAMMPS_name`` is what you set for ``[machine]`` when installing ``lammps_[machine].so``.
By default it is what you set for ``machine`` when compiling ``LAMMPS``, unless the library was renamed when installing.

### Basic instructions for ``LAMMPS`` versions Nov 2020 or before: 

Create an appropriate ``LAMMPS`` **parallel** makefile, and compile with 

- ``make [machine] mode=shlib``

Copy ``lammps_top_dir/python/lammps.py`` to a directory set in your ``PYTHONPATH``.

Copy ``lammps_top_dir/src/liblammps_[machine].so`` to the same place where you copied ``lammps.py``.

The input file variable ``LAMMPS_name`` is what you set for ``[machine]`` when installing ``lammps_[machine].so``.
By default it is what you set for ``machine`` when compiling ``LAMMPS``, unless the library was renamed when installing.
Set the input variable ``LAMMPS_serial=F``. 


### Support for GMC within LAMMPS

Copy the two GMC-related files ``ns_run_dir/lammps_patches/fix_gmc.*`` to the ``LAMMPS`` directory ``lammps_top_dir/src/`` 
before compiling, and set ``LAMMPS_fix_gmc=T`` in the input file.

### Support for molecules

For using molecules use the latest version of LAMMPS. If you need, copy the improper related files ``ns_run_dir/lammps_patches/create_*`` to the ``LAMMPS`` directory ``lammps_top_dir/src/`` 
before compiling.  See the file ``example_inputs/inputs.test.cluster.MD.lammps.polymer`` for an example.

### Mixed ``MPI-OpenMP``

It is possible to use ``OpenMP`` to parallelize each ``LAMMPS`` task.  This has been tested to run, but not for correctness or efficiency.

- ``cd lammps_top_dir/src``
- ``make yes-user-omp``

Add openmp enabling flag (e.g. ``-fopenmp`` for gfortran) to ``CCFLAGS`` in the ``MAKE/MINE/Makefile.[machine]``, then compile and install
as above.

When running:

- Set ``OMP_NUM_THREADS`` environment variable for number of ``OpenMP`` threads per task, and
- add ``LAMMPS_header_extra='package omp 0'`` input file argument.
- Use ``LAMMPS`` pair style that activates omp, e.g. ``pair_style lj/cut/omp 3.00``.
- Pass flags to ``mpirun`` to tell it to run fewer MPI tasks than total number of cores assigned to entire job so that cores are available for OpenMP parallelization.
- Example for OpenMPI, on 8 nodes, with 16 cores each, OpenMP parallelizing each MPI task's ``LAMMPS`` work over all 16 cores:

     - ``export OMP_NUM_THREADS=16``

     - ``mpirun -np 8 -x OMP_NUM_THREADS --map-by slot:pe=$OMP_NUM_THREADS ns_run < inputs``

Note: the ``-np 8`` may not be needed, depending on your queueing system. 

### Other notes

The ``LAMMPS ASE`` interface (``ns_run_dir/lammpslib.py``) is a heavily modified version of

<https://svn.fysik.dtu.dk/projects/ase-extra/trunk/ase/calculators/lammpslib.py>

For more information on how the interface works, see the :any:`lammpslib`.

### For versions of ``mpi4py`` older than 2.0.0

If you have ``mpi4py`` version older than 2.0.0, you will need to patch LAMMPS as follows.

Apply the communicator patch to the ``LAMMPS`` source by doing

- ``cd lammps_top_dir/src``
- ``patch < ns_run_dir/lammps_patches/communicator_self.patch``

where ``ns_run_dir`` is the directory where ``ns_run`` is, and ``lammps_top_dir`` is the ``LAMMPS`` directory.
Create a Makefile for **parallel** lammps in ``lammps_top_dir/src/MAKE``. 
Define ``-DLIBRARY_MPI_COMM_WORLD=MPI_COMM_SELF`` in the ``LMP_INC`` makefile variable, then compile
as above.

### For older versions of ``LAMMPS``

**Important note:** Check the ``lammps.py`` file as the path definition used to have a bug in the line:

``else: self.lib = CDLL(join(modpath,"/liblammps_%s.so" % name),RTLD_GLOBAL)`` 

You HAVE TO delete the ``/`` before ``liblammps`` otherwise it is interpreted as an absolute path!!!


Running 
------------------------------------------------------------------------------

To start a nested sampling run type

   ``ns_run < input``

When running, it is strongly recommendded you set the ``OMP_NUM_THREADS=1`` environment variable (e.g. in your jobscript) to avoid
multiple ``OpenMP`` threads starting which can seriosly slow down the calculations (unless you have compiled ``LAMMPS`` to be used
with mixed ``MPI-OpenMP``). 

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
-------------------------------------------------------------------------------

This assumes that QUIP is installed (``structure_analysis_traj`` and ``mean_var_correl`` are part of QUIP).

Merge configurations using
   ``ns_process_traj -t``

Do analysis on output of ``ns_process_traj`` using ``structure_analysis_traj``.

Add T-dependent weights to analyses using ``ns_set_analysis_weights``.  This will write new analysis files, one per temperature per analysis, with ``do_weights`` set in the header and each data line prepended by the weight.

Finally, use ``mean_var_correl`` to calculated the weighted mean of each analysis at each temperature.

**Automatic script using QUIP ``ns_process_traj`` and ``structure_analysis_traj``:**

``make_thermal_average_xrd_rdfd_lenhisto.py`` is a script for calculating thermally averaged powder spectra (``(...)_xrd``), radial distribution functions (``(...)_rdfd``), which are currently disabled (see below), and histograms of lattice vector lengths (``(...)_lattice_len_histo``).
RDFDs and XRDs are calculated for reference structures and safed under ``$STRUCTURE_NAME_V_mean_of_$TRAJ_signifpart_$SIGNIFICANT_PART.T_$T_xrd`` and ``$STRUCTURE_NAME_V_mean_of_$TRAJ_signifpart_$SIGNIFICANT_PART.T_$T_rdfd``.
It calculates the weights on its own and can deal with single trajectory files as well as combined trajectory files.

Before using, QUIP and quippy need to be installed and the variable ``QUIP_path`` in ``make_thermal_average_xrd_rdfd_lenhisto.py`` line 28 must be set to the QUIP build directory.

**Important note:** Only one script can be active in a single folder at a given time. Otherwise, temporary files will be overwritten and the results incorrect.

The script is called via:

``python make_thermal_average_xrd_rdfd_lenhisto.py -fn traj.extxyz -Ts "600 800 1000" -nc 8 -nw 1920 -sn "bcc fcp hcp" -sc "test_struc_1.xyz test_struc_2.xyz``

- ``-fn`` is the file name. traj.extxyz can be a combined or a single trajectory.
- ``-Ts`` are the different temperatures (which are transformed to integers) in the format "T_1 T_2 ... T_N-1 T_N".
- ``-nc`` is the number of culled walkers per iteration.
- ``-nw`` is the number of walkers.
- ``-sn`` are the names of structures (defined in misc_calc_lib.py) for xrd spectrum identification in format 'struc_name_1 struc_name_2 ... struc_name_N-1 struc_name_N'. Only works for single species configurations.
- ``-sc`` are the paths to the ``.extxyz``/``.xyz`` files of reference structures in format 'path_1 path_2 ... path_N-1 path_N'.

The following variables set in the script may be intersting:

**significant_part**

The parameter ``significant_part`` controls how much of the sampled structures we actually consider. It follows the name
``_signifpart_`` in the filename. For example, if it was set to 0.25 we would only consider the ca 25% most likely structures. (Due to discrete weight steps, this number is not exact.) The default value of ``significant_part`` is 0.95. This ignores irrelevant structures and especially excludes high volume systems when we consider the solid phases. (To speed up the calculations one could go lower, but without further experimentation, no clear recommendations can be made with regards to this.)

**do_rdfd**

``do_rdfd = False`` controls whether radial density functions are calculated. RDFs in QUIP are not using periodic cells. This makes it very hard to compare different cells of the same structure. Hence, it is turned off. If set to ``True``, the script uses a 6x6x6 supercell for the comparison structures.


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


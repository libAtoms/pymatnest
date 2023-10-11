import re, math, time, os
import pprint
import numpy as np, ase, ase.io
import ns_rng
import stacktrace
from copy import deepcopy
import pick_interconnected_clump
try:
    import matscipy.neighbours
except:
    pass
import collections
from traceback import print_exception

import check_memory
from ase.md.verlet import VelocityVerlet
import importlib

print_prefix=""

def usage():
    """ Print help to the standard output about the usage of the code and input parameters. The current list of parameters is the following:

.. glossary::

    ``max_volume_per_atom=float``
       | Maximum volume per atom allowed during the run.
       | default: 1.0e3

    ``min_volume_per_atom=float``
       | Minimum volume per atom allowed during the run.
       | default: 1.0

    ``start_species=int int [ float ] [, int int [ float ] ... ]``
       | MANDATORY (or start_config_file)
       | Atomic number; multiplicity; [ not recommended: mass (amu) ]. Info repeated for each species, separated by commas, mass is optional and not recommended.

    ``start_config_file=str``
       | MANDATORY (or start_species)
       | Name of file to read in for atom information
       | default: ''

    ``keep_atoms_fixed=int ``
       | Number of atoms to fix in surface simulation. 

    ``apply_Z_wall=[T|F] ``
       | Whether to have a boundary in the Z direction to keep free particles reaching the other side of "surface" layer due to periodic boundary conditions.
       | This functionality is not fully tested! Recommended use with MC evaluator and fortran. If constructing a surface layer, make it parallel to the XY plane and set the X dimension with the "wall" taken into account, no atoms should violate the wall restriction initially!!!
       | If True, the wall is set at 10 angstroms below the top of the simulation cell.

    ``restart_file=path_to_file``
       | File for restart configs. Mutually exclusive with ``start_*``, one is required. The file should contain the state of the walkers to continue from along with the restart iteration number. Normally such a file can be the concatenated snapshot files.

    ``n_walkers=int``
       | MANDATORY
       | Total number of walkers, i.e. the size of the live set used. It has to be the multiple of the number of processors used and it also has to be larger than the number of processors.

    ``n_cull=int``
       | Number of walkers to kill at each NS iteration. Use of default 1 is strongly recommended.
       | deafult: 1

    ``n_extra_walk_per_task=int``
       | default: 0

    ``n_iter_times_fraction_killed=int``
       | MANDATORY
       | Number of nested sampling iteration cycles performed per walker. Thus the total number of iterations will be ``n_iter_times_fraction_killed`` / ``(n_cull/n_walkers)``. Either this or ``converge_down_to_T`` is required.

    ``converge_down_to_T=flot``
       | MANDATORY
       | temperature down to which Z(T) should be converged.  Either this or ``n_iter_times_fraction_killed`` is required.

    ``T_estimate_finite_diff_lag=int``
       | default: 1000
       | Lag (in iterations) for doing d log Omega / d E finite difference derivative

    ``min_Emax=float``
       | Termination condition based on Emax: if this value is reached, the iteration will stop.
       | No default.

    ``out_file_prefix=str``
       | String used as prefix for the different output files.
       | No default.

    ``energy_calculator= ( ASE | lammps | internal | fortran)``
       | Energy calculator.
       | default: fortran

    ``n_extra_data=int``
       | Amount of extra data per atom to pass around.
       | default: 0

    ``KEmax_max_T=float``
       | If > 0, maximum temperature for estimating KEmax.  For constant *P* ensemble 3/2 P Vmax will be used if KEmax_max_T is < 0.
       | default: -1

    ``kB=float``
       | Boltzmann constant
       | default: 8.6173324e-5 (eV/A)

    ``start_energy_ceiling_per_atom=float``
       | Maximum potential energy per atom for initial configurations.  P*Vmax is added to this automatically in case of NpT runs.
       | default: 1.0e9

    ``start_energy_ceiling=float``
       | DEPRECATED: use ``start_energy_ceiling_per_atom``. Maximum potential energy for initial configurations.  P*Vmax is added to this automatically in case of NpT runs.
       | default: 1.0e9

    ``random_init_max_n_tries=int``
       | Maximum number of tries to create initial random atomic configuration
       | default 100

    ``n_model_calls_expected=int``
       | Total number of model calls performed during the random walk, expected by the user. The actual number of model calls performed during runtime can be different: this number is divided up among the processors in case of running parallel, and in order to make sure different processors perform about the same number of calls, it can be increased for some. (It is recommended to use n_model_calls_expected/[number of processors] > 20) Either this or the keyword ``n_model_calls`` is mandatory, but the use of this keyword is strongly recommended. If < 0, value will be set to _sum_ of each type of steps (accounting for atom_traj_len), so ``n_*_steps`` will be setting the _number_ of calls, rather than just the ratios.
       | default: 0

    ``n_model_calls=int``
       | Number of model calls performed during the random walk. This is the actual number of calls performed, ideally the program sets its value depending on n_model_calls_expected and the number of processors used. Either this or the keyword ``n_model_calls_expected`` is mandatory, though the usage of the latter is strongly recommended.
       | default: 0

    ``do_good_load_balance=[T | F]``
       | Whether to do steps in groups with good load balance
       | default: F

    ``n_atom_steps=int``
       | Ratio of atomic trajectories will be determined by ``n_atom_steps``/SUM(``n_*_steps``).
       | default: 1

    ``atom_traj_len=int``
       | Length of atomic trajectory (MD steps or MC sweeps) in each atomic type step.
       | default: 8

    ``atom_traj_len_cost_multiplier=int``
       | Multiplier for cost of an atomic step (set to 1 for MD, TE-HMC, and SP-MC with O(1) cost moves, N for naive SP-MC)
       | default: 1

    ``break_up_atom_traj=[T | F]``
       | Whether to intersperse ``n_atom_steps`` atomic sub-trajectories with other types of steps.
       | default: F

    ``n_cell_volume_steps=int``
       | Ratio of cell volume steps will be determined by ``n_cell_volume_steps``/SUM(``n_*_steps``).
       | default: 1

    ``n_cell_shear_steps=int``
       | Ratio of cell shear steps will be determined by ``n_cell_shear_steps``/SUM(``n_*_steps``).
       | default: 1

    ``n_cell_stretch_steps=int``
       | Ratio of cell stretch steps will be determined by ``n_cell_stretch_steps``/SUM(``n_*_steps``).
       | default: 1

    ``n_swap_steps=int``
       | Ratio of species swap steps will be determined by ``n_swap_steps``/SUM(``n_*_steps``). It has to be set other than zero for a multicomponent system.
       | default: 0

    ``swap_max_cluster=int``
       | Maximum size of interconnected cluster to try to swap.
       | default: 1

    ``swap_r_cut=float``
       | Cutoff radius for defining connected atoms for cluster.
       | default: 2.5

    ``swap_cluster_probability_increment=float``
       | Factor between prob. of picking increasing larger clusters.
       | default: 0.75

    ``swap_velo=[T | F]``
       | If true, swap velocities when swapping atoms, breaking coherence a bit.
       | default: F

    ``no_swap_velo_fix_mag_alt=[T | F]``
       | If true, use alternate method for correcting velocity magnitudes when not swapping velocities.
       | default: F

    ``n_semi_grand_steps=int``
       | Ratio of species type-change steps will be determined by ``n_semi_grand_steps``/SUM(``n_*_steps``). Only makes sense for a multicomponent system.
       | default: 0

    ``semi_grand_potentials=Z1 : mu1 [, Z2 : mu2 ...]``
       | Chemical potentials for semi-grand canonical species-type transmutations
       | default: ''

    ``velo_traj_len=int``
       | Number of MC steps in (optional) explicit velocity MC trajectory.
       | default: 0

    ``random_energy_perturbation=float``
       | default: 1.0e-12

    ``atom_algorithm=[MC | MD | GMC]``
       | MANDATORY
       | Use either Monte Carlo or Molecular dynamics to explore. GMC is an alias for atom_algorithm=MC, MC_atom_Galilean=True.

    ``MC_atom_velocities=[T | F]``
       | This keyword is supported only for energy_calculator=fortran.
       | default: F

    ``MC_atom_velocities_pre_perturb=[T | F]``
       | Perturb velocities (rejection free) before MC + velocities walk.
       | default: F

    ``MC_atom_step_size=float``
       | Initial atom move step size in units of (max_volume_per_atom * N_atoms)^(1/3)
       | default: 1.0

    ``MC_atom_step_size_max=float``
       | Maximum atom step size in units of (max_volume_per_atom * N_atoms)^(1/3).
       | default: 1.0

    ``MC_atom_uniform_rv=[T | F]``
       | default: F

    ``MC_atom_Galilean=[T | F]``
       | default: F

    ``GMC_no_reverse=[T | F]``
       | default: T

    ``MD_atom_velo_pre_perturb=[T | F]``
       | Perturb velocities before MD trajectory
       | default: F

    ``MD_atom_velo_post_perturb=[T | F]``
       | Perturb velocities after MD trajectory
       | default: T

    ``MD_atom_velo_flip_accept=[T | F]``
       | default: F

    ``MD_atom_timestep=float``
       | default: 0.1 (ASE time units)

    ``MD_atom_timestep_max=float``
       | default: 2.0 (ASE time units)

    ``MD_atom_energy_fuzz=float``
       | Tolerance for rejecting non-energy conserving trajectories, as fraction of Kinetic Energy
       | default: 1.0e-2

    ``MD_atom_reject_energy_violation=[ T | F ]``
       | Use energy conservation violation (exceeding MD_atom_energy_fuzz * KE) to reject MD trajectories.
       | default: F

    ``python_MD=[ T | F ]``
       | Do MD using python code rather than underlying driver
       | default: F

    ``atom_velo_rej_free_fully_randomize=[T | F]``
       | If true, randomize velocities completely rather than just perturbing.
       | default: F

    ``atom_velo_rej_free_perturb_angle=float``
       | Max angle in radians for random rotations.
       | default: 0.3

    ``MC_atom_velo_step_size=float``
       | default: 50.0

    ``MC_atom_velo_step_size_max=float``
       | default: 10000.0

    ``MC_atom_velo_walk_rej_free=[T | F]``
       | default: T. If true, use rejection free algorithm for MC_atom_walk

    ``MC_cell_P=float``
       | Pressure value to be used. The unit of pressure depends on both the energy calculator and on the potential model used. Note that ASE uses
       | eV as energy and Angstrom as distance units everywhere, thus this case the pressure in the input have to be eV/Angstrom^3 (in case of using LAMMPS
       | the lammpslib module will convert the units for lammps to whatever units are set by the potential type)
       | default: 0.0

    ``MC_cell_flat_V_prior=[T | F]``
       | Flat prior for V (use with MC_cell_P=0 and cell moves, requires reweighting configurations when analyzing)
       | POORLY TESTED, DO NOT TRUST (YET).
       | default: F

    ``MC_cell_volume_per_atom_step_size=float``
       | Initial volume stepsize for volume change.
       | default: 5% of the maximum allowed volume

    ``MC_cell_volume_per_atom_step_size_max=float``
       | Maximum allowed volume step size.
       | default: 50% of the maximum allowed volume

    ``MC_cell_volume_per_atom_prob=float``
       | default: 1.0

    ``MC_cell_stretch_step_size=float``
       | default: 0.1

    ``MC_cell_stretch_step_size_max=float``
       | default: 1.0

    ``MC_cell_stretch_prob=float``
       | default: 1.0

    ``MC_cell_shear_step_size=float``
       | default: 0.1, in units of (max_volume_per_atom * N_atoms)^(1/3)

    ``MC_cell_shear_step_size_max=float``
       | default: 1.0, in units of (max_volume_per_atom * N_atoms)^(1/3)

    ``MC_cell_shear_prob=float``
       | default: 1.0

    ``MC_cell_min_aspect_ratio=float``
       | Smallest allowed distance between parallel faces for cell normalised to unit volume. A higher value of ``MC_cell_min_aspect_ratio`` restricts the system to more cube-like cell shapes, while a low value allows the system to become essentially flat. In case of 64 atoms the use of ``MC_cell_min_aspect_ratio`` < 0.65 *does* effect the melting transition.
       | default: 0.8

    ``cell_shape_equil_steps=int``
       | default: 1000

    ``full_auto_step_sizes=[T | F]``
       | If true (T), automatically calibrate all sizes by performing additional short explorations, including at the start of run. If false (F), use initial input step sizes and make small adjustments to these during the run.
       | default: T

    ``monitor_step_interval_times_fraction_killed=float``
     | Divided by ``n_cull/n_walkers`` to get actual monitoring interval in iterations, negative for only using last iteration, 0 for no monitoring
     | default: 1

    ``adjust_step_interval_times_fraction_killed=float``
     | Divided by ``n_cull/n_walkers`` to get actual step adjustment interval in iterations, negative for only using last iteration, 0 for no adjust.
     | default: 1

    ``MC_adjust_step_factor=float``
       |  default: 1.1

    ``MC_adjust_min_rate=float``
       |  default: 0.25

    ``MC_adjust_max_rate=float``
       |  default: 0.75

    ``GMC_adjust_min_rate=float``
       |  default: 0.25

    ``GMC_adjust_max_rate=float``
       |  default: 0.75

    ``GMC_dir_perturb_angle=float``
       |  default: -1.0

    ``GMC_dir_perturb_angle_during=float``
       |  default: 0.0

    ``MD_adjust_step_factor=float``
       |  default: 1.1

    ``MD_adjust_min_rate=float``
       |  default: 0.5

    ``MD_adjust_max_rate=float``
       |  default: 0.95

    ``ASE_calc_module=str``
       |  MANDATORY if energy_calculator=ASE

    ``FORTRAN_model=str``
       |  MANDATORY if energy_calculator=fortran

    ``FORTRAN_model_params=str``
       |  parameters (optional) for fortran model

    ``LAMMPS_fix_gmc=[T | F]``
       | default: F

    ``LAMMPS_init_cmds=str``
       |  MANDATORY if energy_calculator=lammps

    ``LAMMPS_name=str``
       |  '', arch name for lammps shared object file

    ``LAMMPS_serial=[T | F]``
       |  default True, lammps library is serial so do not pass COMM_SELF as a communicator

    ``LAMMPS_header=str``
       |  lammpslib.py value default, override lammpslib.py header commands for energy_calculator=lammps

    ``LAMMPS_header_extra=str``
       |  '', extra lammpslib.py header commands for energy_calculator=lammps

    ``LAMMPS_atom_types=str``
       |  MANDATORY if energy_calculator=lammps
       |  atomic_symbol1 lammps_type1 [, atomic_symbol2 lammps_type2, ...]
       |  mapping between atom species and lammps potential types

    ``config_file_format=str``
       | File format in which configurations are printed, e.g. for trajectory and snapshot files.
       | default: extxyz

    ``rng=( numpy | internal | rngstream )``
       | Random number generator.
       | default: numpy

    ``profile=rank_to_profile``
       | default: -1

    ``2D=[ T | F ]``
       | Perform 2D simulation. Some functionality may be limited with this option!
       | The 2D area is in the xy plane, the z axis stays constant during the simulation (set its value using the keyword ``Z_cell_axis``).
       | The maximum allowed volume is still defined as volume in this case!
       | If constant pressure is set, the p*V term is now p*A*Z_cell_axis, thus set the MC_cell_P parameter multiplied by 1/Z_cell_axis.
       | default: F

    ``Z_cell_axis=float``
       | Only used if 2D=T. Value of the Z cell axis, perpendicular to the area in case of a 2D simulation.
       | This stays constant during the simulation.
       | default: 10.0

    ``debug=int``
       | Verbosity level used in the output file. The larger its value the more info is printed.
       | default: 0

    ``snapshot_interval=int``
       | Iteration interval at which a snapshot is created: every process prints out its current walkers in extended xyz format. If it is set <=0, no snapshots will be printed except the final positions at the end of the nested sampling run. Note that when new snapshots are printed, the previous set is deleted. The snapshot files are convenient source to see how the sampling progresses, but these are also the basis to restart a sampling! When using restart, the walkers will be read from these files.
       | default: 10000

    ``snapshot_time=int``
       | Max time between snapshots in seconds
       | default: 3600

    ``snapshot_per_parallel_task=[ T | F ]``
       | Save a separate snapshot file for each parallel task.  Faster I/O (in parallel), but less convenient in some cases.
       | default: T

    ``snapshot_clean=[ T | F ]``
       | If true, delete previous iteration snapshot files
       | default: T

    ``random_initialise_pos=[ T | F ]``
       | If true, randomize the initial positions
       | default: T

    ``random_initialise_cell=[ T | F ]``
       | If true, randomize the initial cell (currently by random walk)
       | default: T

    ``LAMMPS_molecular_info=[ T | F ]``
       | If true, create lammps bonds and other molecular MM features from initial atoms config (e.g. can be read from a file)
       | default: T

    ``initial_walk_N_walks=int``
       | number of initial random walks to apply to all walkers
       | default: ''

    ``initial_walk_adjust_interval=int``
       | in initial walks, interval between walks for adjusting step size
       | default: ''

    ``initial_walk_Emax_offset_per_atom=float``
       | offset (per atom) to increase Emax during initial random walk applied to all walkers
       | default: 1.0

    ``initial_walk_only=[T | F]``
       | do initial walk only, then quit
       | default: F

    ``traj_interval=int``
     |  Iteration interval at which the currently culled configuration(s) is/are printed to the trajectory output, in the set format. If it is set <=0, no trajectory files will be printed at all. Useful option for larger runs as the trajectory files can become huge.
     |  default: 100

    ``delta_random_seed=int``
     | Random number seed to be used in the run. If smaller than 0, a seed from /dev/urandom is used.
     | default: -1

    ``no_extra_walks_at_all=[ T | F ]``
     | default: F

    ``track_configs=[ T | F ]``
     | Track configrations across all walks/clones.
     | default: F

    ``track_configs_write=[ T | F ]``
     | Write tracked configrations to a file (if false, tracking info will still be in snapshot and traj files)
     | default: F

    ``ns_run_analyzers=string``
     | Analyzers to apply during run.  String consists of semicolon separated pairs of module name and intervals. Module names correspond to analysis modules in NS_PATH/ns_run_analyzers (see there for examples) or elsewhere in PYTHONPATH
     | Positive interval refers to NS loop, negative to initial walks
     | default:''

    """
    sys.stderr.write("Usage: %s [ -no_mpi ] < input\n" % sys.argv[0])
    sys.stderr.write("input:\n")
    sys.stderr.write("max_volume_per_atom=float (1e3)\n")
    sys.stderr.write("min_volume_per_atom=float (1)\n")
    sys.stderr.write("start_species=int int [ float ] [, int int [ float ] ... ] (MANDATORY, this or start_config_file required. atomic_number multiplicity [not recomended: mass (amu)]. Info repeated for each species, separated by commas, mass is optional and not recommended.\n")
    sys.stderr.write("start_config_file=str (MANDATORY, this or start_species required. if set filename to read initial atom information from (instead of creating them)\n")
    sys.stderr.write("keep_atoms_fixed=int (0, no atoms to be fixed)\n")
    sys.stderr.write("apply_Z_wall=[T | F] (F, supported only for fortran MC)\n")
    sys.stderr.write("restart_file=path_to_file (file for restart configs. Mutually exclusive with start_*, one is required)\n")
    sys.stderr.write("n_walkers=int (MANDATORY)\n")
    sys.stderr.write("n_cull=int (1, number of walkers to kill at each NS iteration)\n")
    sys.stderr.write("n_extra_walk_per_task=int (0)\n")
    sys.stderr.write("n_iter_times_fraction_killed=int (MANDATORY, this or converge_down_to_T required)\n")
    sys.stderr.write("converge_down_to_T=float (MANDATORY, this or n_iter_times_fraction_killed required)\n")
    sys.stderr.write("T_estimate_finite_diff_lag=int (1000, lag for doing finite difference in current T estimate\n")
    sys.stderr.write("min_Emax=float (None.  Termination condition based on Emax)\n")
    sys.stderr.write("out_file_prefix=str (None)\n")
    sys.stderr.write("energy_calculator= ( ASE | lammps | internal | fortran) (fortran)\n")
    sys.stderr.write("n_extra_data=int (0, amount of extra data per atom to pass around)\n")
    sys.stderr.write("ns_run_analyzers=str ('', analyzers to apply during run\n")
    sys.stderr.write("\n")
    sys.stderr.write("KEmax_max_T=float (-1, if > 0 maximum temperature for estimating KEmax)\n")
    sys.stderr.write("kB=float (8.6173324e-5, Boltzmann constant)\n")
    sys.stderr.write("start_energy_ceiling_per_atom=float (1.0e9, max potential energy per atom for initial configs.  P*Vmax is added to this automatically)\n")
    sys.stderr.write("start_energy_ceiling=float (DEPRECATED: use start_energy_ceiling_per_atom. 1.0e9, max potential energy for initial configs.  P*Vmax is added to this automatically)\n")
    sys.stderr.write("random_init_max_n_tries=int (100, maximum number of tries for random initial positions under energy ceiling\n")
    sys.stderr.write("\n")
    sys.stderr.write("n_model_calls_expected=int (0, one of these is required)\n")
    sys.stderr.write("n_model_calls=int (0, one of these is required)\n")
    sys.stderr.write("do_good_load_balance=[T | F] (F, whether to do steps in groups with good load balance\n")
    sys.stderr.write("\n")
    sys.stderr.write("n_atom_steps=int (1, number of atomic trajectories \n")
    sys.stderr.write("atom_traj_len=int (8, length of atomic trajectory (MD steps or MC sweeps) in each step\n")
    sys.stderr.write("atom_traj_len_cost_multiplier=int (1, multiplier for cost of an atomic step (set to 1 for MD, TE-HMC, and SP-MC with O(1) cost moves, N for naive SP-MC)\n")
    sys.stderr.write("break_up_atom_traj=[T | F] (F, whether to intersperse n_atom_steps atomic sub-trajectories with other types of steps\n")
    sys.stderr.write("\n")
    sys.stderr.write("n_cell_volume_steps=int (1, number of cell MC volume steps )\n")
    sys.stderr.write("n_cell_shear_steps=int (1, number of cell MC shear steps )\n")
    sys.stderr.write("n_cell_stretch_steps=int (1, number of cell MC stretch steps )\n")
    sys.stderr.write("n_swap_steps=int (0, number of atom swaps)\n")
    sys.stderr.write("swap_max_cluster=int (1, maximum size of interconnected cluster to try to swap)\n")
    sys.stderr.write("swap_r_cut=float (2.5, cutoff radius for defining connected atoms for cluster)\n")
    sys.stderr.write("swap_cluster_probability_increment=float (0.75, factor between prob. of picking increasing larger clusters)\n")
    sys.stderr.write("swap_velo=[T | F] (F, if true, swap velocities when swapping atoms, breaking coherence a bit\n")
    sys.stderr.write("swap_velo_fix_mag_alt=[T | F] (F, if true, use alternate correction for velocity magnitudes when not swapping velocities\n")
    sys.stderr.write("\n")
    sys.stderr.write("n_semi_grand_steps=int (0, number of species type changes)\n")
    sys.stderr.write("semi_grand_potentials={ int : float [ , int : float ... ] } (None, chemical potentials for species type changes)\n")
    sys.stderr.write("\n")
    sys.stderr.write("velo_traj_len=int (0, number of MC sweeps in each velocity MC segement)\n")
    sys.stderr.write("\n")
    sys.stderr.write("random_energy_perturbation=float (1.0e-12)\n")
    sys.stderr.write("atom_algorithm=[MC | MD | GMC] (MANDATORY)\n")
    sys.stderr.write("\n")
    sys.stderr.write("MC_atom_velocities=[T | F] (F, supported only for energy_calculator=fortran)\n")
    sys.stderr.write("MC_atom_velocities_pre_perturb=[T | F] (F, Perturb velocities (rejection free) before MC + velocities walk)\n")
    sys.stderr.write("MC_atom_step_size=float (1.0, in units of (max_volume_per_atom * N_atoms)^(1/3) )\n")
    sys.stderr.write("MC_atom_step_size_max=float (1.0, in units of (max_volume_per_atom * N_atoms)^(1/3) )\n")
    sys.stderr.write("MC_atom_uniform_rv=[T | F] (F)\n")
    sys.stderr.write("MC_atom_Galilean=[T | F] (F)\n")
    sys.stderr.write("GMC_no_reverse=[T | F] (T)\n")
    sys.stderr.write("\n")
    sys.stderr.write("MD_atom_velo_pre_perturb=[T | F] (F. Perturb velocities before MD trajectory\n")
    sys.stderr.write("MD_atom_velo_post_perturb=[T | F] (T. Perturb velocities after MD trajectory\n")
    sys.stderr.write("MD_atom_velo_flip_accept=[T | F] (F)\n")
    sys.stderr.write("MD_atom_timestep=float (0.1, in ASE time units)\n")
    sys.stderr.write("MD_atom_timestep_max=float (2.0, in ASE time units)\n")
    sys.stderr.write("MD_atom_energy_fuzz=float (1.0e-2. Tolerance for rejecting non-energy conserving trajectories, as fraction of KE)\n")
    sys.stderr.write("MD_atom_reject_energy_violation=[ T | F ] (F, use energy conservation violation (exceeding MD_atom_energy_fuzz * KE) to reject MD trajectories)\n")
    sys.stderr.write("\n")
    sys.stderr.write("atom_velo_rej_free_fully_randomize=[T | F] (F. If true, randomize velocities completely rather than just perturbing.\n")
    sys.stderr.write("atom_velo_rej_free_perturb_angle=float (0.3. Max angle in radians for random rotations.)\n")
    sys.stderr.write("MC_atom_velo_step_size=float (50.0)\n")
    sys.stderr.write("MC_atom_velo_step_size_max=float (10000.0)\n")
    sys.stderr.write("MC_atom_velo_walk_rej_free=[T | F] (T. If true, use rejection free algorithm for MC_atom_walk\n")
    sys.stderr.write("\n")
    sys.stderr.write("\n")
    sys.stderr.write("MC_cell_P=float (0.0)\n")
    sys.stderr.write("MC_cell_flat_V_prior=[T | F] (F). POORLY TESTED, DO NOT TRUST (YET)\n")
    sys.stderr.write("MC_cell_volume_per_atom_step_size=float (5% of the maximum allowed volume)\n")
    sys.stderr.write("MC_cell_volume_per_atom_step_size_max=float (50% of the maximum allowed volume)\n")
    sys.stderr.write("MC_cell_volume_per_atom_prob=float (1.0)\n")
    sys.stderr.write("MC_cell_stretch_step_size=float (0.1)\n")
    sys.stderr.write("MC_cell_stretch_step_size_max=float (1.0)\n")
    sys.stderr.write("MC_cell_stretch_prob=float (1.0)\n")
    sys.stderr.write("MC_cell_shear_step_size=float (0.1 in units of (max_volume_per_atom * N_atoms)^(1/3))\n")
    sys.stderr.write("MC_cell_shear_step_size_max=float (1.0 in units of (max_volume_per_atom * N_atoms)^(1/3))\n")
    sys.stderr.write("MC_cell_shear_prob=float (1.0)\n")
    sys.stderr.write("MC_cell_min_aspect_ratio=float (0.8)\n")
    sys.stderr.write("cell_shape_equil_steps=int (1000)\n")
    sys.stderr.write("\n")
    sys.stderr.write("monitor_step_interval_times_fraction_killed=float (1, divided by n_cull/n_walkers to get actual monitoring interval in iterations, negative for only using last iteration, 0 for no monitoring)\n")
    sys.stderr.write("adjust_step_interval_times_fraction_killed=float (5, divided by n_cull/n_walkers to get actual adjustment interval in iterations, negative for only using last iteration, 0 for no adjust)\n")
    sys.stderr.write("full_auto_step_sizes=[T | F] (T) (T. Automatically calibrate all sizes by performing additional short explorations, including at start of run. F. Use initial input step sizes and make small adjustments to step sizes during run.)\n")
    sys.stderr.write("MC_adjust_step_factor=float (1.1)\n")
    sys.stderr.write("MC_adjust_min_rate=float (0.25)\n")
    sys.stderr.write("MC_adjust_max_rate=float (0.75)\n")
    sys.stderr.write("GMC_adjust_min_rate=float (0.25)\n")
    sys.stderr.write("GMC_adjust_max_rate=float (0.75)\n")
    sys.stderr.write("GMC_dir_perturb_angle=float (-1.0)\n")
    sys.stderr.write("GMC_dir_perturb_angle_during=float (0.0)\n")
    sys.stderr.write("MD_adjust_step_factor=float (1.1)\n")
    sys.stderr.write("MD_adjust_min_rate=float (0.5)\n")
    sys.stderr.write("MD_adjust_max_rate=float (0.95)\n")
    sys.stderr.write("\n")
    sys.stderr.write("ASE_calc_module=str (MANDATORY if energy_calculator=ASE)\n")
    sys.stderr.write("FORTRAN_model=str (MANDATORY if energy_calculator=fortran)\n")
    sys.stderr.write("FORTRAN_model_params=str (parameters for energy_calculator=fortran)\n")
    sys.stderr.write("LAMMPS_fix_gmc=[T | F]\n")
    sys.stderr.write("LAMMPS_init_cmds=str (MANDATORY if energy_calculator=lammps)\n")
    sys.stderr.write("LAMMPS_name=str ('', arch name for lammps shared object file)\n")
    sys.stderr.write("LAMMPS_serial=[T | F] (T, is lammps library serial)\n")
    sys.stderr.write("LAMMPS_header=str (lammpslib.py value default, override lammpslib.py header commands for energy_calculator=lammps)\n")
    sys.stderr.write("LAMMPS_header_extra=str ('', extra lammpslib.py header commands for energy_calculator=lammps)\n")
    sys.stderr.write("LAMMPS_atom_types=symbol int [, symbol int ] ... ('', mapping from atomic symbols to type numbers for LAMMPS ASE interface)\n")
    sys.stderr.write("\n")
    sys.stderr.write("config_file_format=str (extxyz)\n") # julia
    sys.stderr.write("rng=( numpy | internal | rngstream ) (numpy)\n") # julia
    sys.stderr.write("profile=rank_to_profile (-1)\n")
    sys.stderr.write("2D=[ T | F ] (F)\n")
    sys.stderr.write("Z_cell_axis=float (10.0)\n")
    sys.stderr.write("debug=debug_level (0, <= 0 for no debugging tests/prints)\n")
    sys.stderr.write("snapshot_interval=int (10000, <=0 for no snapshots except final positions)\n")
    sys.stderr.write("snapshot_time=int (3600, <=0 for no time based snapshots)\n")
    sys.stderr.write("snapshot_per_parallel_task=[T | F] (T, if true save separate snapshot file for each parallel task)\n")
    sys.stderr.write("snapshot_clean=[T | F] (T, if true clean previous iter snapshots\n")
    sys.stderr.write("random_initialise_pos=[T | F] (T, if true randomize the initial positions\n")
    sys.stderr.write("random_initialise_cell=[T | F] (T, if true randomize the initial cell\n")
    sys.stderr.write("LAMMPS_molecular_info=[T | F] (T, if true create lammps molecular info\n")
    sys.stderr.write("initial_walk_N_walks=int (0 number of rounds for initial walk) \n")
    sys.stderr.write("initial_walk_adjust_interval=int (10 interval (in walks) between adjustments of steps during initial walk) \n")
    sys.stderr.write("initial_walk_Emax_offset_per_atom=float (1.0, offset (per atom) to add to Emax for initial walks) \n")
    sys.stderr.write("initial_walk_only=[T | F] (F, quit after doing initial walk) \n")
    sys.stderr.write("traj_interval=int (100, <=0 for no trajectory)\n")
    sys.stderr.write("delta_random_seed=int (-1, < 0 for seed from /dev/urandom)\n")
    sys.stderr.write("no_extra_walks_at_all=[ T | F ] (F)\n")
    sys.stderr.write("track_configs=[ T | F ] (F)\n")
    sys.stderr.write("track_configs_write=[ T | F ] (F)\n")

def excepthook_mpi_abort(exctype, value, tb):
    print( print_prefix,'Uncaught Exception Type:', exctype)
    print( print_prefix,'Value:', value)
    print( print_prefix,'Traceback:', tb)
    # print_tb(tb)
    print_exception(exctype, value, tb)
    print( print_prefix, "Aborting")
    sys.stdout.flush()
    try:
        MPI.COMM_WORLD.Abort(1)
    except:
        pass
    sys.exit(1)

def exit_error(message, stat):
    sys.stderr.write(message)
    sys.stderr.flush()
    try:
        MPI.COMM_WORLD.Abort(stat)
    except:
        pass
    sys.exit(stat)

def str_to_logical(x):
    if x in [ "t", "T", "true", "True", "TRUE", "1" ]:
        return True
    if x in [ "f", "F", "false", "False", "FALSE", "0" ]:
        return False
    exit_error("Failed to parse logical\n",5)

def internal_Vlj(dr_vec, Eshift=0.0):
    return (np.sum(dr_vec**-12 - dr_vec**-6) - Eshift*len(dr_vec))

# pos is Nx3, l is a scalar box size
def energy_internal_pos(pos, l):
    n_at=np.size(pos,0)
    dr_list = np.zeros( (n_at*(n_at-1)/2,3) )
    offset=0
    # calc relative pos vector
    for i in range(n_at):
        dr_list[offset:offset+(n_at-i-1),0] = pos[i,0] - pos[i+1:,0]
        dr_list[offset:offset+(n_at-i-1),1] = pos[i,1] - pos[i+1:,1]
        dr_list[offset:offset+(n_at-i-1),2] = pos[i,2] - pos[i+1:,2]
        offset += n_at-i-1
    # apply min image
    dr_list[:,:] -= l*np.floor(dr_list[:,:]/l+0.5)
    # magnitude of vectors
    dr_list_mag = np.sqrt(np.sum(dr_list*dr_list,axis=1))

    return np.sum(internal_Vlj(dr_list_mag[np.where(dr_list_mag < internal_cutoff)],Eshift))

def energy_internal(at):
    return energy_internal_pos(at.get_positions(), at.get_cell()[0,0])

def eval_energy_PE(at):
    if do_calc_ASE or do_calc_lammps:
        if do_calc_lammps:
            #NB only MD can make crazy positions, so maybe just do this after MD propagation?
            at.wrap()
        energy = at.get_potential_energy()
    elif do_calc_internal:
        energy = energy_internal(at)
    elif do_calc_fortran:
        energy = f_MC_MD.eval_energy(at)
    else:
        sys.stderr.write("No way to eval_energy_PE()\n", 5)

    return energy

def eval_energy_KE(at):
    if at.has('momenta') and at.has('masses'):
        energy = at.get_kinetic_energy()
    else:
        energy = 0.0

    return energy

def eval_energy_PV(at):
    return movement_args['MC_cell_P']*at.get_volume()

def eval_energy_mu(at):
    if movement_args['semi_grand_potentials'] is not None:
        mu = movement_args['semi_grand_potentials']
        energy = np.sum([ mu[Z] for Z in at.get_atomic_numbers() ])
    else:
        energy = 0.0

    return energy

def eval_energy(at, do_KE=True, do_PE=True):
    """Calls appropriate functions to calculate the potential energy, kinetic energy and the p*V term.
    """
    energy = eval_energy_PV(at) + eval_energy_mu(at)

    if do_PE:
        energy += eval_energy_PE(at)

    if do_KE:
        energy += eval_energy_KE(at)

    return energy

def eval_forces(at):
    if do_calc_ASE or do_calc_lammps:
        if do_calc_lammps:
            #NB only MD can make crazy positions, so maybe just do this after MD propagation?
            at.wrap()
        forces = at.get_forces()
    elif do_calc_internal:
        exit_error('no forces for do_calc_internal', 10)
    elif do_calc_fortran:
        forces = np.zeros( at.get_positions().shape )
        f_MC_MD.eval_forces(at, forces)
    else:
        sys.stderr.write("No way to eval_forces()\n", 5)

    return forces

def propagate_NVE_ASE(at, dt, n_steps):
    at.calc = pot

    # dt is being converted from ASE units to fs
    if ns_args['debug'] >= 10:
        if movement_args['keep_atoms_fixed'] > 0: #LIVIA
            from ase.constraints import FixAtoms
            c = FixAtoms(indices=[atom.index for atom in at if atom.index < movement_args['keep_atoms_fixed']])
            at.set_constraint(c)
        vv = VelocityVerlet(at, timestep=dt, logfile='-')

    else:
        if movement_args['keep_atoms_fixed'] > 0: #LIVIA
            from ase.constraints import FixAtoms
            c = FixAtoms(indices=[atom.index for atom in at if atom.index < movement_args['keep_atoms_fixed']])
            at.set_constraint(c)
        vv = VelocityVerlet(at, timestep=dt)

    vv.run(n_steps)


def propagate_lammps(at, dt, n_steps, algo, Emax=None):
    if pot.first_propagate:
        pot.first_propagate = False
    else:
        pot.lmp.command('unfix 1 ')

    # set appropriate fix
    if algo == 'NVE':
        if movement_args['keep_atoms_fixed'] > 0: #LIVIA
            #pot.lmp.command('group mobile id > ' + str(movement_args['keep_atoms_fixed']))
            pot.lmp.command('fix 1 mobile nve')
        else:
            pot.lmp.command('fix 1 all nve')
    elif algo == 'GMC':
        # GMC does not work with fixed atoms # Ray Yang
        if movement_args['keep_atoms_fixed'] > 0:
            exit_error("propagate_lammps got algo '%s', which does not work with fixed atoms\n" % algo)
        # Hard coded value needed for LAMMPS. Let's hope our RNG maximum is at
        # least this large.
        pot.lmp.command('fix 1 all ns/gmc {} {} '.format(rng.int_uniform(1, 900000000), Emax))
    else:
        exit_error(
            "propagate_lammps got algo '%s', neither NVE nor GMC\n" % algo)

    # NB maybe not do this every time? Just _after_ MD, since that's the only
    # way position can become crazy?
    at.wrap()

    # do propagation, recover from exception if needed
    try:
        if algo == 'NVE':
            pot.propagate(at, properties=['energy', 'forces'],
                          system_changes=['positions'], n_steps=n_steps, dt=dt)
        else:
            pot.propagate(at, properties=['energy', 'forces'],
                          system_changes=['positions'], n_steps=n_steps, dt=dt,
                          dt_not_real_time=True)
    except Exception as err:
        # clean up and return failure
        if ns_args['debug'] >= 4:
            print("propagate_lammps got exception ", err)
        pot.restart_lammps(at)
        pot.first_propagate = True
        return False

    return True


def velo_rv_mag(n):
    if movement_args['2D']:
        #unit_rv[:,2] = 0.0 # commented out from here, as unit_rv is unknown to this function at this point
        nDOF = 2.0
    else:
        nDOF = 3.0
    # In 3D rv_mag should have prob distrib p(r) = r^(3N-1).
    # Using transformation rule p(y) = p(x) |dx/dy|, with p(y) = y^{3N-1} and p(x) = 1,
    #       one gets dx/dy = y^{3N-1}
    #                x = y^{3N}
    #                y = x^{1/3N}
    return rng.float_uniform(0.0, 1.0)**(1.0/(nDOF*n))


def velo_unit_rv(n):
    unit_rv = rng.normal(1.0, (n, 3))
    unit_rv /= np.linalg.norm(unit_rv)
    return unit_rv


def gen_random_velo(at, KEmax, unit_rv=None):
    if unit_rv is None:
        unit_rv = velo_unit_rv(len(at))
    if movement_args['2D']:
        unit_rv[:, 2] = 0.0

    # propose right magnitude of randomness consistent with mobile atoms
    if movement_args['keep_atoms_fixed'] > 0:
        rv_mag = velo_rv_mag(len(at)-movement_args['keep_atoms_fixed'])
    else:
        rv_mag = velo_rv_mag(len(at))
    #rv_mag = velo_rv_mag(len(at))

    # from Baldock thesis Eq. 11.10
    #     p^{**} = r \mathbf{S} \hat{\mathbf{r}}
    # and 11.11
    #     S_{ij} = \delta_{ij} (2 m_i [ E_{lim} - U(q))])^{1/2}
    # p_i = r (2 m_i)^{1/2} (E-U)^{1/2} \hat{r}_i
    # v_i = p_i / m_i = r (2/m)^{1/2} (E-U)^{1/2} \hat{r}_i

    masses = at.get_masses()

    # the "np.sqrt(2.0/np.array([masses,]*3).transpose()) * np.sqrt(KEmax)" is filling it up with the same term
    velocities = rv_mag * np.sqrt(2.0/np.array([masses,]*3).transpose()) * np.sqrt(KEmax) * unit_rv
    # zero velocities of fixed atoms (these were not taken into account in the nDOF)
    if movement_args['keep_atoms_fixed'] > 0:
        velocities[:movement_args['keep_atoms_fixed']]=0.0

    return velocities


def pairwise(iterable):
    a = iter(iterable)
    return zip(a, a)


def rotate_dir_3N(vec, max_ang):
    if max_ang <= 0.0:
        return
    # apply random rotations
    indices = [ (int(i/3), i%3) for i in range(3*vec.shape[0]) ]
    rng.shuffle_in_place(indices)
    for ((ind_1_a,ind_1_c), (ind_2_a,ind_2_c)) in pairwise(indices):
        ang = rng.float_uniform(-max_ang,max_ang)
        c_ang = np.cos(ang)
        s_ang = np.sin(ang)
        v_1 =  vec[ind_1_a,ind_1_c] * c_ang + vec[ind_2_a,ind_2_c] * s_ang
        v_2 = -vec[ind_1_a,ind_1_c] * s_ang + vec[ind_2_a,ind_2_c] * c_ang
        vec[ind_1_a,ind_1_c] = v_1
        vec[ind_2_a,ind_2_c] = v_2

def rotate_dir_3N_2D(vec, max_ang):
    if max_ang <= 0.0:
        return
    # apply random rotations
    indices = [ (int(i/2), i%2) for i in range(2*vec.shape[0]) ]
    rng.shuffle_in_place(indices)
    for ((ind_1_a,ind_1_c), (ind_2_a,ind_2_c)) in pairwise(indices):
        ang = rng.float_uniform(-max_ang,max_ang)
        c_ang = np.cos(ang)
        s_ang = np.sin(ang)
        v_1 =  vec[ind_1_a,ind_1_c] * c_ang + vec[ind_2_a,ind_2_c] * s_ang
        v_2 = -vec[ind_1_a,ind_1_c] * s_ang + vec[ind_2_a,ind_2_c] * c_ang
        vec[ind_1_a,ind_1_c] = v_1
        vec[ind_2_a,ind_2_c] = v_2


def rej_free_perturb_velo(at, Emax, KEmax, rotate=True):
#DOC
#DOC ``rej_free_perturb_velo``

    if not at.has('momenta') or not at.has('masses'):
        return

    KEmax_use = None
    if Emax is not None:
        initial_KE = eval_energy_KE(at)
        KEmax_use = Emax - (at.info['ns_energy'] - initial_KE)
    if KEmax > 0.0:
        if Emax is None:
            KEmax_use = KEmax
        else:
            KEmax_use = min(KEmax_use, KEmax)
    if KEmax_use is None:
        exit_error("rej_free_perturb_velo() called with Emax and KEmax both None\n", 9)

    orig_KE = at.get_kinetic_energy()
    #DOC \item if atom\_velo\_rej\_free\_fully\_randomize, pick random velocities consistent with Emax
    if movement_args['atom_velo_rej_free_fully_randomize']:
        # randomize completely, fixed atoms are taken care of
        at.set_velocities(gen_random_velo(at, KEmax_use))  

    #DOC \item else perturb velocities
    else:  # perturb
        velo = at.get_velocities()

        velo_mag = np.linalg.norm(velo)
        #DOC \item if current velocity=0, can't rescale, so pick random velocities consistent with Emax
        if velo_mag == 0.0:
            # generate random velocities, fixed atoms are taken care of
            at.set_velocities(gen_random_velo(at, KEmax_use))  

        #DOC \item else, pick new random magnitude consistent with Emax, random rotation of current direction with angle uniform in +/- atom\_velo\_rej\_free\_perturb\_angle
        else:
            # using default settings we will do this. 
            # pick new random magnitude - count on dimensionality to make change small
            # WARNING: check this for variable masses

            sqrt_masses_2D = np.sqrt(at.get_masses().reshape((len(at), 1)))
            scaled_vel = gen_random_velo(at, KEmax_use, velo/velo_mag) * sqrt_masses_2D

            if rotate:
                if movement_args['2D']:
                    rotate_dir_3N_2D(scaled_vel, movement_args['atom_velo_rej_free_perturb_angle'])
                else:
                    if movement_args['keep_atoms_fixed'] > 0:  # this is the new line
                        rotate_dir_3N(scaled_vel[movement_args['keep_atoms_fixed']:, :], movement_args['atom_velo_rej_free_perturb_angle'])
                    else:
                        rotate_dir_3N(scaled_vel, movement_args['atom_velo_rej_free_perturb_angle'])

            at.set_velocities(scaled_vel / sqrt_masses_2D)
            if movement_args['keep_atoms_fixed'] > 0:  # for surface simulations
                random_velo = at.get_velocities()
                random_velo[:movement_args['keep_atoms_fixed'], :] = 0
                at.set_velocities(random_velo)

    new_KE = at.get_kinetic_energy()

    # rej_free_perturb_velo expects at.info['ns_energy'] to be set correctly initially
    at.info['ns_energy'] += new_KE-orig_KE


def do_MC_atom_velo_walk(at, movement_args, Emax, nD, KEmax):
#DOC
#DOC ``do\_MC\_atom\_velo\_walk``

    #DOC \item Else if MC\_atom\_velo\_walk\_rej\_free
    if movement_args['MC_atom_velo_walk_rej_free']:
        #DOC \item call rej\_free\_perturb\_velo()
            rej_free_perturb_velo(at, Emax, KEmax)
            return {}
    #DOC \item else do MC pertubation to velocities

    n_steps = movement_args['velo_traj_len']
    step_size = movement_args['MC_atom_velo_step_size']

    n_try = n_steps*len(at)

    initial_KE = eval_energy_KE(at)
    KEmax_use = Emax - (at.info['ns_energy'] - initial_KE) # expecting ns_energy = KE + PE (+PV)
    if KEmax > 0.0 and KEmax < KEmax_use:
        KEmax_use = KEmax

    if do_calc_fortran:
        (n_try, n_accept, final_KE) = f_MC_MD.MC_atom_walk_velo(at, n_steps, step_size, nD, KEmax_use)
        at.info['ns_energy'] += final_KE-initial_KE
    else:
        masses = at.get_masses()
        velocities = at.get_velocities()
        KE = initial_KE
        n_accept = 0
        for i_step in range(n_steps):
            at_list = list(range(len(at)))
            rng.shuffle_in_place(at_list)
            d_vel = rng.normal(step_size, (len(at), 3))
            if movement_args['keep_atoms_fixed'] > 0:  # for surface simulations
                d_vel[:movement_args['keep_atoms_fixed'], :] = 0
            for i_at in at_list:
                d_KE = 0.5*masses[i_at]*(np.sum((velocities[i_at,:]+d_vel[i_at,:])**2) - np.sum(velocities[i_at,:]**2))
                if KE + d_KE < KEmax_use:
                    velocities[i_at,:] += d_vel[i_at,:]
                    KE += d_KE
                    n_accept += 1
        at.set_velocities(velocities)
        at.info['ns_energy'] += KE-initial_KE

    return {'MC_atom_velo': (n_try, n_accept)}


def do_MD_atom_walk(at, movement_args, Emax, KEmax):
#DOC
#DOC ``do_MD_atom_walk``

    """ perform MD walk on the configuration """
    orig_E = at.info['ns_energy']
    nD = 3
    if movement_args['2D']:
       nD = 2

    if orig_E >= Emax:
        print(print_prefix, ": WARNING: orig_E =", orig_E, " >= Emax =", Emax)

    #DOC \item if MD\_atom\_velo\_pre\_perturb, call do\_MC\_atom\_velo\_walk() for magnitude and rotation
    if movement_args['MD_atom_velo_pre_perturb']:
        do_MC_atom_velo_walk(at, movement_args, Emax, nD, KEmax)

    pre_MD_pos = at.get_positions()
    pre_MD_velo = at.get_velocities()
    if ns_args['n_extra_data'] > 0:
        pre_MD_extra_data = at.arrays['ns_extra_data'].copy()

    pre_MD_E = at.info['ns_energy']

    #DOC \item propagate in time atom\_traj\_len time steps of length MD\_atom\_timestep
    if movement_args['python_MD']:
        forces = eval_forces(at)
        final_E = None
        timestep = movement_args['MD_atom_timestep']
        for i in range(movement_args['atom_traj_len']):
            at.set_momenta(at.get_momenta()+forces*0.5*timestep)
            at.set_positions(at.get_positions()+at.get_momenta()*timestep)
            try:
                forces = eval_forces(at)
            except:
                final_E = 2.0*abs(Emax)
                break
            at.set_momenta(at.get_momenta()+forces*0.5*timestep)
        if final_E is None:  # didn't abort due to exception in eval_forces()
            final_E = eval_energy(at)
    else:
        if do_calc_ASE:
            propagate_NVE_ASE(at, dt=movement_args['MD_atom_timestep'], n_steps=movement_args['atom_traj_len'])
            final_E = eval_energy(at)
        elif do_calc_lammps:

            if propagate_lammps(at, dt=movement_args['MD_atom_timestep'], n_steps=movement_args['atom_traj_len'], algo='NVE'):
                final_E = pot.results['energy'] + eval_energy(at, do_PE=False)
            else:  # propagate returned success == False
                final_E = 2.0*abs(Emax)
                ## print("error in propagate_lammps NVE, setting final_E = 2*abs(Emax) =" , final_E)

        elif do_calc_fortran:
            # Fortran MD code does not support fixed atoms # Ray Yang
            if movement_args['keep_atoms_fixed'] > 0:
                exit_error("do_MD_atom_walk() called with keep_atoms_fixed > 0, but no way to do that with fortran code\n", 3)
            final_E = f_MC_MD.MD_atom_NVE_walk(
                at, n_steps=movement_args['atom_traj_len'],
                timestep=movement_args['MD_atom_timestep'],
                debug=ns_args['debug'])
            final_E += eval_energy(at, do_PE=False, do_KE=False)
        else:
            exit_error("Need some non-ASE, non-fortran, non-lammps way of doing MD\n",3)

    reject_fuzz = False
    final_KE = eval_energy_KE(at)
    #DOC \item If MD\_atom\_reject\_energy\_violation is set, accept/reject entire move on E deviating by less than MD\_atom\_energy\_fuzz times kinetic energy
    if abs(final_E-pre_MD_E) > movement_args['MD_atom_energy_fuzz'] * final_KE:
        if movement_args['MD_atom_reject_energy_violation']:
            reject_fuzz = True
        # else:
            # print(print_prefix, ": WARNING: MD energy deviation > fuzz*final_KE. Pre-MD, post-MD, difference, final_KE ", pre_MD_E, final_E, final_E-pre_MD_E, final_KE)

    #DOC \item accept/reject entire move on E < Emax and KE < KEmax
    reject_Emax = (final_E >= Emax)
    reject_KEmax = (KEmax > 0.0 and final_KE >= KEmax)

    #DOC \item if reject
    if reject_fuzz or reject_Emax or reject_KEmax:  # reject
        #DOC \item set positions, velocities, energy back to value before perturbation (maybe should be after?)
        # print(print_prefix, ": WARNING: reject MD traj Emax ", Emax, " initial E ", orig_E, " velo perturbed E ", pre_MD_E, " final E ",final_E, " KEmax ", KEmax, " KE ", final_KE)
        at.set_positions(pre_MD_pos)
        if movement_args['MD_atom_velo_flip_accept']:
            at.set_velocities(pre_MD_velo)  # TODO: should we be redundant in zeroing?
        else:
            at.set_velocities(-pre_MD_velo)
        if ns_args['n_extra_data'] > 0:
            at.arrays['ns_extra_data'][...] = pre_MD_extra_data
        at.info['ns_energy'] = pre_MD_E
        n_accept = 0
    #DOC \item else
    else:  # accept
        #DOC \item flip velocities if MD\_atom\_velo\_flip\_accept
        # remember to reverse velocities on acceptance to preserve detailed
        # balance, since velocity is (sometimes) being perturbed, not completely
        # randomized
        if movement_args['MD_atom_velo_flip_accept']:
            # TODO: negative of 0 is 0
            at.set_velocities(-at.get_velocities())  # is there a faster way of doing this in ASE?  Can you do at.velocities?
        at.info['ns_energy'] = final_E
        n_accept = 1

    #DOC \item if MD\_atom\_velo\_post\_perturb, call do\_MC\_atom\_velo\_walk() for magnitude and rotation
    if movement_args['MD_atom_velo_post_perturb']:
        do_MC_atom_velo_walk(at, movement_args, Emax, nD, KEmax)

    return {'MD_atom': (1, n_accept)}


def do_MC_atom_walk(at, movement_args, Emax, KEmax):
#DOC
#DOC ``do_MC_atom_walk``

    n_steps = movement_args['atom_traj_len']
    step_size = movement_args['MC_atom_step_size']
    step_size_velo = movement_args['MC_atom_velo_step_size']
    n_try = n_steps*(len(at)-movement_args['keep_atoms_fixed']) # does this do anything? actual n_try is set within fortran code 
    n_accept=0
    n_accept_velo = None

    nD=3
    if movement_args['2D']:
       nD=2

    #DOC \item if MC\_atom\_velocities and MC\_atom\_velocities\_pre\_perturb, call do\_MC\_atom\_velo\_walk() to perturb velocities, magnitude and and rotation
    if movement_args['MC_atom_velocities'] and movement_args['MC_atom_velocities_pre_perturb']:
        do_MC_atom_velo_walk(at, movement_args, Emax, nD, KEmax)

    if movement_args['MC_atom_Galilean']:
        if movement_args['GMC_dir_perturb_angle'] < 0.0 or np.linalg.norm(at.arrays['GMC_direction']) == 0.0:
            # completely random orientation, magnitude 1
            at.arrays['GMC_direction'][:,:] = rng.normal(1.0, (len(at), 3))
            if (movement_args['keep_atoms_fixed'] > 0):
                at.arrays['GMC_direction'][:movement_args['keep_atoms_fixed'],:] = 0.0
            at.arrays['GMC_direction'] /= np.linalg.norm(at.arrays['GMC_direction']) # what is this line doing?

        elif movement_args['GMC_dir_perturb_angle'] > 0.0:
            # apply random rotations
            if movement_args['keep_atoms_fixed'] > 0:
                rotate_dir_3N(at.arrays['GMC_direction'][movement_args['keep_atoms_fixed']:,:], movement_args['GMC_dir_perturb_angle'])

            else:
                rotate_dir_3N(at.arrays['GMC_direction'], movement_args['GMC_dir_perturb_angle'])

    #DOC \item if using fortran calculator and not reproducible
    if do_calc_fortran and not ns_args['reproducible']:
        #DOC \item call fortran MC code f\_MC\_MD.MC\_atom\_walk
        at.wrap()
        if (ns_args['debug'] >= 10 and movement_args['2D']): 
            pp=at.get_positions()
            if (any(pp[:,2])>0.00001 and nD==2): #LIVIA
                exit_error("Not a 2D system anymore\n",3)
 
        if movement_args['MC_atom_velocities']:
            (n_try, n_accept, n_accept_velo, final_E) = f_MC_MD.MC_atom_walk(at, n_steps, step_size, Emax-eval_energy(at, do_PE=False, do_KE=False), nD, movement_args['keep_atoms_fixed'], KEmax, step_size_velo)
            at.info['ns_energy'] = final_E + eval_energy(at, do_PE=False, do_KE=False)
        else:
            if movement_args['MC_atom_Galilean']:
                (n_try, n_accept, final_E) = f_MC_MD.GMC_atom_walk(at, n_steps, step_size, Emax-eval_energy(at, do_PE=False), no_reverse=movement_args['GMC_no_reverse'], pert_ang=movement_args['GMC_dir_perturb_angle_during'])
            else:
                # fixed atoms MC works by including the number of atoms to be kept fixed - LIVIA
                (n_try, n_accept, final_E) = f_MC_MD.MC_atom_walk(at, n_steps, step_size, Emax-eval_energy(at, do_PE=False), nD, movement_args['keep_atoms_fixed'], movement_args['wall_dist'])
            at.info['ns_energy'] = final_E + eval_energy(at, do_PE=False, do_KE=True)

    elif (do_calc_lammps and movement_args['MC_atom_Galilean'] and ns_args['LAMMPS_fix_gmc']):
        orig_pos = at.get_positions()
        if (ns_args['debug'] >= 10 and movement_args['2D']): 
            if (any(orig_pos[:,2])>0.00001): #LIVIA
                exit_error("Not a 2D system anymore\n",3)
        orig_energy = at.info['ns_energy']
        extra_term = eval_energy(at, do_PE=False, do_KE=False)
        if propagate_lammps(at, step_size, n_steps, algo='GMC', Emax=Emax-extra_term ):
            energy1 = pot.results['energy'] + extra_term
            if energy1 < Emax:
                at.info['ns_energy'] = energy1
                n_accept = n_steps*len(at)
            else:
                at.info['ns_energy'] = orig_energy
                at.set_positions(orig_pos)
                n_accept = 0
        else: # got an exception, reject traj
            at.info['ns_energy'] = orig_energy
            at.set_positions(orig_pos)
            n_accept = 0

    #DOC \item else
    else:
        #DOC \item do python MC
        if movement_args['MC_atom_velocities']:
            exit_error("MC_atom_velocities only supported for FORTRAN calculator\n", 8)
        dz=0.0
        #DOC \item if MC_atom_Galilean
        if movement_args['MC_atom_Galilean']:
            #DOC \item go Galilean MC in python

            do_no_reverse = movement_args['GMC_no_reverse']

            if do_no_reverse:
                last_good_pos = at.get_positions()
            else:
                n_reverse = 0

            n_reflect = 0
            pos = at.get_positions()
            d_pos = step_size*at.arrays['GMC_direction']

            for i_MC_step in range(n_steps):
                if not do_no_reverse:
                    last_good_pos = at.get_positions()
                    last_good_d_pos = d_pos.copy()

                # perturb direction if needed
                rotate_dir_3N(d_pos, movement_args['GMC_dir_perturb_angle_during'])

                # step and evaluate
                pos += d_pos
                at.set_positions(pos)
                E = eval_energy(at)
                cur_E_is_correct = True

                if E >= Emax: # reflect
                    Fhat = eval_forces(at)
                    Fhat /= np.linalg.norm(Fhat)
                    d_pos -= 2.0*Fhat*np.sum(Fhat*d_pos)

                    n_reflect += 1

                    if not do_no_reverse: # do reflection step and check for reverse
                        # step on reflected traj and evaluate
                        pos += d_pos
                        at.set_positions(pos)
                        E = eval_energy(at)
                        cur_E_is_correct = True

                        if E >= Emax: # reverse
                            pos[:,:] = last_good_pos
                            at.set_positions(pos)
                            cur_E_is_correct = False
                            d_pos[:,:] = -last_good_d_pos
                            n_reflect -= 1
                            n_reverse += 1

            if do_no_reverse: # accept/reject
                n_try = 1
                if E < Emax: # accept
                    # save new E, pos, direction
                    at.info['ns_energy'] = E
                    at.set_positions(pos)
                    at.arrays['GMC_direction'][:,:] = d_pos/np.linalg.norm(d_pos)
                    n_accept = 1
                else: # reject
                    # E and GMC_direction in at were never overwritten, no need to restore, but do need to reverse dir
                    at.set_positions(last_good_pos)
                    at.arrays['GMC_direction'] *= -1.0
                    n_accept = 0
            else:
                if n_reverse > 0 and not cur_E_is_correct:
                    E = eval_energy(at)
                at.info['ns_energy'] = E
                at.set_positions(pos)
                at.arrays['GMC_direction'][:,:] = d_pos/np.linalg.norm(d_pos)
                n_try = n_reflect + n_reverse
                n_accept = n_reflect

        #DOC \item else
        else:
            #DOC \item loop atom\_traj\_len times
            for i_MC_step in range(n_steps):
                #DOC \item loop over atoms in random order
                at_list=list(range(len(at)))
                rng.shuffle_in_place(at_list)
                for i_at in at_list:
                    #DOC \item propose single atom move
                    if movement_args['MC_atom_uniform_rv']:
                        dx = rng.float_uniform(-step_size,step_size)
                        dy = rng.float_uniform(-step_size,step_size)
                        dz = 0.0
                        if not movement_args['2D']:
                            dz = rng.float_uniform(-step_size,step_size)
                    else:
                        dx = rng.normal(step_size)
                        dy = rng.normal(step_size)
                        dz = 0.0
                        if not movement_args['2D']:
                            dz = rng.normal(step_size)
                    orig_energy = at.info['ns_energy']
                    orig_pos = at.get_positions()
                    if ns_args['n_extra_data'] > 0:
                        orig_extra_data = at.arrays['ns_extra_data'].copy()
                    new_pos = orig_pos.copy()
                    new_pos[i_at,:] += (dx, dy,dz)
                    at.set_positions(new_pos)
                    #DOC \item accept/reject on E < Emax
                    energy = eval_energy(at)
                    if energy >= Emax:
                        # reject move
                        at.set_positions(orig_pos)
                        if ns_args['n_extra_data'] > 0:
                            at.arrays['ns_extra_data'][...] = orig_extra_data
                        energy = orig_energy
                    else:
                        n_accept += 1
                    at.info['ns_energy'] = energy

    out = {}
    if n_accept_velo is not None:
        out['MC_atom_velo'] = (n_try, n_accept_velo)
    out['MC_atom'] = (n_try, n_accept)

    return out

def propose_volume_step(at, step_size, flat_V_prior):
    dV = rng.normal(step_size*len(at))
    orig_V = at.get_volume()
    new_V = orig_V+dV
    if new_V < 0: # negative number cannot be raised to fractional power, so this is only to avoid fatal error during the run
        new_V=abs(new_V)
        # print("Warning, the step_size for volume change might be too big, resulted in negative new volume", step_size, dV, orig_V+dV)
    #print("TRANSFORM", new_V, orig_V, dV, step_size, len(at))
    transform = np.identity(3)*(new_V/orig_V)**(1.0/3.0)
    if flat_V_prior:
        p_accept = 1.0
    else:
        p_accept = min(1.0, (new_V/orig_V)**len(at))
    return (p_accept, transform)

def propose_area_step(at, step_size, flat_V_prior):
    dA = rng.normal(step_size*len(at))
    orig_cell = at.get_cell()
    orig_A = at.get_volume() / orig_cell[2,2]
    new_A = orig_A+dA
    if new_A < 0: # negative number cannot be raised to fractional power, so this is only to avoid fatal error during the run
        new_A=abs(new_A)
        # print("Warning, the step_size for volume change might be too big, resulted in negative new volume", step_size, dV, orig_V+dV)
    #print("TRANSFORM", new_V, orig_V, dV, step_size, len(at))
    transform = np.identity(3)*np.sqrt(new_A/orig_A)
    transform[2,2] = 1.0 # no change in the Z direction
    if flat_V_prior:
        p_accept = 1.0
    else:
        p_accept = min(1.0, (new_A/orig_A)**len(at))
    return (p_accept, transform)

def propose_shear_step(at, step_size):
   if movement_args['2D']:
       # pick random vector
       rnd_vec_ind = rng.int_uniform(0, 2)
       # turn other two into orthonormal pair
       other_vec_ind = list(range(3))
       other_vec_ind.remove(rnd_vec_ind)
       orig_cell = at.get_cell()
       v1 = orig_cell[other_vec_ind[0],:].copy()
       v1 /= np.sqrt(np.dot(v1,v1))
       #pick random magnitude
       rv1 = rng.normal(step_size)
       # create new cell and transformation matrix
       new_cell = orig_cell.copy()
       new_cell[rnd_vec_ind,:] += rv1*v1

   else:
       # pick random vector
       rnd_vec_ind = rng.int_uniform(0, 3)
       # turn other two into orthonormal pair
       other_vec_ind = list(range(3))
       other_vec_ind.remove(rnd_vec_ind)
       orig_cell = at.get_cell()
       v1 = orig_cell[other_vec_ind[0],:].copy()
       v2 = orig_cell[other_vec_ind[1],:].copy()
       v1 /= np.sqrt(np.dot(v1,v1))
       v2 -= v1*np.dot(v1,v2)
       v2 /= np.sqrt(np.dot(v2,v2))
       # pick random magnitudes
       rv1 = rng.normal(step_size)
       rv2 = rng.normal(step_size)
       # create new cell and transformation matrix
       new_cell = orig_cell.copy()
       new_cell[rnd_vec_ind,:] += rv1*v1 + rv2*v2
   transform = np.dot(np.linalg.inv(orig_cell), new_cell)
   return (1.0, transform)

def propose_stretch_step(at, step_size):
    if movement_args['2D']:
        rnd_v1_ind = 0
        rnd_v2_ind = 1
    else:
        rnd_v1_ind = rng.int_uniform(0, 3)
        rnd_v2_ind = rng.int_uniform(0, 3)
        if rnd_v1_ind == rnd_v2_ind:
            rnd_v2_ind = (rnd_v2_ind+1) % 3

    rv = rng.normal(step_size)
    transform = np.identity(3)
    transform[rnd_v1_ind,rnd_v1_ind] = np.exp(rv)
    transform[rnd_v2_ind,rnd_v2_ind] = np.exp(-rv)
    return (1.0, transform)

def min_aspect_ratio(vol, cell):
    min_aspect_ratio = sys.float_info.max
    nD=3
    if movement_args['2D']:
        # in 2D, do not check the Z direction
        nD=2
    for i in range(nD):
        vi = cell[i,:]
        vnorm_hat = np.cross(cell[(i+1)%3,:],cell[(i+2)%3,:])
        vnorm_hat /= np.sqrt(np.dot(vnorm_hat,vnorm_hat))
        min_aspect_ratio = min(min_aspect_ratio, abs(np.dot(vnorm_hat,vi)))

    if movement_args['2D']:
        # in 2D divide by the square root of area
        min_aspect_ratio /= (vol/np.sqrt(np.dot(cell[2,:],cell[2,:])))**(1.0/2.0)
    else:
        min_aspect_ratio /= vol**(1.0/3.0)

    return min_aspect_ratio

def do_cell_step(at, Emax, p_accept, transform):
    if p_accept < 1.0 and rng.float_uniform(0.0,1.0) > p_accept:
        return False

    # save old config and apply transformation
    orig_cell = at.get_cell()
    new_cell = np.dot(orig_cell,transform)
    new_vol = abs(np.dot(new_cell[0,:],np.cross(new_cell[1,:],new_cell[2,:])))

    # check size and shape constraints
    if new_vol > ns_args['max_volume_per_atom']*len(at) or new_vol < ns_args['min_volume_per_atom']*len(at):
        return False
    if min_aspect_ratio(new_vol, new_cell) < movement_args['MC_cell_min_aspect_ratio']:
        return False

    # orig_cell already saved previous value
    orig_pos = at.get_positions()
    if ns_args['n_extra_data'] > 0:
        extra_data = at.arrays['ns_extra_data'].copy()

    # set new positions and velocities
    at.set_cell(new_cell, scale_atoms=True)

    if Emax is None:
        return

    # calculate new energy
    try:
        new_energy = eval_energy(at)
    except Exception as err:
        if ns_args['debug'] >= 4:
            print( "eval_energy got exception ", err)
        new_energy = 2.0*abs(Emax)
        #print("error in eval_energy setting new_energy = 2*abs(Emax)=" , new_energy)

    # accept or reject
    if new_energy < Emax: # accept
        at.info['ns_energy'] = new_energy
        return True
    else: # reject and revert
        at.set_cell(orig_cell,scale_atoms=False)
        at.set_positions(orig_pos)
        if ns_args['n_extra_data'] > 0:
            at.arrays['ns_extra_data'][...] = extra_data
        return False

def do_cell_shape_walk(at, movement_args):
    possible_moves = {
        'MC_cell_shear': propose_shear_step,
        'MC_cell_stretch': propose_stretch_step }
    items = list(possible_moves.items())

    for i in range(movement_args['cell_shape_equil_steps']):
        rng.shuffle_in_place(items)
        for key, propose_step_func in items:
            (p_accept, transform) = propose_step_func(at, movement_args[key+"_step_size"])
            do_cell_step(at, None, p_accept, transform)

def do_MC_semi_grand_step(at, movement_args, Emax, KEmax):
    global Z_list

    Z = at.get_atomic_numbers()
    n_types = len(movement_args['semi_grand_potentials'])
    at_i = rng.int_uniform(0,len(at))
    type_i = rng.int_uniform(0,n_types)
    while Z_list[type_i] == Z[at_i]:
        type_i = rng.int_uniform(0,n_types)
    Z_new = Z.copy()
    Z_new[at_i] = Z_list[type_i]
    at.set_atomic_numbers(Z_new)

    new_energy = eval_energy(at)
    new_KE = eval_energy_KE(at)
    if new_energy < Emax and (KEmax < 0.0 or new_KE < KEmax): # accept
        # update energy
        at.info['ns_energy'] = new_energy
        n_accept = 1
    else:
        # undo
        at.set_atomic_numbers(Z)
        n_accept = 0

    return (1, {'MC_semi_grand' : (1, n_accept) })

def do_MC_swap_step(at, movement_args, Emax, KEmax):
#DOC
#DOC ``do_MC_swap_step``
    Z = at.get_atomic_numbers()

    #DOC \item return if all atomic numbers are identical
    if (Z[:] == Z[0]).all():
        # don't try to swap when all atoms are the same
        return (0, {})

    r_cut = movement_args['swap_r_cut']
    #DOC \item randomly pick a desired cluster size
    cluster_size = np.where(rng.float_uniform(0,1) < movement_args['swap_probs'])[0][0]+1
    if cluster_size > 1:
        (i_list, j_list) = matscipy.neighbours.neighbour_list('ij', at, r_cut)
    else:
        i_list = None
        j_list = None
    #DOC \item pick two clusters with distinct atomic numbers, backing off to smaller clusters on failure to find appropriate larger clusters, but always pick at least a pair of atoms to swap
    c1 = None
    c2 = None
    while (c1 is None or c2 is None or np.all(Z[c1] == Z[c2])):
        # print(print_prefix, ": do_MC_swap try to find cluster ", cluster_size)
        c1 = pick_interconnected_clump.pick_interconnected(rng, len(at), i_list, j_list, cluster_size, r_cut)
        c2 = pick_interconnected_clump.pick_interconnected(rng, len(at), i_list, j_list, cluster_size, r_cut)
        # print(print_prefix, ": do_MC_swap got ", c1, c2)
        # decrement cluster size
        cluster_size -= 1
        if cluster_size < 1:
            cluster_size = 1

    # failed to find appropriate
    if c1 is None or c2 is None or np.all(Z[c1] == Z[c2]):
        # print(print_prefix, ": do_MC_swap giving up on cluster ",c1,c2)
        # return 1 for number of model calls so walk will finish, even if no clusters are ever found
        return (1, {})

    # print(print_prefix, ": do_MC_swap trying cluster", c1, c2)
    p_1_orig = at.positions[c1,:].copy()
    p_2_orig = at.positions[c2,:].copy()
    at.positions[c1,:] = p_2_orig
    at.positions[c2,:] = p_1_orig
    if not movement_args['swap_velo']:
        # TODO: don't do this for the time being, but we should be able to modify this
        # If we just swap particles by swapping positions, velocities will follow to the new position, and so will effectively be swapped
        # Therefore, if we _don't_ want to actually swap velocities, we have to (apparently) swap them
        velocities = at.get_velocities()
        v_1_orig = velocities[c1,:].copy()
        v_2_orig = velocities[c2,:].copy()
        if movement_args['no_swap_velo_fix_mag_alt']:
            # make velocity have the same direction, but magnitude from the previous velocity of the other particle so local KE and momentum are preserved
            n1 = np.linalg.norm(v_1_orig)
            n2 = np.linalg.norm(v_2_orig)
            velocities[c1,:] = (v_2_orig.T * n1/n2).T
            velocities[c2,:] = (v_1_orig.T * n2/n1).T
        else: # this will be executed (unnecessarily) even if all masses are the same
            # make velocity have same direction, and magnitude rescaled to keep local kinetic energy the same
            masses = at.get_masses()
            velocities[c1,:] = (v_2_orig.T * np.sqrt(masses[c2]/masses[c1])).T
            velocities[c2,:] = (v_1_orig.T * np.sqrt(masses[c1]/masses[c2])).T
        at.set_velocities(velocities)
    if ns_args['n_extra_data'] > 0:
        extra_data_1_orig = at.arrays['ns_extra_data'][c1,...].copy()
        extra_data_2_orig = at.arrays['ns_extra_data'][c2,...].copy()
        at.arrays['ns_extra_data'][c1,...] = extra_data_2_orig
        at.arrays['ns_extra_data'][c2,...] = extra_data_1_orig

    #DOC \item accept swap if energy < Emax
    new_energy = eval_energy(at)
    new_KE = eval_energy_KE(at)

    if new_energy < Emax and (KEmax < 0.0 or new_KE < KEmax): # accept
        at.info['ns_energy'] = new_energy
        accept_n = 1
    else: # reject
        at.positions[c1,:] = p_1_orig
        at.positions[c2,:] = p_2_orig
        if not movement_args['swap_velo']:
            velocities[c1,:] = v_1_orig
            velocities[c2,:] = v_2_orig
            at.set_velocities(velocities)
        if ns_args['n_extra_data'] > 0:
            at.arrays['ns_extra_data'][c1,...] = extra_data_1_orig
            at.arrays['ns_extra_data'][c2,...] = extra_data_2_orig
        accept_n = 0


    return (1, {('MC_swap_%d' % len(c1)) : (1, accept_n) })

def do_MC_cell_volume_step(at, movement_args, Emax, KEmax):
#DOC
#DOC ``do_MC_cell_volume_step``
    step_rv = rng.float_uniform(0.0, 1.0)
    if step_rv > movement_args['MC_cell_volume_per_atom_prob']:
        return (0, {})
    if movement_args['2D']:
        (p_accept, transform) = propose_area_step(at, movement_args['MC_cell_volume_per_atom_step_size'], movement_args['MC_cell_flat_V_prior'])
    else:
        (p_accept, transform) = propose_volume_step(at, movement_args['MC_cell_volume_per_atom_step_size'], movement_args['MC_cell_flat_V_prior'])

    if do_cell_step(at, Emax, p_accept, transform):
        return (1, {'MC_cell_volume_per_atom' : (1, 1) })
    else:
        return (1, {'MC_cell_volume_per_atom' : (1, 0) })

def do_MC_cell_shear_step(at, movement_args, Emax, KEmax):
#DOC
#DOC ``do_MC_cell_shear_step``
    step_rv = rng.float_uniform(0.0, 1.0)
    if step_rv > movement_args['MC_cell_shear_prob']:
        return (0, {})
    (p_accept, transform) = propose_shear_step(at, movement_args['MC_cell_shear_step_size'])
    if do_cell_step(at, Emax, p_accept, transform):
        return (1, {'MC_cell_shear' : (1, 1) })
    else:
        return (1, {'MC_cell_shear' : (1, 0) })

def do_MC_cell_stretch_step(at, movement_args, Emax, KEmax):
#DOC
#DOC ``do_MC_cell_stretch_step``
    step_rv = rng.float_uniform(0.0, 1.0)
    if step_rv > movement_args['MC_cell_stretch_prob']:
        return (0, {})
    (p_accept, transform) = propose_stretch_step(at, movement_args['MC_cell_stretch_step_size'])
    if do_cell_step(at, Emax, p_accept, transform):
        return (1, {'MC_cell_stretch' : (1, 1) })
    else:
        return (1, {'MC_cell_stretch' : (1, 0) })


def do_atom_walk(at, movement_args, Emax, KEmax):
#DOC
#DOC ``do_atom_walk``
    n_reps = movement_args['n_atom_steps_per_call']
    out = {}
    #DOC \item loop n\_atom\_steps\_per\_call times, calling do\_MC\_atom\_walk() or do\_MD\_atom\_walk()
    for i in range(n_reps):
        if movement_args['atom_algorithm'] == 'MC':
            accumulate_stats(out, do_MC_atom_walk(at, movement_args, Emax, KEmax))
        elif movement_args['atom_algorithm'] == 'MD':
            accumulate_stats(out, do_MD_atom_walk(at, movement_args, Emax, KEmax))
        else:
            exit_error("do_atom_walk got unknown 'atom_algorithm' = '%s'\n" % movement_args['atom_algorithm'], 5)

    return (n_reps*movement_args['atom_traj_len']*movement_args['atom_traj_len_cost_multiplier'], out)

def rand_perturb_energy(energy, perturbation, Emax=None):
    if Emax is None:
        if abs(energy) <= 1.0:
            energy += rng.float_uniform(-1.0, 0.0)*perturbation
        else:
            energy *= (1.0+rng.float_uniform(-1.0, 0.0)*perturbation)
    else:
        if abs(energy) <= 1.0:
            pert = rng.float_uniform(-1.0, 0.0)*perturbation
            n_tries = 0
            while energy+pert >= Emax and n_tries < 100:
                pert = rng.float_uniform(-1.0, 0.0)*perturbation
                n_tries += 1
            if energy+pert >= Emax:
                print(print_prefix, "WARNING: failed to do random energy perturbation below Emax ", energy, Emax)
            energy += pert
        else:
            pert = 1.0 + rng.float_uniform(-1.0, 0.0)*perturbation
            n_tries = 0
            while energy*pert >= Emax and n_tries < 100:
                pert = 1.0 + rng.float_uniform(-1.0, 0.0)*perturbation
                n_tries += 1
            if energy*pert >= Emax:
                print(print_prefix, "WARNING: failed to do random energy perturbation below Emax ", energy, Emax)
            energy *= pert

    return energy

####################################################################################################

def add_to_list(list, ind, costs, nums, walk_len_avail, first_free_slot):

    if nums[0] == 0:
        return (walk_len_avail, first_free_slot)

    n_steps_r = float(walk_len_avail) * (float(costs[0]*nums[0]) / float(sum(costs*nums))) / float(costs[0])

    n_steps_r_fract = n_steps_r - math.floor(n_steps_r)
    rv = np.random.rand()
    if rv < n_steps_r_fract:
        n_steps = int(math.floor(n_steps_r)+1)
    else:
        n_steps = int(math.floor(n_steps_r))

    list[first_free_slot:first_free_slot+n_steps] = [ind]*n_steps
    first_free_slot += n_steps
    walk_len_avail -= n_steps*costs[0]

    return (walk_len_avail, first_free_slot)

def create_list(costs, nums, walk_len):

    list = []
    walk_len_avail = walk_len
    first_free_slot = 1
    for i in range(len(costs)):
        (walk_len_avail, first_free_slot) = add_to_list(list, i, costs[i:], nums[i:], walk_len_avail, first_free_slot)

    return list

####################################################################################################

def walk_single_walker(at, movement_args, Emax, KEmax):
    """Do random walk on a single atoms object."""
#DOC
#DOC ``walk_single_walker``

    out = {}

    if movement_args['do_good_load_balance']:
        possible_moves = np.array([do_atom_walk,
                                   do_MC_cell_volume_step,
                                   do_MC_cell_shear_step,
                                   do_MC_cell_stretch_step,
                                   do_MC_swap_step,
                                   do_MC_semi_grand_step])
        nums = np.array([movement_args['n_atom_steps_n_calls'],
                         movement_args['n_cell_volume_steps'],
                         movement_args['n_cell_shear_steps'],
                         movement_args['n_cell_stretch_steps'],
                         movement_args['n_swap_steps'],
                         movement_args['n_semi_grand_steps']])
        costs = np.array([movement_args['atom_traj_len']*movement_args['atom_traj_len_cost_multiplier'],
                          1,
                          1,
                          1,
                          1,
                          1])

        list = create_list(costs, nums, movement_args['n_model_calls'])
        for move_i in list:
            (t_n_model_calls, t_out) = possible_moves[move_i](at, movement_args,
                                                              Emax, KEmax)
            accumulate_stats(out, t_out)

    else:
        #DOC \item create list
                            #DOC \item do\_atom\_walk :math:`*` n\_atom\_step\_n\_calls
        possible_moves = ([do_atom_walk] * movement_args['n_atom_steps_n_calls'] +
                            #DOC \item do\_cell\_volume\_step :math:`*` n\_cell\_volume\_steps
                          [do_MC_cell_volume_step] * movement_args['n_cell_volume_steps'] +
                            #DOC \item do\_cell\_shear\_step :math:`*` n\_cell\_shear\_steps
                          [do_MC_cell_shear_step] * movement_args['n_cell_shear_steps'] +
                            #DOC \item do\_cell\_stretch\_step :math:`*` n\_cell\_stretch\_steps
                          [do_MC_cell_stretch_step] * movement_args['n_cell_stretch_steps'] +
                            #DOC \item do\_swap\_step :math:`*` n\_swap\_steps
                          [do_MC_swap_step] * movement_args['n_swap_steps'] +
                            #DOC \item do\_semi\_grand\_step :math:`*` n\_semi\_grand\_steps
                          [do_MC_semi_grand_step] * movement_args['n_semi_grand_steps'])

        out = {}
        n_model_calls_used = 0

        #DOC \item loop while n\_model\_calls\_used < n\_model\_calls
        while n_model_calls_used < movement_args['n_model_calls']:
            #DOC \item pick random item from list
            move = possible_moves[rng.int_uniform(0, len(possible_moves))]
            #DOC \item do move
            (t_n_model_calls, t_out) = move(at, movement_args, Emax, KEmax)
            n_model_calls_used += t_n_model_calls
            accumulate_stats(out, t_out)

    #DOC \item perturb final energy by random\_energy\_perturbation
    # perturb final energy
    at.info['ns_energy'] = rand_perturb_energy(
        at.info['ns_energy'], ns_args['random_energy_perturbation'], Emax)

    # DEBUG print("walk_single_walker end ", eval_energy(at, do_PE=False),
    # eval_energy(at) #DEBUG)

    return out


def max_energy(walkers, n, kinetic_only=False):
    """Collect the current energies of the walkers from all the processes and
    chooses the right number of highest energies to be culled"""
    # do local max
    if kinetic_only:
        energies_loc = np.array([eval_energy_KE(at) for at in walkers])
    else:
        energies_loc = np.array([ at.info['ns_energy'] for at in walkers])
    volumes_loc = np.array([ at.get_volume() for at in walkers])
    if comm is not None:
        energies = np.zeros( (comm.size*len(energies_loc)) )
        volumes = np.zeros( (comm.size*len(volumes_loc)) )
        # comm.barrier() #BARRIER
        comm.Allgather( [ energies_loc, MPI.DOUBLE ], [ energies, MPI.DOUBLE ] )
        energies = energies.flatten()
        comm.Allgather( [ volumes_loc, MPI.DOUBLE ], [ volumes, MPI.DOUBLE ] )
        volumes = volumes.flatten()
    else:
        energies = energies_loc
        volumes = volumes_loc

    # n is n_cull
    Emax_ind = energies.argsort()[-1:-n-1:-1]
    Emax = energies[Emax_ind]
    Vmax = volumes[Emax_ind]
    # WARNING: assumes that each node has equal number of walkers
    rank_of_max = np.floor(Emax_ind/len(walkers)).astype(int)
    ind_of_max = np.mod(Emax_ind,len(walkers))

    return (Emax, Vmax, rank_of_max, ind_of_max)

# WARNING: full_auto_set_stepsizes shouldn't really depend on walk_stats from a real walk, since it does all its own pilot walks
# right now it uses this data structure to figure out what kind of steps actually occur, but it should probably get this
# information right from movement_args, rather than relying on what happened to have happened in the real walks
def full_auto_set_stepsizes(walkers, walk_stats, movement_args, comm, Emax, KEmax, size_n_proc):
    """Automatically set all step sizes. Returns the time (in seconds) taken for the routine to run."""
#DOC
#DOC ``full_auto_set_stepsizes``
    #DOC \item Step sizes for each (H)MC move are set via a loop which performs additional exploration moves, calibrating each step size to obtain an acceptance rate inside a specified range.

    global print_prefix

    orig_prefix = print_prefix
    print_prefix = "full_auto"
    if comm is not None:
        print_prefix = "%d %s" % (comm.rank, print_prefix)

    full_auto_start_time = time.time()
    n_samples_per_move_type=200 # total number of (H)MC moves used to set each step length

    if comm is not None:
        walk_n_walkers = int(np.ceil(float(n_samples_per_move_type)/size_n_proc))
        # in each trial we will evolve walk_n_walkers configurations
    else:
        walk_n_walkers = n_samples_per_move_type
    #DOC \item The routine is MPI parallelised, so that the wall time goes as 1/num\_of\_processes

    key_list=[]
    for key, value in walk_stats.items():
        key_list.append(key)

    if (comm is not None):
        # Identify move types that are being used
        # Not all processes may have used the same move types
        key_ints=[]
        for key in key_list:
            if (key == "MD_atom"):
                i = 0
            elif (key == "MC_atom"):
                i = 1
            elif (key == "MC_atom_velo"):
                i = 2
            elif (key == "MC_cell_shear"):
                i = 3
            elif (key == "MC_cell_stretch"):
                i = 4
            elif (key == "MC_cell_volume_per_atom"):
                i = 5
            elif (key[:7] == "MC_swap"):
                i = 6
            else:
                i = -1
            key_ints.append(i)

        key_flags =  1*np.asarray([ i in key_ints for i in range(7)])

        totalkeys = np.zeros( (len(key_flags)), dtype=int)
        comm.Allreduce([key_flags, MPI.INT], [totalkeys, MPI.INT], MPI.SUM)

        key_list = []
        for i in range(7):
            if (totalkeys[i]>0):
                if (i==0):
                    key_list.append("MD_atom")
                if (i==1):
                    key_list.append("MC_atom")
                if (i==2):
                    key_list.append("MC_atom_velo")
                if (i==3):
                    key_list.append("MC_cell_shear")
                if (i==4):
                    key_list.append("MC_cell_stretch")
                if (i==5):
                    key_list.append("MC_cell_volume_per_atom")
                if (i==6):
                    key_list.append("MC_swap_")

    #DOC \item For each (H)MC move type the following is performed
    for key in key_list:

        #DOC \item Set ''movement\_args'' parameters so that only one (H)MC call is made at a time
        # reprogram n_atom_steps_n_calls, n_cell_volume_steps, n_cell_shear_steps, n_cell_stretch_steps, n_swap_steps according to key
        # possible key values:
        # MD_atom # atoms HMC
        # MC_atom # MC atom sweeps
        # MC_atom_velo # MC velocity sweeps
        # MC_cell_shear # cell shear move
        # MC_cell_stretch # cell stretch move
        # MC_cell_volume_per_atom # volume move

        exploration_movement_args = deepcopy(movement_args)
        # turn all step types off to begin with
        exploration_movement_args['n_atom_steps_n_calls'] = 0
        exploration_movement_args['n_cell_volume_steps'] = 0
        exploration_movement_args['n_cell_shear_steps'] = 0
        exploration_movement_args['n_cell_stretch_steps'] = 0
        exploration_movement_args['n_swap_steps'] = 0
        exploration_movement_args['n_semi_grand_steps'] = 0
        exploration_movement_args['MC_atom_velocities']=False

        nD=3
        if exploration_movement_args['2D']:
           nD=2

        # check that the total number of attempts for this key is not zero
        if key in walk_stats:
            (n_try, n_accept) = walk_stats[key]
        else:
            n_try=0
        n_try_g = np.zeros( (1), dtype=int)
        if (comm is not None):
            n_try_s = np.array( [n_try], dtype = int)
            comm.Allreduce([n_try_s, MPI.INT], [n_try_g, MPI.INT], MPI.SUM)
        else:
            n_try_g[0] = n_try

        if (n_try_g[0]==0):
            continue # skip this key - it is not used

        if (key=="MD_atom" or key=="MC_atom"):
            exploration_movement_args['n_atom_steps_n_calls'] = 1
            # one call to do_atom_walk per walk_single_walker call

            exploration_movement_args['n_atom_steps_per_call'] = 1
            # do_atom_walk makes one do_MC_walk/do_MD_walk per call

            if key == "MC_atom" and not movement_args['MC_atom_Galilean']:
                exploration_movement_args['atom_traj_len'] = 1
                # 1 atom sweep per do_MC_walk call
                exploration_movement_args['n_model_calls'] = 1 * exploration_movement_args['atom_traj_len_cost_multiplier']
                # The do_atom_walk routine reports that it has performed
                # #model_calls = the number of complete MC sweeps performed,
                # rather than single point evaluations.
            else:
                exploration_movement_args['n_model_calls'] = exploration_movement_args['atom_traj_len']*exploration_movement_args['atom_traj_len_cost_multiplier']
                # The do_atom_walk routine reports that it has performed
                # #model_calls = number of single point evaluations (time steps)

        elif (key=="MC_atom_velo"):
            exploration_movement_args['velo_traj_len']=1
            exploration_movement_args['MC_atom_velo_walk_rej_free']=False
            # one velocities sweep

        elif (key=="MC_cell_shear"):
            exploration_movement_args['n_cell_shear_steps'] = 1
            exploration_movement_args['n_model_calls'] = 1
            # one call to do_MC_cell_shear_step per walk_single_walker call

        elif (key=="MC_cell_stretch"):
            exploration_movement_args['n_cell_stretch_steps'] = 1
            exploration_movement_args['n_model_calls'] = 1
            # one call to do_MC_cell_stretch_step per walk_single_walker call

        elif (key=="MC_cell_volume_per_atom"):
            exploration_movement_args['n_cell_volume_steps'] = 1
            exploration_movement_args['n_model_calls'] = 1
            # one call to do_MC_cell_volume_step per walk_single_walker call
        elif (key[:7] == "MC_swap"):
            # skip swap moves, since they have no step size
            break

        else:
            exit_error("full_auto_set_stepsizes got key '%s', unkown to this routine\n" % key, 5)

        #DOC \item Min and max acceptance rates are copied from parameters MC\_adjust\_min\_rate / MD\_adjust\_min\_rate and MC\_adjust\_max\_rate / MD\_adjust\_max\_rate

        if key.find("MC") == 0:
            if key.find("MC_atom_step_size") and movement_args['MC_atom_Galilean']:
                min_rate = movement_args['GMC_adjust_min_rate']
                max_rate = movement_args['GMC_adjust_max_rate']
            else:
                min_rate = movement_args['MC_adjust_min_rate']
                max_rate = movement_args['MC_adjust_max_rate']
            suffix="step_size"
        elif key.find("MD") == 0:
            min_rate = movement_args['MD_adjust_min_rate']
            max_rate = movement_args['MD_adjust_max_rate']
            suffix="timestep"
        else:
            exit_error("full_auto_set_stepsizes got key '%s', neither MC nor MD\n" % key, 5)

        first_walker = rng.int_uniform(0, len(walkers)) # random starting point for config selection
        first_time = True # we will make at least two tries. Logical flag ensures this.

        steplength_store = movement_args[key+"_"+suffix]
        # protects against possible future bugs that would be hard to detect


        #DOC \item Step size calibration loop:
        dir = None
        while True:
            stats = {} # clear acceptance / trial stats for new step size
            stats_cumul = {} # clear acceptance / trial stats for new step size
            #DOC \item Repeat the following 200/num\_of\_MPI\_processes times:
                #DOC \item Copy a configuration from the live set (each MPI process chooses a different configuration)
            for i in range(first_walker,first_walker + walk_n_walkers):

                k = i%len(walkers) # cycle through walkers array
                buf = walkers[k].copy() # copy config k into buffer "buf" for walking (walkers array unchanged)
                buf.calc = walkers[k].calc # copy calculator

                #DOC \item Each MPI processes performs one (H)MC move on its cloned configuration
                # build up stats from walkers
                if (not key=="MC_atom_velo"):
                    stats = walk_single_walker(buf, exploration_movement_args, Emax, KEmax)
                else:
                    stats = do_MC_atom_velo_walk(buf, exploration_movement_args, Emax, nD, KEmax)

                  #DOC     running statistics for the number of accepted/rejected moves on each process are recorded
                accumulate_stats(stats_cumul, stats)

            first_walker = first_walker + walk_n_walkers # cycle through local samples
            (n_try, n_accept) = stats_cumul[key]

            if comm is not None:
                n_try_s = np.array( [n_try], dtype = int)
                n_accept_s = np.array( [n_accept], dtype = int)
                n_try_g = np.zeros( (1), dtype=int)
                n_accept_g = np.zeros( (1), dtype=int)
                # comm.barrier() #BARRIER
                comm.Allreduce([n_try_s, MPI.INT], [n_try_g, MPI.INT], MPI.SUM)
                comm.Allreduce([n_accept_s, MPI.INT], [n_accept_g, MPI.INT], MPI.SUM)
                n_try = n_try_g[0]
                n_accept = n_accept_g[0]

            rate = float(n_accept)/float(n_try)
            #DOC \item The total number of accepted/rejected moves for this step size (summed across all MPI processes) are estabilshed

            if ((comm is None or comm.rank == 0) and (ns_args['debug'] >= 1)):
                print( print_prefix, "trial stepsize and accept rate for %s = %e , %f (%d)" % (key, movement_args[key+"_"+suffix], rate, n_try))

            if (rate>min_rate and rate<max_rate):
                #DOC \item If the total acceptance rate is within the desired range, return this stepsize
                if (comm is None or comm.rank == 0):
                    print( print_prefix, "full_auto_set_stepsizes adjusted %s to %f" % (key+"_"+suffix, movement_args[key+"_"+suffix]))
                break
            else:
                if( not first_time ): # dodge this the first time round
                    # Check whether rate and rate_store are on different sides
                    # of interval
                    if ((min(rate,rate_store) < min_rate) and (max(rate,rate_store)>max_rate)):
                        #DOC \item If this is NOT the first time the loop has been performed for this (H)MC move AND we previously obtained an acceptance rate on one side of the desired range, and now find an acceptance rate on the other side of the desired range
                            #DOC \item Return the step size that gave an acceptance rate closest to the middle of the desired range.

                        # check which gives a accept_rate closer to centre of acceptable window
                        # and take that
                        target = 0.5*(min_rate+max_rate)
                        if (abs(rate-target)<abs(rate_store-target)):
                            # take current step length
                            if (comm is None or comm.rank == 0):
                                print( print_prefix, "full_auto_set_stepsizes adjusted %s to %f" % (key+"_"+suffix , movement_args[key+"_"+suffix]))
                            break
                        else:
                            # take saved step length
                            movement_args[key+"_"+suffix] = steplength_store
                            exploration_movement_args[key+"_"+suffix] = steplength_store
                            rate = rate_store
                            if (comm is None or comm.rank == 0):
                                print( print_prefix, "full_auto_set_stepsizes adjusted %s to %f" % (key+"_"+suffix, movement_args[key+"_"+suffix]))
                            break

                else: # this is the first time
                    first_time = False

                #DOC \item Store step length and acceptance rate
                # save current step length and acceptance rate
                steplength_store = movement_args[key+"_"+suffix]
                rate_store = rate

                #DOC \item update step length, by :math:`*` or :math:`/` by MC\_adjust\_step\_factor, to improve acceptance rate
                #update step length
                dir = None
                if rate < min_rate:
                    exp = -1.0
                    dir = "down"
                elif rate >= max_rate:
                    exp = 1.0
                    dir = "up"
                else:
                    exp = None

                # try to adjust
                if dir is not None:
                    movement_args[key+"_"+suffix] *= movement_args['MC_adjust_step_factor']**exp
                    exploration_movement_args[key+"_"+suffix] *= exploration_movement_args['MC_adjust_step_factor']**exp
                    if (comm is None or comm.rank == 0) and (ns_args['debug'] >= 1):
                        print( print_prefix, "new trial step size for %s = %e" % (key, movement_args[key+"_"+suffix]))

                #DOC \item Check that step size is not larger than max allowed value (specified by user), and also that step size is not smaller than 10\ :sup:`-20`\ (useful for detecting errors).
                # if exceeded maximum, cap change
                if movement_args[key+"_"+suffix] > movement_args[key+"_"+suffix+"_max"]:
                    movement_args[key+"_"+suffix] = movement_args[key+"_"+suffix+"_max"]
                    exploration_movement_args[key+"_"+suffix] = exploration_movement_args[key+"_"+suffix+"_max"]
                # Error check:
                if (movement_args[key+"_"+suffix] < 1.0e-20):
                    exit_error("full_auto_set_stepsizes got '%s'_'%s' = '%e': too small. Is everything correct?\n" % (key, suffix, movement_args[key+"_"+suffix]), 25)

                # if was already at max and didn't really change, break
                if (movement_args[key+"_"+suffix] == steplength_store):
                    dir = None
                    if (comm is None or comm.rank == 0):
                        print( print_prefix, "full_auto_set_stepsizes adjusted %s to %f" % (key+"_"+suffix, movement_args[key+"_"+suffix]))
                    break

    #DOC \item Return step sizes and time taken for routine to run

    print_prefix = orig_prefix

    full_auto_end_time = time.time()
    duration = full_auto_end_time - full_auto_start_time
    return duration

def adjust_step_sizes(walk_stats, movement_args, comm, do_print_rate=True, monitor_only=False):
    """
    Adjust step size to keep the acceptance ratio at the desired level.
    """
    for key in walk_stats:
        (n_try, n_accept) = walk_stats[key]
        if comm is not None:
            n_try_s = np.array( [n_try], dtype = int)
            n_accept_s = np.array( [n_accept], dtype = int)
            n_try_g = np.zeros( (1), dtype=int)
            n_accept_g = np.zeros( (1), dtype=int)
            # comm.barrier() #BARRIER
            comm.Allreduce([n_try_s, MPI.INT], [n_try_g, MPI.INT], MPI.SUM)
            comm.Allreduce([n_accept_s, MPI.INT], [n_accept_g, MPI.INT], MPI.SUM)
            n_try = n_try_g[0]
            n_accept = n_accept_g[0]

        if n_try > 0:
            rate = float(n_accept)/float(n_try)

            if do_print_rate and comm is None or comm.rank == 0:
                print( print_prefix, "accept rate for %s = %f (%d)" % (key, rate, n_try))

            if monitor_only:
                continue

            if key.find("MC") == 0:
                if key.find("MC_atom_step_size") and movement_args['MC_atom_Galilean']:
                    min_rate = movement_args['GMC_adjust_min_rate']
                    max_rate = movement_args['GMC_adjust_max_rate']
                else:
                    min_rate = movement_args['MC_adjust_min_rate']
                    max_rate = movement_args['MC_adjust_max_rate']
                suffix="step_size"
            elif key.find("MD") == 0:
                min_rate = movement_args['MD_adjust_min_rate']
                max_rate = movement_args['MD_adjust_max_rate']
                suffix="timestep"
            else:
                if (comm is None or comm.rank == 0):
                    print( "WARNING: adjust_step_size got key '%s', neither MC nor MD\n" % key)
                continue

            if key+"_"+suffix not in movement_args:
                continue

            dir = None
            exp = 0.0
            if rate < min_rate:
                exp = -1.0
                dir = "down"
            elif rate >= max_rate:
                exp = 1.0
                dir = "up"

            orig_value = movement_args[key+"_"+suffix]

            if dir is not None:
                # try to adjust
                movement_args[key+"_"+suffix] *= movement_args['MC_adjust_step_factor']**exp
                # if exceeded maximum, cap change
                if movement_args[key+"_"+suffix] > movement_args[key+"_"+suffix+"_max"]:
                    movement_args[key+"_"+suffix] = movement_args[key+"_"+suffix+"_max"]

            # if was already at max and didn't really change, unset dir
            if movement_args[key+"_"+suffix] == orig_value:
                dir=None

            if dir is not None and (comm is None or comm.rank == 0):
                print( print_prefix, "adjust_step_sizes adjusted %s %s to %f" % (key, dir, movement_args[key+"_"+suffix]))


def zero_stats(d, movement_args):
    for key in movement_args:
        m = re.search('(.*)_(step_size|timestep)$', key)
        if m is not None:
            d[m.group(1)] = (0, 0)


def accumulate_stats(d_cumul, d):
    for key in d:
        if key in d_cumul:
            d_cumul[key] = tuple([i1+i2 for i1, i2 in zip(d[key], d_cumul[key])])
        else:
            d_cumul[key] = d[key]


# figure out n_steps to walk on each iteration to get correct expected number
def set_n_from_expected(prop):
    if movement_args[prop+'_expected'] > 0:
        if movement_args[prop] > 0:
            exit_error("Got both "+prop+" and "+prop+"_expected, conflict\n", 5)

        if rank == 0:
            print( "Calculating %s from %s_expected=%d" % (prop, prop, movement_args[prop+'_expected']))

        if max_n_cull_per_task*size == n_cull and n_extra_walk_per_task == 0: # no extra walkers
            movement_args[prop] = movement_args[prop+'_expected']
            if rank == 0:
                print( "No extra walkers (n_cull mod n_tasks == 0), trivial, so average n_walks at kill is 1, and %s=%s_expected" % (prop, prop))
        else:
            # f_c = n_c/n_t [ fraction of total that are culled (and walked once)]
            # f_e = n_e/(n_t-n_c) [fraction of ones that aren't culled that are also walked ]

            # n(1) <- f_c n_t + (1-f_c)*(1-f_e) n(1)
            # n(l >= 2) <- (1-f_c)(1-f_e) n(l) + (1-f_c) f_e n(l-1)

            # f = (1-f_c)f_e / (f_c+f_e - f_c f_e)

            # n(1)/n_t = f_c/(1-(1-f_c)(1-f_e)) = f_c/(f_c+f_e-f_c f_e)
            # n(l >= 2)/n_t = f n(l-1) = f^(l-1) n(1)/n_t
            # n_walks = n(1)/n_t + sum_{l=2}^\infty l n(l)/n_t
            #         = n(1)/n_t (1 + sum_{l=2}^\infty l f^{l-1})
            #         = n(1)/n_t (1 + \sum{l=1}^\infty (l+1) f^l)
            #         = n(1)/n_t (1 + \sum f^l + sum l f^l)
            #         = n(1)/n_t (1 + f/(1-f) + f/(1-f)**2

            f_cull = float(n_cull)/float(ns_args['n_walkers'])
            f_extra = float( (max_n_cull_per_task*size-n_cull) + n_extra_walk_per_task*size ) / float(ns_args['n_walkers']-n_cull)
            f = (1.0-f_cull)*f_extra/(f_cull+f_extra-f_cull*f_extra)

            n_walks = f_cull/(f_cull+f_extra-f_cull*f_extra) * (1.0 + f/(1.0-f) + f/(1.0-f)**2)
            if rank == 0:
                print("f_cull" ,f_cull,"f_extra", f_extra, "f", f)
                print( "Calculated average n_walks at kill = ",n_walks, " from f_cull ",f_cull," f_extra (due to otherwise idle processors doing walks and explicitly requested extra walks) ",f_extra)
                if n_walks < 1:
                    print( "WARNING: got average n_walks < 1, but will always do at least 1 walk, so effective %s_expected will be higher than requested" % prop)
                print( "Setting %s = ceiling(%s_expected/n_walks)" % (prop, prop))
            movement_args[prop] = int(math.ceil(movement_args[prop+'_expected']/n_walks))

    else:
        if movement_args[prop] > 0 and rank == 0:
            print( "WARNING: using absolute number of "+prop)

    if rank == 0:
        print( "Final value of %s=%d" % (prop, movement_args[prop]))

def additive_init_config(at, Emax):
    if do_calc_lammps:
        exit_error("python additive_init_config doesn't work with LAMMPS, since it varies list of atoms\n", 10)
    pos = at.get_positions()
    for i_at in range(1,len(at)):
        at_new = at[0:i_at+1]
        if do_calc_ASE or do_calc_lammps:
            at_new.calc = at.calc
        at_new.set_positions(pos[0:i_at+1,:])
        success = False
        for i_try in range(10):
            pos[i_at,:] = np.dot(at_new.get_cell(), rng.float_uniform(0.0, 1.0, (3) ))
            at_new.set_positions(pos[0:i_at+1,:])
            if movement_args['2D']: # zero the Z coordiates in a 2D simulation
                 temp_at=at_new.get_positions()
                 temp_at[:,2]=0.0
                 at_new.set_positions(temp_at)
            energy = eval_energy(at_new)
            if energy < Emax:
                success = True
                break
        if not success:
            exit_error("Failed 10 times to insert atom %d with Emax %f" % (i_at, Emax), 7)
    at.set_positions(pos)
    return energy

def save_snapshot(snapshot_id):
    """
    Save the current walker configurations' as a snapshot in the file ``out_file_prefix.iter.rank.config_file_format``
    """

    if comm is not None:
        comm.barrier() # if parallel, ensure that we are always in sync, so snapshots are always a consistent set

    if ns_args['snapshot_per_parallel_task']:
        rank_id = "%d" % rank
    else:
        rank_id = "ALL"

    if ns_args['snapshot_per_parallel_task'] or rank == 0:
        try:
            snapshot_io = open(ns_args['out_file_prefix']+'snapshot.%s.%s.%s' % (snapshot_id, rank_id, ns_args['config_file_format']), "w")
        except:
            snapshot_io = open(ns_args['out_file_prefix']+'snapshot.%d.%s.%s' % (snapshot_id, rank_id, ns_args['config_file_format']), "w")

        root_walkers_write_t0 = time.time()
        for at in walkers:
            at.info['volume'] = at.get_volume()
            at.info['iter']=snapshot_id
            ase.io.write(snapshot_io, at, parallel=False, format=ns_args['config_file_format'])
        print( "root walkers write time ", time.time() - root_walkers_write_t0)

    if not ns_args['snapshot_per_parallel_task']:
        if comm is not None: # gather other walkers to do I/O on master task
            if rank == 0: # I/O on master task
                for r in range(1,size):
                    remote_walkers_recv_t0 = time.time()
                    remote_walkers = comm.recv(source=r, tag=2)
                    print( "save_snapshot remote walkers recv time ", r, time.time() - remote_walkers_recv_t0)
                    remote_walkers_write_t0 = time.time()
                    for at in remote_walkers:
                        at.info['volume'] = at.get_volume()
                        at.info['iter']=snapshot_id
                        ase.io.write(snapshot_io, at, format=ns_args['config_file_format'])
                    print( "save_snapshot remote walkers write time ", r, time.time() - remote_walkers_write_t0)
            else: # not master task
                comm.send(walkers, dest=0, tag=2)

    if ns_args['snapshot_per_parallel_task'] or rank == 0:
        snapshot_io.close()

def clean_prev_snapshot(iter):
    if iter is not None and ns_args['snapshot_clean']:
        snapshot_file=ns_args['out_file_prefix']+'snapshot.%d.%d.%s' % (iter, rank, ns_args['config_file_format'])
        try:
            os.remove(snapshot_file)
        except:
            print( print_prefix, ": WARNING: Failed to delete '%s'" % snapshot_file)


def do_ns_loop():
    """
    This is the main nested sampling loop, doing the iterations.
    """
    global print_prefix
    global cur_config_ind

    # cant keep atoms fixed and change the simulation cell at the moment
    if (movement_args['keep_atoms_fixed'] > 0
            and movement_args['n_cell_volume_steps'] > 0):
        exit_error("cant keep atoms fixed and change the simulation cell\n", 11)

    nD = 3  # dimensionality of the system
    if movement_args['2D']:  # simulate a 2D system only
        nD = 2
    if rank == 0:
        if energy_io.tell() == 0:  # tell returns the current stream position
            if movement_args['do_velocities']:
                nExtraDOF = 0
            else:
                if movement_args['keep_atoms_fixed'] > 0:
                    nExtraDOF = (n_atoms-movement_args['keep_atoms_fixed'])*nD
                else:
                    nExtraDOF = n_atoms*nD
            energy_io.write("%d %d %d %s %d true\n" % (ns_args['n_walkers'], ns_args['n_cull'], nExtraDOF, movement_args['MC_cell_flat_V_prior'], n_atoms))

    ## print(print_prefix, ": random state ", np.random.get_state())
    ## if rank == 0:
        ## print(print_prefix, ": common random state ", common_random_state)

    if ns_args['debug'] >= 10 and size <= 1:
        for at in walkers:
            at.info['n_walks'] = 0

    for at in walkers:
        at.info['KEmax'] = KEmax
        if movement_args['MC_cell_P'] > 0:
            print(rank, ": initial enthalpy ", at.info['ns_energy'], " PE ", eval_energy_PE(at), " KE ", eval_energy_KE(at), " PV ", eval_energy_PV(at), " mu ", eval_energy_mu(at), " vol ", at.get_volume())
        else:
            print(rank, ": initial enthalpy ", at.info['ns_energy'], " PE ", eval_energy_PE(at), " KE ", eval_energy_KE(at), " mu ", eval_energy_mu(at), " vol ",at.get_volume())
    sys.stdout.flush()

    # stats for purpose of adjusting step size
    walk_stats_adjust = {}
    # stats for purpose of monitoring acceptance rates
    walk_stats_monitor = {}
    zero_stats(walk_stats_adjust, movement_args)
    zero_stats(walk_stats_monitor, movement_args)

    initial_time = time.time()
    prev_time = initial_time
    step_size_setting_duration = 0.0
    total_step_size_setting_duration = 0.0

    Emax_of_step = None
    Emax_save = []
    i_ns_step_save = []
    traj_walker_list = []
    E_dump_list = []
    E_dump_list_times = []

    verbose = False

    # to avoid errors of unassigned values, if in case of a restart the final number of iter is the same as the starting, stop.
    if start_first_iter == ns_args['n_iter']:
        print("WARNING: Increase the n_iter_times_fraction_killed variable in the input if you want NS cycles to be performed.")
        exit_error("starting iteration and the total number of required iterations are the same,hence no NS cycles will be performed\n", 11)

    last_log_X_n = 0.0
    i_range_mod_n_cull = np.array(range(ns_args['n_cull']))
    i_range_plus_1_mod_n_cull = np.mod(np.array(range(ns_args['n_cull']))+1, ns_args['n_cull'])
    log_X_n_term = np.log(ns_args['n_walkers']-i_range_mod_n_cull) - np.log(ns_args['n_walkers']+1-i_range_mod_n_cull)
    log_X_n_term_cumsum = np.cumsum(log_X_n_term)
    log_X_n_term_cumsum_modified = log_X_n_term_cumsum - np.log(ns_args['n_walkers']+1-i_range_plus_1_mod_n_cull)
    log_X_n_term_sum = log_X_n_term_cumsum[-1]
    if ns_args['converge_down_to_T'] > 0:
        converge_down_to_beta = 1.0/(ns_args['kB']*ns_args['converge_down_to_T'])
        log_Z_term_max = np.NINF

    prev_snapshot_iter = None
    pprev_snapshot_iter = None
    last_snapshot_time = time.time()

    # for estimating current temperature from d log Omega / d E
    if ns_args['T_estimate_finite_diff_lag'] > 0:
        log_alpha = np.log(float(ns_args['n_walkers']+1-ns_args['n_cull'])/float(ns_args['n_walkers']+1))
        Emax_history = collections.deque(maxlen=ns_args['T_estimate_finite_diff_lag'])

    if ns_analyzers is not None:
        for (ns_analyzer, ns_analyzer_interval) in ns_analyzers:
            ns_analyzer.analyze(walkers, -1, "NS_loop start")

    # START MAIN LOOP
    i_ns_step = start_first_iter
    while ns_args['n_iter'] < 0 or i_ns_step < ns_args['n_iter']:

        if ns_args['debug'] == -5:
            print(i_ns_step, rank, " ".join(["{:.2f}".format(eval_energy(x))
                                             for x in walkers]))

        check_memory.check_memory("start_ns_main_loop")
        print_prefix = "%d NS %d" % (rank, i_ns_step)

        if ns_args['debug'] >= 4 and ns_args['track_configs']:
            for at in walkers:
                print(print_prefix, "INFO: 10 config_ind ", at.info['config_ind'], " from ", at.info['from_config_ind'], " at ", at.info['config_ind_time'])

        if movement_args['adjust_step_interval'] < 0:
            zero_stats(walk_stats_adjust, movement_args)
        if movement_args['monitor_step_interval'] < 0:
            zero_stats(walk_stats_monitor, movement_args)

        if ns_args['debug'] >= 20:
            print(print_prefix, "%30s" % ": LOOP_TE START 00 ", i_ns_step, ["%.10f" % eval_energy(at) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_PE START 01 ", i_ns_step, ["%.10f" % eval_energy(at, do_KE=False) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_X START 02 ", i_ns_step, ["%.10f" % at.positions[0, 0] for at in walkers])

        # get list of highest energy configs
        (Emax, Vmax, cull_rank, cull_ind) = max_energy(walkers, n_cull)
        Emax_next = Emax[-1]
        if rank == 0 and Emax_of_step is not None and Emax[0] > Emax_of_step:
            print(print_prefix, ": WARNING: energy above Emax ", Emax_of_step, " bad energies: ", Emax[np.where(Emax > Emax_of_step)], cull_rank[np.where(Emax > Emax_of_step)], cull_ind[np.where(Emax > Emax_of_step)])
            # comm.barrier()
            # exit_error("Energy above Emax\n", 5)

        if rank == 0 and (i_ns_step > start_first_iter and Emax_next >= Emax_of_step):
            print("WARNING: Emax not decreasing ", Emax_of_step, Emax_next)
        Emax_of_step = Emax_next

        if ns_args['min_Emax'] is not None and Emax_of_step < ns_args['min_Emax']:
            if rank == 0:
                # if the termination was set by a minimum energy, and it is reached, stop.
                print("Leaving loop because Emax=", Emax_of_step, " < min_Emax =", ns_args['min_Emax'])
            i_ns_step += 1  # add one so outside loop when one is subtracted to get real last iteration it's still correct
            break

        if rank == 0:
            cur_time = time.time()
            output_this_iter = (cur_time > prev_time+60 or i_ns_step == 0 or i_ns_step == ns_args['n_iter'] or (ns_args['n_iter'] > 0 and i_ns_step % max(int(ns_args['n_iter']/1000), 1) == 0))
        else:
            output_this_iter = False

        if ns_args['converge_down_to_T'] > 0:
            # see ns_analyse.py calc_log_a() for math
            log_a = log_X_n_term_sum*i_ns_step + log_X_n_term_cumsum_modified
            # DEBUG if rank == 0:
                # DEBUG for ii in range(len(log_a)):
                    # DEBUG print i_ns_step, "log_a, beta, Es, beta*Es ", log_a[ii], beta, Emax[ii], beta*Emax[ii]
            log_Z_term_max = max(log_Z_term_max, np.amax(log_a - converge_down_to_beta * Emax))
            log_Z_term_last = log_a[-1]-converge_down_to_beta*Emax[-1]
            if output_this_iter:
                print("log_Z_term max ", log_Z_term_max, "last ", log_Z_term_last, "diff ", log_Z_term_max-log_Z_term_last)
            if log_Z_term_last < log_Z_term_max - 10.0:
                if rank == 0:
                    print(print_prefix, "Leaving loop because Z(%f) is converged" % ns_args['converge_down_to_T'])
                i_ns_step += 1  # add one so outside loop when one is subtracted to get real last iteration it's still correct
                break

        if ns_args['T_estimate_finite_diff_lag'] > 0:
            Emax_history.append(Emax_of_step)
        if output_this_iter:
            if ns_args['T_estimate_finite_diff_lag'] > 0 and len(Emax_history) > 1:
                beta_estimate = (len(Emax_history)-1)*log_alpha/(Emax_history[-1]-Emax_history[0])
                T_estimate = 1.0/(ns_args['kB']*beta_estimate)
            else:
                T_estimate = -1
            print(i_ns_step, "Emax_of_step ", Emax_of_step, "T_estimate ", T_estimate, " loop time ", cur_time-prev_time-step_size_setting_duration, " time spent setting step sizes: ", step_size_setting_duration)
            sys.stdout.flush()
            prev_time = cur_time
            step_size_setting_duration = 0.0

        entries_for_this_rank = np.where(cull_rank == rank)[0]
        cull_list = cull_ind[entries_for_this_rank]
        if rank == 0 and ns_args['debug'] >= 4 and len(cull_ind[entries_for_this_rank]) > 0:
            print(print_prefix, "INFO: 20 cull ", cull_ind[entries_for_this_rank], " on ", rank)

        # record Emax walkers energies
        if rank == 0:
            for ii, (E, V) in enumerate(zip(Emax, Vmax)):
                energy_io.write("%d %.50f %.50f %d\n" % (i_ns_step, E, V, ns_args['n_walkers'] - ii))
            energy_io.flush()

            ## Save the energies and corresponding iteration numbers in a list then print(them out only when printing a snapshot)
            #Emax_save.extend(Emax)
            #i_ns_step_save.extend(n_cull*[i_ns_step])
            ## if it is time to print((i.e. at the same iteration when a snapshot is written, or at every iter if no snapshots - for smooth restarts))
            #if ns_args['snapshot_interval'] < 0 or i_ns_step % ns_args['snapshot_interval'] == ns_args['snapshot_interval']-1:
            #    for istep,E in zip(i_ns_step_save,Emax_save):
            #        energy_io.write("%d %.60f\n" % (istep, E))
            #    energy_io.flush()
            #    #empty the save lists, so they are ready for the next bunch of saved energies
            #    Emax_save[:]=[]
            #    i_ns_step_save[:]=[]

        # record Emax walkers configurations
        if cull_list is not None:
            for (i, global_n_offset) in zip(cull_list, entries_for_this_rank):
                if ns_args['debug'] >= 10 and size <= 1:
                    print(print_prefix, "walker killed at age ", walkers[i].info['n_walks'])
                # store culled config in list to be written (when
                # snapshot_interval has passed) every traj_interval steps
                global_n = i_ns_step*n_cull + global_n_offset
                if ns_args['traj_interval'] > 0 and global_n % ns_args['traj_interval'] == ns_args['traj_interval']-1:
                    walker_copy = walkers[i].copy()
                    walker_copy.info['volume'] = walker_copy.get_volume()
                    walker_copy.info['ns_P'] = movement_args['MC_cell_P']
                    walker_copy.info['iter'] = i_ns_step
                    walker_copy.info['config_n_global'] = global_n
                    if walker_copy.has('masses') and walker_copy.has('momenta'):
                        walker_copy.info['ns_KE'] = walker_copy.get_kinetic_energy()

                    traj_walker_list.append(walker_copy)

                # if tracking all configs, save this one that has been culled
                if track_traj_io is not None:
                    at = walkers[i].copy()
                    at.info['culled'] = True
                    ase.io.write(track_traj_io, at, parallel=False, format=ns_args['config_file_format'])

        if ns_args['E_dump_interval'] > 0 and i_ns_step % ns_args['E_dump_interval'] == 0:  # ns_args['E_dump_interval']-1:
            if walkers[0].has('masses') and walkers[0].has('momenta'):
                E_dump_list.append([w.info['ns_energy'] - w.get_kinetic_energy() for w in walkers])
            else:
                E_dump_list.append([w.info['ns_energy'] for w in walkers])
            E_dump_list_times.append(i_ns_step)

        if ns_args['traj_interval'] > 0:
            for at in traj_walker_list:
                ase.io.write(traj_io, at, parallel=False, format=ns_args['config_file_format'])
            traj_io.flush()
            traj_walker_list=[]

        # print(the recorded Emax walkers configurations to output file)
        if (ns_args['snapshot_interval'] < 0 or i_ns_step % ns_args['snapshot_interval'] == ns_args['snapshot_interval']-1 or
            (ns_args['snapshot_seq_pairs'] and i_ns_step > 0 and i_ns_step%ns_args['snapshot_interval'] == 0) ) :
            ##NB if ns_args['traj_interval'] > 0:
                ##NB for at in traj_walker_list:
                    ##NB ase.io.write(traj_io, at, parallel=False, format=ns_args['config_file_format'])
                ##NB traj_io.flush()
                ##NB traj_walker_list=[]
            if ns_args['E_dump_interval'] > 0:
                if comm is not None:
                    E_dump_list_all = np.array(comm.allgather(E_dump_list))
                else:
                    E_dump_list_all = np.array(E_dump_list)
                if rank == 0:
                    for i in range(E_dump_list_all.shape[1]):
                        E_dump_io.write("step %d\n" % E_dump_list_times[i])
                        if len(E_dump_list_all.shape) == 3:
                            np.savetxt(E_dump_io, E_dump_list_all[:,i,:])
                        else:
                            np.savetxt(E_dump_io, E_dump_list_all[i,:])
                    E_dump_io.flush()
                E_dump_list = []
                E_dump_list_all = None
                E_dump_list_times = []

        # calculate how many will be culled on each rank
        n_cull_of_rank = np.array([sum(cull_rank == r) for r in range(size)])

        # label configs to be culled
        status = np.empty((size, n_walkers), np.object_)
        status[:, :] = ''
        for r in range(size):
            status[r, cull_ind[np.where(cull_rank == r)[0]]] = 'c_t'

        if ns_args['debug'] >= 10:
            initial_PE_loc = [eval_energy(at, do_KE=False) for at in walkers]
            initial_E_loc = [eval_energy(at) for at in walkers]
            if comm is not None:
                initial_PE = np.array(comm.allgather(initial_PE_loc)).flatten()
                initial_E = np.array(comm.allgather(initial_E_loc)).flatten()
            else:
                initial_PE = np.array(initial_PE_loc)
                initial_E = np.array(initial_E_loc)
            initial_changed = initial_PE[np.where(status.flatten() == 'c_t')]
            initial_unchanged = initial_PE[np.where(status.flatten() == '')]

        if ns_args['debug'] >= 30:
            for r in range(len(status)):
                print(print_prefix, ": initial status ", r,
                      [s for s in status[r, :]])

        # find load balance by cloning on top of excess maxima
        recv_ind = []
        recv_rank = []
        send_ind = []
        send_rank = []
        cull_inds_to_remove = []

        if n_cull > 1:  # send/recv for fixing load balance
            # CHECK FOR RANDOMNESS ISSUES AND WHICH NODES ARE USED FOR CLONES
            for r in range(size):
                # maybe remote_r should be chosen completely randomly, rather
                # than close to task of extra culled configs
                for dr in np.array(list(zip(np.array(range(1, size)), -np.array(range(1, size))))).flatten():
                    if n_cull_of_rank[r] <= max_n_cull_per_task: # not too many that need to be culled on this rank
                        break
                    # this rank has too many to cull, must receive replacement from another node
                    remote_r = (r+dr) % size
                    if n_cull_of_rank[remote_r] < max_n_cull_per_task: # send from r+dr to r
                        n_transfer = min(n_cull_of_rank[r]-max_n_cull_per_task, max_n_cull_per_task-n_cull_of_rank[remote_r])
                        recv_rank.extend([r]*n_transfer)
                        send_rank.extend([remote_r]*n_transfer)
                        local_ind = np.where(status[r, :] == 'c_t')[0][0:n_transfer]
                        recv_ind.extend(local_ind)
                        remote_ind = np.where(status[remote_r, :] == '')[0][0:n_transfer]
                        send_ind.extend(remote_ind)
                        status[r, local_ind] = 'c_s'
                        status[remote_r, remote_ind] = 'c_t_a'
                        n_cull_of_rank[r] -= n_transfer
                        n_cull_of_rank[remote_r] += n_transfer

        # save local random state, and switch to common one
        rng.switch_to_common()

        # select clones
        for r in range(size):
            list_clone_target = np.where(status[r, :] == 'c_t')[0]
            # assign clones
            n_remaining_clones = len(list_clone_target)
            while n_remaining_clones > 0:
                remote_r = rng.int_uniform(0, size)
                n_avail_remote = sum(status[remote_r, :] == '')
                if n_avail_remote > 0:  # something is available on remote_r
                    # send from random avail walker on remote_r to clone_target on r
                    n_transfer = min(n_remaining_clones, n_avail_remote)

                    # set ranks
                    send_rank.extend([remote_r]*n_transfer)
                    recv_rank.extend([r]*n_transfer)

                    # set indices
                    r_is = []
                    for ii in range(n_transfer):
                        r_i = rng.int_uniform(0, n_walkers)
                        while status[remote_r, r_i] != '':
                            r_i = rng.int_uniform(0, n_walkers)
                        # now r_i should be something with status ''
                        status[remote_r, r_i] = 'c_s'
                        r_is.append(r_i)
                    send_ind.extend(r_is)

                    status[r, list_clone_target[0:n_transfer]] = 'c_t_a'
                    recv_ind.extend(list_clone_target[0:n_transfer])

                    if n_transfer < len(list_clone_target):
                        list_clone_target = list_clone_target[n_transfer:]
                    n_remaining_clones -= n_transfer

        if ns_args['debug'] >= 20:
            print(print_prefix, "%30s" % ": LOOP_TE POST_LOC_CLONE 15 ", i_ns_step, ["%.10f" % eval_energy(at) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_PE POST_LOC_CLONE 16 ", i_ns_step, ["%.10f" % eval_energy(at, do_KE=False) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_X POST_LOC_CLONE 17 ", i_ns_step, ["%.10f" % at.positions[0, 0] for at in walkers])

        # make into numpy arrays so that mathematical operations will work
        send_rank = np.array(send_rank)
        send_ind = np.array(send_ind)
        recv_rank = np.array(recv_rank)
        recv_ind = np.array(recv_ind)

        if ns_args['debug'] >= 10:
            if rank == 0:
                for i in range(len(send_rank)):
                    print(print_prefix, "send from ", send_rank[i], send_ind[i],
                          " to ", recv_rank[i], recv_ind[i])

        # save new common state, and restore to local state
        rng.switch_to_local()

        if n_cull == 1:
            if send_rank[0] == recv_rank[0] and send_rank[0] == rank:  # local copy
                walkers[recv_ind[0]].set_positions(walkers[send_ind[0]].get_positions())
                walkers[recv_ind[0]].set_cell(walkers[send_ind[0]].get_cell())
                if movement_args['do_velocities']:
                    walkers[recv_ind[0]].set_velocities(walkers[send_ind[0]].get_velocities())
                if movement_args['do_GMC']:
                    walkers[recv_ind[0]].arrays['GMC_direction'][:, :] = walkers[send_ind[0]].arrays['GMC_direction']
                if ns_args['n_extra_data'] > 0:
                    walkers[recv_ind[0]].arrays['ns_extra_data'][...] = walkers[send_ind[0]].arrays['ns_extra_data']
                if ns_args['swap_atomic_numbers']:
                    walkers[recv_ind[0]].set_atomic_numbers(walkers[send_ind[0]].get_atomic_numbers())
                    if movement_args['do_velocities']:
                        walkers[recv_ind[0]].set_masses(walkers[send_ind[0]].get_masses())
                if ns_args['track_configs']:
                    walkers[recv_ind[0]].info['config_ind'] = walkers[send_ind[0]].info['config_ind']
                    walkers[recv_ind[0]].info['from_config_ind'] = walkers[send_ind[0]].info['from_config_ind']
                    walkers[recv_ind[0]].info['config_ind_time'] = walkers[send_ind[0]].info['config_ind_time']
                walkers[recv_ind[0]].info['ns_energy'] = eval_energy(walkers[recv_ind[0]])
                if ns_args['debug'] >= 10 and size <= 1:
                    walkers[recv_ind[0]].info['n_walks'] = 0
            else:  # need send/recv
                n_send = 3*(n_atoms + 3)
                if movement_args['do_velocities']:
                    n_send += 3*n_atoms
                if movement_args['do_GMC']:
                    n_send += 3*n_atoms
                if ns_args['n_extra_data'] > 0:
                    n_send += ns_args['n_extra_data']*n_atoms
                if ns_args['swap_atomic_numbers']:
                    n_send += n_atoms  # Z
                    if movement_args['do_velocities']:
                        n_send += n_atoms  # mass
                if ns_args['track_configs']:
                    n_send += 3
                buf = np.zeros(n_send)
                if send_rank[0] == rank:  # only one config is sent/received
                    buf_o = 0
                    buf[buf_o:buf_o+3*n_atoms] = walkers[send_ind[0]].get_positions().reshape((3*n_atoms)); buf_o += 3*n_atoms
                    buf[buf_o:buf_o+3*3] = walkers[send_ind[0]].get_cell().reshape((3*3)); buf_o += 3*3
                    if movement_args['do_velocities']:
                        buf[buf_o:buf_o+3*n_atoms] = walkers[send_ind[0]].get_velocities().reshape((3*n_atoms)); buf_o += 3*n_atoms
                    if movement_args['do_GMC']:
                        buf[buf_o:buf_o+3*n_atoms] = walkers[send_ind[0]].arrays['GMC_direction'].reshape((3*n_atoms)); buf_o += 3*n_atoms
                    if ns_args['n_extra_data'] > 0:
                        buf[buf_o:buf_o+ns_args['n_extra_data']*n_atoms] = walkers[send_ind[0]].arrays['ns_extra_data'].reshape((ns_args['n_extra_data']*n_atoms)); buf_o += ns_args['n_extra_data']*n_atoms
                    if ns_args['swap_atomic_numbers']:
                        buf[buf_o:buf_o+n_atoms] = walkers[send_ind[0]].get_atomic_numbers(); buf_o += n_atoms
                        if movement_args['do_velocities']:
                            buf[buf_o:buf_o+n_atoms] = walkers[send_ind[0]].get_masses(); buf_o += n_atoms
                    if ns_args['track_configs']:
                        buf[buf_o] = walkers[send_ind[0]].info['config_ind']; buf_o += 1
                        buf[buf_o] = walkers[send_ind[0]].info['from_config_ind']; buf_o += 1
                        buf[buf_o] = walkers[send_ind[0]].info['config_ind_time']; buf_o += 1
                    comm.Send([buf,  MPI.DOUBLE], dest=recv_rank[0], tag=100)
                elif recv_rank[0] == rank:
                    comm.Recv([buf, MPI.DOUBLE], source=send_rank[0], tag=100)
                    buf_o = 0
                    walkers[recv_ind[0]].set_positions(buf[buf_o:buf_o+3*n_atoms].reshape((n_atoms, 3))); buf_o += 3*n_atoms
                    walkers[recv_ind[0]].set_cell(buf[buf_o:buf_o+3*3].reshape((3, 3))); buf_o += 3*3
                    if movement_args['do_velocities']:
                        walkers[recv_ind[0]].set_velocities(buf[buf_o:buf_o+3*n_atoms].reshape((n_atoms, 3))); buf_o += 3*n_atoms
                    if movement_args['do_GMC']:
                        walkers[recv_ind[0]].arrays['GMC_direction'][:, :] = buf[buf_o:buf_o+3*n_atoms].reshape((n_atoms, 3)); buf_o += 3*n_atoms
                    if ns_args['n_extra_data'] > 0:
                        walkers[recv_ind[0]].arrays['ns_extra_data'][...] = buf[buf_o:buf_o+3*n_atoms].reshape(walkers[recv_ind[0]].arrays['ns_extra_data'].shape); buf_o += ns_args['n_extra_data']*n_atoms
                    if ns_args['swap_atomic_numbers']:
                        walkers[recv_ind[0]].set_atomic_numbers(buf[buf_o:buf_o+n_atoms].astype(int)); buf_o += n_atoms
                        if movement_args['do_velocities']:
                            walkers[recv_ind[0]].set_masses(buf[buf_o:buf_o+n_atoms]); buf_o += n_atoms
                    if ns_args['track_configs']:
                        walkers[recv_ind[0]].info['config_ind'] = int(buf[buf_o]); buf_o += 1
                        walkers[recv_ind[0]].info['from_config_ind'] = int(buf[buf_o]); buf_o += 1
                        walkers[recv_ind[0]].info['config_ind_time'] = int(buf[buf_o]); buf_o += 1
                    walkers[recv_ind[0]].info['ns_energy'] = eval_energy(walkers[recv_ind[0]])

        else:  # complicated construction of sending/receiving buffers
            # figure out how much is sent per config
            n_data_per_config = 1+3*(n_atoms + 3)
            if movement_args['do_velocities']:
                n_data_per_config += 3*n_atoms
            if movement_args['do_GMC']:
                n_data_per_config += 3*n_atoms
            if ns_args['n_extra_data'] > 0:
                n_data_per_config += ns_args['n_extra_data']*n_atoms
            if ns_args['swap_atomic_numbers']:
                n_send += n_atoms # Z
                if movement_args['do_velocities']:
                    n_send += n_atoms # mass
            if ns_args['track_configs']:
                n_data_per_config += 3

            # figure out send counts
            send_count = [0] * size
            for i in np.where(send_rank == rank)[0]:
                r_recv = recv_rank[i]
                send_count[r_recv] += n_data_per_config

            # figure out send displacements
            send_displ = [0] * size
            send_displ[0] = 0
            for i in range(1,size):
                send_displ[i] = send_displ[i-1] + send_count[i-1]

            # create empty buffer for sending
            send_count_tot = sum(send_count)
            send_data = np.zeros(send_count_tot)

            # copy data to be sent to buffer
            send_displ_t = list(send_displ)
            for i in np.where(send_rank == rank)[0]:
                r_recv = recv_rank[i]
                i_send = send_ind[i]

                data_o = send_displ_t[r_recv]
                send_data[data_o] = walkers[i_send].info['ns_energy']; data_o += 1
                send_data[data_o:data_o+3*n_atoms] = walkers[i_send].get_positions().reshape( (3*n_atoms) ); data_o += 3*n_atoms
                send_data[data_o:data_o+3*3] = walkers[i_send].get_cell().reshape( (3*3) ); data_o += 3*3
                if movement_args['do_velocities']:
                    send_data[data_o:data_o+3*n_atoms] = walkers[i_send].get_velocities().reshape( (3*n_atoms) ); data_o += 3*n_atoms
                if movement_args['do_GMC']:
                    send_data[data_o:data_o+3*n_atoms] = walkers[i_send].arrays['GMC_direction'].reshape( (3*n_atoms) ); data_o += 3*n_atoms
                if ns_args['n_extra_data'] > 0:
                    send_data[data_o:data_o+ns_args['n_extra_data']*n_atoms] = walkers[i_send].arrays['ns_extra_data'].reshape( (ns_args['n_extra_data']*n_atoms) ); data_o += ns_args['n_extra_data']*n_atoms
                if ns_args['swap_atomic_numbers']:
                    send_data[data_o:data_o+n_atoms] = walkers[i_send].get_atomic_numbers(); data_o += n_atoms
                    if movement_args['do_velocities']:
                        send_data[data_o:data_o+n_atoms] = walkers[i_send].get_masses(); data_o += n_atoms
                if ns_args['track_configs']:
                    send_data[data_o] = walkers[i_send].info['config_ind']; data_o += 1
                    send_data[data_o] = walkers[i_send].info['from_config_ind']; data_o += 1
                    send_data[data_o] = walkers[i_send].info['config_ind_time']; data_o += 1
                send_displ_t[r_recv] = data_o

            # figure out recv counts
            recv_count = [0] * size
            for i in np.where(recv_rank == rank)[0]:
                r_send = send_rank[i]
                recv_count[r_send] += n_data_per_config

            # figure out recv displacements
            recv_displ = [0] * size
            recv_displ[0] = 0
            for i in range(1,size):
                recv_displ[i] = recv_displ[i-1] + recv_count[i-1]

            # create empty buffer for receiving
            recv_count_tot = sum(recv_count)
            recv_data = np.zeros(recv_count_tot)

            # do communications
            if comm is not None:
                send_buf = [send_data, send_count, send_displ, MPI.DOUBLE]
                recv_buf = (recv_data, recv_count, recv_displ, MPI.DOUBLE)
                comm.Alltoallv(send_buf, recv_buf)
            else:
                recv_data = send_data.copy()

            # copy data from recv buffer to walkers
            recv_displ_t = list(recv_displ)
            for i in np.where(recv_rank == rank)[0]:
                r_send = send_rank[i]
                i_recv = recv_ind[i]

                data_o = recv_displ_t[r_send]
                walkers[i_recv].info['ns_energy'] = recv_data[data_o]; data_o += 1
                walkers[i_recv].set_positions( recv_data[data_o:data_o+3*n_atoms].reshape( (n_atoms, 3) )); data_o += 3*n_atoms
                walkers[i_recv].set_cell( recv_data[data_o:data_o+3*3].reshape( (3, 3) )); data_o += 3*3
                if movement_args['do_velocities']:
                    walkers[i_recv].set_velocities( recv_data[data_o:data_o+3*n_atoms].reshape( (n_atoms, 3) )); data_o += 3*n_atoms
                if movement_args['do_GMC']:
                    walkers[i_recv].arrays['GMC_direction'][:,:] = recv_data[data_o:data_o+3*n_atoms].reshape( (n_atoms, 3) ); data_o += 3*n_atoms
                if ns_args['n_extra_data'] > 0:
                    walkers[i_recv].arrays['ns_extra_data'][...] = recv_data[data_o:data_o+ns_args['n_extra_data']*n_atoms].reshape( walkers[i_recv].arrays['ns_extra_data'].shape ); data_o += ns_args['n_extra_data']*n_atoms
                if ns_args['swap_atomic_numbers']:
                    walkers[i_recv].set_atomic_numbers(recv_data[data_o:data_o+n_atoms].astype(int)); data_o += n_atoms
                    if movement_args['do_velocities']:
                        walkers[i_recv].set_masses(recv_data[data_o:data_o+n_atoms]); data_o += n_masses
                if ns_args['track_configs']:
                    walkers[i_recv].info['config_ind'] = int(recv_data[data_o]); data_o += 1
                    walkers[i_recv].info['from_config_ind'] = int(recv_data[data_o]); data_o += 1
                    walkers[i_recv].info['config_ind_time'] = int(recv_data[data_o]); data_o += 1
                recv_displ_t[r_send] = data_o

        if ns_args['debug'] >= 20:
            print(print_prefix, "%30s" % ": LOOP_TE POST_CLONE 20 ", i_ns_step, ["%.10f" % eval_energy(at) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_PE POST_CLONE 21 ", i_ns_step, ["%.10f" % eval_energy(at, do_KE=False) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_X POST_CLONE 22 ", i_ns_step, ["%.10f" % at.positions[0, 0] for at in walkers])

        if ns_args['track_configs']:
            # loop over _all_ clone targets and increment cur_config_ind,
            # setting appropriate configs' new config_ind as needed
            for r in range(size):
                clone_walk_ind = np.where(status[r, :] == 'c_t_a')[0]
                for i_at in clone_walk_ind:
                    if r == rank:
                        walkers[i_at].info['from_config_ind'] = walkers[i_at].info['config_ind']
                        walkers[i_at].info['config_ind'] = cur_config_ind
                        walkers[i_at].info['config_ind_time'] = i_ns_step
                    cur_config_ind += 1
        # move cloned walkers

        if (i_ns_step == start_first_iter and movement_args['full_auto_step_sizes']):
            # set initial step sizes. Performed here since this is the first time all the arrays are in place
            conf_pre=walkers[0].copy()
            conf_pre.calc = walkers[0].calc
            move_args_pre=deepcopy(movement_args)
            walk_stats_pre=walk_single_walker(conf_pre, move_args_pre, Emax_of_step, KEmax)
            delta_step_size_setting_duration = full_auto_set_stepsizes(walkers, walk_stats_pre, movement_args, comm, Emax_of_step, KEmax, size)
            total_step_size_setting_duration += delta_step_size_setting_duration
            step_size_setting_duration += delta_step_size_setting_duration
            del(walk_stats_pre)
            del(move_args_pre)
            del(conf_pre)

        sys.stdout.flush()
        # walk clone targets
        if ns_args['debug'] >= 4:
            for i in np.where(status[rank, :] == 'c_s')[0]:
                print(print_prefix, "INFO: 30 clone source ", rank, i)
        clone_walk_ind = np.where(status[rank, :] == 'c_t_a')[0]
        for i_at in clone_walk_ind:
            if ns_args['debug'] >= 4:
                print(print_prefix, "INFO: 40 WALK clone_target ", rank, i_at)
            walk_stats = walk_single_walker(walkers[i_at], movement_args,
                                            Emax_of_step, KEmax)
            walkers[i_at].info['last_walked_iter_clone'] = i_ns_step
            # if tracking all configs, save this one that has been walked
            if track_traj_io is not None:
                walkers[i_at].info['iter'] = i_ns_step
                ase.io.write(track_traj_io, walkers[i_at], format=ns_args['config_file_format'])
            #print("WALK on rank ", rank, "at iteration ", i_ns_step, " walker ", i_at )
            if ns_args['debug'] >= 10 and size <= 1:
                walkers[i_at].info['n_walks'] += movement_args['n_model_calls']
            accumulate_stats(walk_stats_adjust, walk_stats)
            accumulate_stats(walk_stats_monitor, walk_stats)
        sys.stdout.flush()

        if ns_args['debug'] >= 20:
            print(print_prefix, "%30s" % ": LOOP_TE POST_CLONE_WALK 25 ", i_ns_step, ["%.10f" % eval_energy(at) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_PE POST_CLONE_WALK 26 ", i_ns_step, ["%.10f" % eval_energy(at, do_KE=False) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_X POST_CLONE_WALK 27 ", i_ns_step, ["%.10f" % at.positions[0, 0] for at in walkers])

        # check that everything that should have been changed has, and things
        # that shouldn't have, haven't
        if ns_args['debug'] >= 10:
            final_PE_loc = [eval_energy(at, do_KE=False) for at in walkers]
            final_E_loc = [eval_energy(at) for at in walkers]
            if comm is not None:
                final_PE = np.array(comm.allgather(final_PE_loc)).flatten()
                final_E = np.array(comm.allgather(final_E_loc)).flatten()
            else:
                final_PE = final_PE_loc
                final_E = final_E_loc
            if rank == 0:
                final_status = status.flatten()
                for e in initial_unchanged:
                    if e not in final_PE:
                        print("initial_PE ", initial_PE)
                        print("final_PE ", final_PE)
                        print("initial_E ", initial_E)
                        print("final_E ", final_E)
                        print("final_status ", final_status)
                        print("WARNING: energy that should have been unchanged ", e," missing from final energies")
                for e in initial_changed:
                    if e in final_PE:
                        print("initial_PE ", initial_PE)
                        print("final_PE ", final_PE)
                        print("initial_E ", initial_E)
                        print("final_E ", final_E)
                        print("final_status ", final_status)
                        print("WARNING: energy that should have been changed ", e," still there in final energies")

        # walk extras
        if not ns_args['no_extra_walks_at_all']:
                r_i = rng.int_uniform(0, n_walkers)
                # WARNING: this may select walkers for extra walks multiple
                # times, yet never re-walk ones that were walked as clone
                # targets
                while status[rank, r_i] != '' and status[rank, r_i] != 'c_s':
                    r_i = rng.int_uniform(0, n_walkers)
                if ns_args['debug'] >= 4:
                    print(print_prefix, "INFO: 50 WALK extra ", rank, r_i)
                walk_stats = walk_single_walker(walkers[r_i], movement_args,
                                                Emax_of_step, KEmax)
                walkers[r_i].info['last_walked_iter_extra'] = i_ns_step
                # if tracking all configs, save this one that has been walked
                if track_traj_io is not None:
                    walkers[i_at].info['iter'] = i_ns_step
                    ase.io.write(track_traj_io, walkers[i_at],
                                 format=ns_args['config_file_format'])
                # print("WALK EXTRA on rank ", rank, "at iteration ", i_ns_step,
                # " walker ", r_i)
                if ns_args['debug'] >= 10 and size <= 1:
                    #walkers[r_i].info['n_walks'] += movement_args['n_steps'] # LIVIA-this gives error, does not exist
                    walkers[r_i].info['n_walks'] += movement_args['atom_traj_len']
                accumulate_stats(walk_stats_adjust, walk_stats)
                accumulate_stats(walk_stats_monitor, walk_stats)

        monitored_this_step = False
        if movement_args['monitor_step_interval'] != 0 and i_ns_step % abs(movement_args['monitor_step_interval']) == abs(movement_args['monitor_step_interval'])-1:
            adjust_step_sizes(walk_stats_monitor, movement_args, comm, monitor_only=True)
            zero_stats(walk_stats_monitor, movement_args)
            monitored_this_step = True

        if movement_args['adjust_step_interval'] != 0 and i_ns_step % abs(movement_args['adjust_step_interval']) == abs(movement_args['adjust_step_interval'])-1:

            if (not movement_args['full_auto_step_sizes']):
                adjust_step_sizes(walk_stats_adjust, movement_args, comm, do_print_rate=(not monitored_this_step))
            else:
                delta_step_size_setting_duration = full_auto_set_stepsizes(walkers, walk_stats_adjust, movement_args, comm, Emax_of_step, KEmax, size)
                total_step_size_setting_duration += delta_step_size_setting_duration
                step_size_setting_duration += delta_step_size_setting_duration
            zero_stats(walk_stats_adjust, movement_args)

        if ns_args['debug'] >= 20:
            print(print_prefix, "%30s" % ": LOOP_TE END 30 ", i_ns_step, ["%.10f" % eval_energy(at) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_PE END 31 ", i_ns_step, ["%.10f" % eval_energy(at, do_KE=False) for at in walkers])
            print(print_prefix, "%30s" % ": LOOP_X END 32 ", i_ns_step, ["%.10f" % at.positions[0, 0] for at in walkers])

        if ns_args['debug'] >= 30:
            for r in range(len(status)):
                print(print_prefix, ": final status ", r, [s for s in status[r, :]])

        if (rank == 0) and ((ns_args['snapshot_interval'] > 0 and i_ns_step > 0 and i_ns_step % ns_args['snapshot_interval'] == 0) or
                            (ns_args['snapshot_seq_pairs'] and i_ns_step > 1 and i_ns_step%ns_args['snapshot_interval'] == 1) or
                            (ns_args['snapshot_time'] > 0 and time.time()-last_snapshot_time > ns_args['snapshot_time'])):
            do_snapshot = True
        else:
            do_snapshot = False
        if comm is not None:
            do_snapshot = comm.bcast(do_snapshot, root=0)
        if do_snapshot:
            save_snapshot(i_ns_step)
            last_snapshot_time = time.time()
            clean_prev_snapshot(pprev_snapshot_iter)
            pprev_snapshot_iter = prev_snapshot_iter
            prev_snapshot_iter = i_ns_step

        if ns_analyzers is not None:
            for (ns_analyzer, ns_analyzer_interval) in ns_analyzers:
                if ns_analyzer_interval > 0 and (i_ns_step+1)%ns_analyzer_interval == 0:
                    ns_analyzer.analyze(walkers, i_ns_step, "NS_loop %d" % i_ns_step)
        i_ns_step += 1
        ### END OF MAIN LOOP

    # flush remaining traj configs
    for at in traj_walker_list:
        ase.io.write(traj_io, at, parallel=False, format=ns_args['config_file_format'])
    traj_io.flush()
    traj_walker_list = []

    if ns_args['E_dump_interval'] > 0:
        if comm is not None:
            E_dump_list_all = np.array(comm.allgather(E_dump_list))
        else:
            E_dump_list_all = np.array(E_dump_list)
        if rank == 0:
            for i in range(E_dump_list_all.shape[1]):
                E_dump_io.write("step %d\n" % E_dump_list_times[i])
                if len(E_dump_list_all.shape) == 3:
                    np.savetxt(E_dump_io, E_dump_list_all[:,i,:])
                else:
                    np.savetxt(E_dump_io, E_dump_list_all[i,:])
            E_dump_io.flush()

    cur_time = time.time()
    if rank == 0:
        print( "LOOP TIME total ",cur_time-initial_time-total_step_size_setting_duration, " per iter ", (cur_time-initial_time-total_step_size_setting_duration)/(i_ns_step+1))
        print( "TIME SPENT SETTING STEP SIZES total ",total_step_size_setting_duration)

    return i_ns_step-1


def main():
        """ Main function """
        global movement_args
        global ns_args, start_first_iter
        global max_n_cull_per_task
        global size, rank, comm, rng, np, sys, ns_analyzers
        global n_cull, n_walkers, n_walkers_per_task
        global n_extra_walk_per_task
        global do_calc_ASE, do_calc_lammps, do_calc_internal, do_calc_fortran
        global energy_io, traj_io, walkers
        global n_atoms, KEmax, pot
        global MPI, f_MC_MD
        global track_traj_io, cur_config_ind
        global E_dump_io
        global print_prefix
        global Z_list

        import sys

        sys.excepthook = excepthook_mpi_abort

        stacktrace.listen()

        if len(sys.argv) != 1 and len(sys.argv) != 2:
            usage()
            sys.exit(1)

        use_mpi = True
        if len(sys.argv) == 2:
            if sys.argv[1] == "-no_mpi":
                use_mpi = False
            else:
                usage()
                sys.exit(1)

        # initialize mpi
        comm = None
        calculator_comm = None
        rank = 0
        size = 1
        if use_mpi:
            print("INFO: use_mpi true, importing mpi4py module")
            try:
                from mpi4py import MPI
            except:
                sys.stderr.write("Failed to import mpi4py\n")
                sys.exit(10)
            comm = MPI.COMM_WORLD
            calculator_comm = MPI.COMM_SELF

        if use_mpi:
            try:
                rank = comm.Get_rank()  # id of process
                size = comm.Get_size()  # number of processes
            except:
                exit_error("Failed to get rank or size\n", 10)

        if comm is not None:
            print("comm ", comm, " size ", size, " rank ", rank)

        # read inputs on root, then bcast
        if rank == 0:
            lines = sys.stdin.readlines()
            if len(lines) == 0:
                try:
                    infile = open("ns_inputs", "r")
                except:
                    exit_error("Failed to read ns_inputs file\n", 1)
                lines = infile.readlines()
            args = {}
            if rank == 0:
                for line in lines:
                    if re.match("\s*(#.*)?$", line):
                        continue
                    matches = re.match("\s*(\S+)\s*=\s*(.*\S)", line)
                    if matches is None:
                        exit_error("Failed to parse line '%s'" % line, 1)
                    args[matches.group(1)] = matches.group(2)
        else:
            args = None
        if comm is not None:
            args = comm.bcast(args, root=0)  # send args to other processes

#DOC ``main``: parse arguments
        # parse args
        ns_args = {}

        ns_args['check_memory'] = str_to_logical(args.pop('check_memory', 'F'))
        if ns_args['check_memory'] and (comm is None or comm.rank == 0):
            check_memory.active = True

        # convert from strings to actual args
        try:
            ns_args['n_walkers'] = int(args.pop('n_walkers'))
        except:
            exit_error("need number of walkers n_walkers\n", 1)

        ns_args['n_cull'] = int(args.pop('n_cull', 1))

        ns_args['n_iter_times_fraction_killed'] = float(args.pop('n_iter_times_fraction_killed', -1))
        if ns_args['n_iter_times_fraction_killed'] > 0:
            ns_args['n_iter'] = int(round(ns_args['n_iter_times_fraction_killed']/(float(ns_args['n_cull'])/float(ns_args['n_walkers']))))
        else:
            ns_args['n_iter'] = -1
        ns_args['converge_down_to_T'] = float(args.pop('converge_down_to_T', -1))
        if ns_args['n_iter'] <= 0 and ns_args['converge_down_to_T'] <= 0:
            exit_error("need either n_iter_times_fraction_killed or converge_down_to_T", 1)

        ns_args['T_estimate_finite_diff_lag'] = int(args.pop('T_estimate_finite_diff_lag', 1000))+1

        try:
            ns_args['min_Emax'] = float(args.pop('min_Emax'))
        except:
            ns_args['min_Emax'] = None

        ns_args['start_species'] = args.pop('start_species', None)
        ns_args['start_config_file'] = args.pop('start_config_file', None)
        if ns_args['start_species'] is None and ns_args['start_config_file'] is None:
            exit_error("always need start_species or start_config_file, even if restart_file is specified\n",1)
        if ns_args['start_species'] is not None and ns_args['start_config_file'] is not None:
            exit_error("can't specify both start_species and start_config_file\n",1)
        ns_args['restart_file'] = args.pop('restart_file', 'AUTO')

        ns_args['max_volume_per_atom'] = float(args.pop('max_volume_per_atom', 1.0e3))
        ns_args['min_volume_per_atom'] = float(args.pop('min_volume_per_atom', 1.0))

        ns_args['out_file_prefix'] = args.pop('out_file_prefix', '')
        if ns_args['out_file_prefix'] != '':
            ns_args['out_file_prefix'] += '.'
        ns_args['profile'] = int(args.pop('profile', -1))
        ns_args['debug'] = int(args.pop('debug', -1))
        ns_args['snapshot_interval'] = int(args.pop('snapshot_interval', -1))
        ns_args['snapshot_time'] = int(args.pop('snapshot_time', 3600))
        ns_args['snapshot_per_parallel_task'] = str_to_logical(args.pop('snapshot_per_parallel_task', 'T'))
        ns_args['snapshot_seq_pairs'] = str_to_logical(args.pop('snapshot_seq_pairs', "F"))
        ns_args['snapshot_clean'] = str_to_logical(args.pop('snapshot_clean', "T"))
        ns_args['random_initialise_pos'] = str_to_logical(args.pop('random_initialise_pos', "T"))
        ns_args['random_initialise_cell'] = str_to_logical(args.pop('random_initialise_cell', "T"))
        ns_args['LAMMPS_molecular_info'] = str_to_logical(args.pop('LAMMPS_molecular_info', "T"))
        ns_args['initial_walk_N_walks'] = int(args.pop('initial_walk_N_walks', 0))
        ns_args['initial_walk_adjust_interval'] = int(args.pop('initial_walk_adjust_interval', 1))
        ns_args['initial_walk_Emax_offset_per_atom'] = float(args.pop('initial_walk_Emax_offset_per_atom', 1))
        ns_args['initial_walk_only'] = str_to_logical(args.pop('initial_walk_only', 'F'))
        ns_args['traj_interval'] = int(args.pop('traj_interval', 100))
        ns_args['E_dump_interval'] = int(args.pop('E_dump_interval', -1))
        ns_args['delta_random_seed'] = int(args.pop('delta_random_seed', -1))
        ns_args['n_extra_walk_per_task'] = int(args.pop('n_extra_walk_per_task', 0))
        ns_args['random_energy_perturbation'] = float(args.pop('random_energy_perturbation', 1.0e-12))
        ns_args['n_extra_data'] = int(args.pop('n_extra_data', 0))
        ns_args['Z_cell_axis'] = float(args.pop('Z_cell_axis', 10.0))

        # surely there's a cleaner way of doing this?
        try:
            ns_args['start_energy_ceiling_per_atom'] = float(args.pop('start_energy_ceiling_per_atom'))
        except:
            ns_args['start_energy_ceiling_per_atom'] = None
        try:
            ns_args['start_energy_ceiling'] = float(args.pop('start_energy_ceiling'))
        except:
            ns_args['start_energy_ceiling'] = None
        if ns_args['start_energy_ceiling_per_atom'] is not None and ns_args['start_energy_ceiling'] is not None:
            # conflict
            exit_error("got both start_energy_ceiling and start_energy_ceiling_per_atom\n", 1)
        elif ns_args['start_energy_ceiling_per_atom'] is None and ns_args['start_energy_ceiling'] is None:
            # neither specified, use default
            ns_args['start_energy_ceiling_per_atom'] = 1.0e9
        elif ns_args['start_energy_ceiling'] is not None and rank == 0:
            # warn on deprecated feature
            sys.stderr.write("WARNING: got DEPRECATED start_energy_ceiling\n")
        ns_args['random_init_max_n_tries'] = int(args.pop('random_init_max_n_tries', 100))

        ns_args['KEmax_max_T'] = float(args.pop('KEmax_max_T', -1))
        ns_args['kB'] = float(args.pop('kB', 8.6173324e-5))  # eV/K

        # parse energy_calculator
        ns_args['energy_calculator'] = args.pop('energy_calculator', 'fortran')
        do_calc_ASE = False
        do_calc_lammps = False
        do_calc_internal = False
        do_calc_fortran = False
        if ns_args['energy_calculator'] == 'ASE':
            try:
                ns_args['ASE_calc_module'] = args.pop('ASE_calc_module')
            except:
                exit_error("need ASE calculator module name ASE_calc_module\n",1)

            do_calc_ASE=True
        elif ns_args['energy_calculator'] == 'lammps':
            try:
                from lammpslib import LAMMPSlib
            except:
                exit_error("energy_calculator=lammps and failed to import lammpslib module\n", 1)
            do_calc_lammps=True
            try:
                ns_args['LAMMPS_init_cmds'] = args.pop('LAMMPS_init_cmds')
            except:
                exit_error("need LAMMPS initialization commands LAMMPS_init_cmds\n", 1)
            ns_args['LAMMPS_name'] = args.pop('LAMMPS_name', os.environ.get('LAMMPS_name', ''))
            ns_args['LAMMPS_serial'] = str_to_logical(args.pop('LAMMPS_serial', 'T'))
            ns_args['LAMMPS_header'] = args.pop('LAMMPS_header', 'units metal; atom_style atomic; atom_modify map array sort 0 0')
            ns_args['LAMMPS_header_extra'] = args.pop('LAMMPS_header_extra', '')
            ns_args['LAMMPS_atom_types'] = None

            ns_args['LAMMPS_fix_gmc'] = str_to_logical(args.pop('LAMMPS_fix_gmc', "F"))
            LAMMPS_atom_types = args.pop('LAMMPS_atom_types', '')
            if len(LAMMPS_atom_types) > 0:
                if LAMMPS_atom_types == 'TYPE_EQUALS_Z':
                    ns_args['LAMMPS_atom_types'] = LAMMPS_atom_types
                else:
                    ns_args['LAMMPS_atom_types'] = {}
                    for type_pair in [s.strip() for s in LAMMPS_atom_types.split(',')]:
                        f = type_pair.split()
                        ns_args['LAMMPS_atom_types'][f[0]] = int(f[1])
            else:
               exit_error("LAMMPS_atom_types is mandatory if calculator type is LAMMPS\n", 1)
        elif ns_args['energy_calculator'] == 'internal':
            do_calc_internal=True
        elif ns_args['energy_calculator'] == 'fortran':
            import fortranMCMDpy
            do_calc_fortran=True
            try:
                ns_args['FORTRAN_model'] = args.pop('FORTRAN_model')
            except:
                exit_error("need FORTRAN model FORTRAN_model\n", 1)
            ns_args['FORTRAN_model_params'] = args.pop('FORTRAN_model_params', '0')
            f_MC_MD = fortranMCMDpy.fortran_MC_MD(ns_args['FORTRAN_model'])
            params = np.array([float(x) for x in ns_args['FORTRAN_model_params'].split()])
            f_MC_MD.init_model(params)
        else:
            exit_error("energy_calculator=%s unknown\n" % ns_args['energy_calculator'], 3)

        ns_args['no_extra_walks_at_all'] = str_to_logical(args.pop('no_extra_walks_at_all', "F"))

        ns_args['track_configs'] = str_to_logical(args.pop('track_configs', "F"))
        ns_args['track_configs_write'] = str_to_logical(args.pop('track_configs_write', "F"))

        ns_args['config_file_format'] = args.pop('config_file_format', 'extxyz')

        ns_args['rng'] = args.pop('rng', 'numpy')

        ns_args['ns_run_analyzers'] = args.pop('ns_run_analyzers', '')

        if ns_args['rng'] == 'numpy':
            rng = ns_rng.NsRngNumpy(ns_args['delta_random_seed'], comm)
        # elif ns_args['rng'] == 'julia':
        #    import julia
        #    j = julia.Julia()
        #    rng = ns_rng.NsRngJulia(j)
        elif ns_args['rng'] == 'rngstream':
            import rngstream
            rng = ns_rng.NsRngStream(ns_args['delta_random_seed'], comm)
        elif ns_args['rng'] == 'internal':
            rng = ns_rng.NsRngInternal(ns_args['delta_random_seed'], comm)
        else:
            exit_error("rng=%s unknown\n" % ns_args['rng'], 3)

        if do_calc_fortran:
            l_seed = f_MC_MD.seed_size()
            seed = np.array([0] * l_seed, dtype=np.int32)
            for i in range(l_seed):
                # maybe we need a better way of getting max possible int
                seed[i] = rng.int_uniform(1, sys.maxsize)
            f_MC_MD.set_seed(seed)

        ns_args['reproducible'] = str_to_logical(args.pop('reproducible', "F"))
        if ns_args['reproducible']:
            # reset seed after using some random numbers to generate fortran
            # seed, so that fortran and non-fortran have the same seed
            if ns_args['rng'] == 'numpy':
                rng = ns_rng.NsRngNumpy(ns_args['delta_random_seed'], comm)
            elif ns_args['rng'] == 'rngstream':
                rng = ns_rng.NsRngStream(ns_args['delta_random_seed'], comm)
            elif ns_args['rng'] == 'internal':
                rng = ns_rng.NsRngInternal(ns_args['delta_random_seed'], comm)
            else:
                exit_error("rng=%s unknown\n" % ns_args['rng'], 3)

        movement_args = {}

        movement_args['n_model_calls_expected'] = int(args.pop('n_model_calls_expected', 0))
        movement_args['n_model_calls'] = int(args.pop('n_model_calls', 0))
        movement_args['do_good_load_balance'] = str_to_logical(args.pop('do_good_load_balance', "F"))

        #DOC \item process n\_atom\_steps
            #DOC \item If break\_up\_atom\_traj
                #DOC \item n\_atom\_steps\_per\_call = 1
                #DOC \item n\_atom\_steps\_n\_calls = n\_atom\_steps
            #DOC \item else
                #DOC \item n\_atom\_steps\_per\_call = n\_atom\_steps
                #DOC \item n\_atom\_steps\_n\_calls = 1

        movement_args['n_atom_steps'] = int(args.pop('n_atom_steps', 1))
        movement_args['atom_traj_len'] = int(args.pop('atom_traj_len', 8))
        movement_args['atom_traj_len_cost_multiplier'] = int(args.pop('atom_traj_len_cost_multiplier', 1))
        movement_args['break_up_atom_traj'] = str_to_logical(args.pop('break_up_atom_traj', "T"))
        if movement_args['n_atom_steps'] == 0:
            movement_args['n_atom_steps_per_call'] = 0
            movement_args['n_atom_steps_n_calls'] = 0
        else:
            if movement_args['break_up_atom_traj']:
                movement_args['n_atom_steps_per_call'] = 1
                movement_args['n_atom_steps_n_calls'] = movement_args['n_atom_steps']
            else:
                movement_args['n_atom_steps_per_call'] = movement_args['n_atom_steps']
                movement_args['n_atom_steps_n_calls'] = 1

        movement_args['n_cell_volume_steps'] = int(args.pop('n_cell_volume_steps', 1))
        movement_args['n_cell_shear_steps'] = int(args.pop('n_cell_shear_steps', 1))
        movement_args['n_cell_stretch_steps'] = int(args.pop('n_cell_stretch_steps', 1))

        movement_args['n_swap_steps'] = int(args.pop('n_swap_steps', 0))
        movement_args['swap_max_cluster'] = int(args.pop('swap_max_cluster', 1))
        movement_args['swap_r_cut'] = float(args.pop('swap_r_cut', 2.5))
        movement_args['swap_cluster_probability_increment'] = float(args.pop('cluster_probability_increment', 0.75))
        movement_args['swap_velo'] = str_to_logical(args.pop('swap_velo', "F"))
        movement_args['no_swap_velo_fix_mag_alt'] = str_to_logical(args.pop('no_swap_velo_fix_mag_alt', "F"))

        movement_args['n_semi_grand_steps'] = int(args.pop('n_semi_grand_steps', 0))
        if movement_args['n_semi_grand_steps'] > 0:
            try:
                m = re.match('\s*{([^}]*)}\s*$', args.pop('semi_grand_potentials'))
                semi_grand_potentials = m.group(1).split(',')
            except:
                exit_error('Got n_semi_grand_steps > 0 but no valid semi_grand_potentials\n', 3)
            movement_args['semi_grand_potentials'] = {}
            for semi_grand_item in semi_grand_potentials:
                (Z_s, mu_s) = semi_grand_item.strip().split(":")
                movement_args['semi_grand_potentials'][int(Z_s)] = float(mu_s)
        else:
            movement_args['semi_grand_potentials'] = None

        if movement_args['n_model_calls_expected'] < 0:
            movement_args['n_model_calls_expected'] = (movement_args['n_atom_steps']*movement_args['atom_traj_len']*movement_args['atom_traj_len_cost_multiplier'] +
                                                       movement_args['n_cell_volume_steps'] +
                                                       movement_args['n_cell_shear_steps'] +
                                                       movement_args['n_cell_stretch_steps'] +
                                                       movement_args['n_swap_steps'])

        # initialize swap cluster size probabilities
        movement_args['swap_probs'] = np.zeros((movement_args['swap_max_cluster']))
        movement_args['swap_probs'][0] = 1.0
        for i in range(1, movement_args['swap_max_cluster']):
            movement_args['swap_probs'][i] = movement_args['swap_probs'][i-1] * movement_args['swap_cluster_probability_increment']
        movement_args['swap_probs'] /= np.sum(movement_args['swap_probs'])
        for i in range(1, movement_args['swap_max_cluster']):
            movement_args['swap_probs'][i] = movement_args['swap_probs'][i] + movement_args['swap_probs'][i-1]

        if (movement_args['n_model_calls_expected'] <= 0 and
            movement_args['n_model_calls'] <= 0):
            exit_error("Got all of n_model_calls* == 0\n", 3)

        if (movement_args['n_atom_steps'] <= 0 and
            movement_args['n_cell_volume_steps'] <= 0 and
            movement_args['n_cell_shear_steps'] <= 0 and
            movement_args['n_cell_stretch_steps'] <= 0 and
            movement_args['n_swap_steps'] <= 0 and
            movement_args['n_semi_grand_steps'] <= 0):
            exit_error("Got all of n_steps_* == 0\n", 3)

        movement_args['velo_traj_len'] = int(args.pop('velo_traj_len', 8))

        try:
            movement_args['atom_algorithm'] = args.pop('atom_algorithm')
        except:
            exit_error("Failed to read algorithm for atom motion atom_algorithm", 1)
        if movement_args['atom_algorithm'] != 'MC' and movement_args['atom_algorithm'] != 'MD' and movement_args['atom_algorithm'] != 'GMC':
            exit_error("Got unknown atom_algorithm '%s'\n" % movement_args['atom_algorithm'], 3)

        movement_args['MC_atom_velocities'] = str_to_logical(args.pop('MC_atom_velocities', "F"))
        movement_args['MC_atom_velocities_pre_perturb'] = str_to_logical(args.pop('MC_atom_velocities_pre_perturb', "F"))
        movement_args['MC_atom_step_size'] = float(args.pop('MC_atom_step_size', 1.0))
        movement_args['MC_atom_step_size_max'] = float(args.pop('MC_atom_step_size_max', 1.0))
        movement_args['MC_atom_velo_step_size'] = float(args.pop('MC_atom_velo_step_size', 50.0))
        movement_args['MC_atom_velo_step_size_max'] = float(args.pop('MC_atom_velo_step_size_max', 10000.0))
        movement_args['MC_atom_uniform_rv'] = str_to_logical(args.pop('MC_atom_uniform_rv', "F"))
        movement_args['MC_atom_Galilean'] = str_to_logical(args.pop('MC_atom_Galilean', "F"))
        movement_args['GMC_no_reverse'] = str_to_logical(args.pop('GMC_no_reverse', "T"))
        movement_args['do_velocities'] = (movement_args['atom_algorithm'] == 'MD' or movement_args['MC_atom_velocities'])
        # atom_algorithm == GMC is just an alias for atom_algorithm = MC, MC_atom_Galilean = True
        if movement_args['atom_algorithm'] == 'GMC':
            movement_args['atom_algorithm'] = 'MC'
            movement_args['MC_atom_Galilean'] = True
        movement_args['python_MD'] = str_to_logical(args.pop('python_MD', "F"))

        movement_args['MD_atom_velo_pre_perturb'] = str_to_logical(args.pop('MD_atom_velo_pre_perturb', "F"))
        movement_args['MD_atom_velo_post_perturb'] = str_to_logical(args.pop('MD_atom_velo_post_perturb', "T"))
        movement_args['MD_atom_velo_flip_accept'] = str_to_logical(args.pop('MD_atom_velo_flip_accept', "F"))
        movement_args['atom_velo_rej_free_fully_randomize'] = str_to_logical(args.pop('atom_velo_rej_free_fully_randomize', "F"))
        movement_args['atom_velo_rej_free_perturb_angle'] = float(args.pop('atom_velo_rej_free_perturb_angle', 0.3))
        movement_args['MC_atom_velo_walk_rej_free'] = str_to_logical(args.pop('MC_atom_velo_walk_rej_free', "T"))

        movement_args['MD_atom_timestep'] = float(args.pop('MD_atom_timestep', 0.1))
        movement_args['MD_atom_timestep_max'] = float(args.pop('MD_atom_timestep_max', 2.0))
        movement_args['MD_atom_energy_fuzz'] = float(args.pop('MD_atom_energy_fuzz', 1.0e-2))
        movement_args['MD_atom_reject_energy_violation'] = str_to_logical(args.pop('MD_atom_reject_energy_violation', "F"))

        movement_args['MC_cell_P'] = float(args.pop('MC_cell_P', 0.0))
        movement_args['MC_cell_flat_V_prior'] = str_to_logical(args.pop('MC_cell_flat_V_prior', "F"))

        default_value = ns_args['max_volume_per_atom']/20.0  # 5% of maximum allowed volume per atom
        movement_args['MC_cell_volume_per_atom_step_size'] = float(args.pop('MC_cell_volume_per_atom_step_size', default_value))
        movement_args['MC_cell_volume_per_atom_step_size_max'] = float(args.pop('MC_cell_volume_per_atom_step_size_max', 10.0*default_value))  # 50% of maximum allowed volume per atom
        movement_args['MC_cell_volume_per_atom_prob'] = float(args.pop('MC_cell_volume_per_atom_prob', 1.0))
        movement_args['MC_cell_stretch_step_size'] = float(args.pop('MC_cell_stretch_step_size', 0.1))
        movement_args['MC_cell_stretch_step_size_max'] = float(args.pop('MC_cell_stretch_step_size_max', 1.0))
        movement_args['MC_cell_stretch_prob'] = float(args.pop('MC_cell_stretch_prob', 1.0))
        movement_args['MC_cell_shear_step_size'] = float(args.pop('MC_cell_shear_step_size', 0.1))
        movement_args['MC_cell_shear_step_size_max'] = float(args.pop('MC_cell_shear_step_size_max', 1.0))
        movement_args['MC_cell_shear_prob'] = float(args.pop('MC_cell_shear_prob', 1.0))

        movement_args['MC_cell_min_aspect_ratio'] = float(args.pop('MC_cell_min_aspect_ratio', 0.8))
        movement_args['cell_shape_equil_steps'] = int(args.pop('cell_shape_equil_steps', 1000))

        try:
            movement_args['monitor_step_interval'] = float(args.pop('monitor_step_interval'))
        except:
            movement_args['monitor_step_interval_times_fraction_killed'] = float(args.pop('monitor_step_interval_times_fraction_killed', 1))
            movement_args['monitor_step_interval'] = int(round(movement_args['monitor_step_interval_times_fraction_killed']/(float(ns_args['n_cull'])/float(ns_args['n_walkers']))))
        try:
            movement_args['adjust_step_interval'] = float(args.pop('adjust_step_interval'))
        except:
            movement_args['adjust_step_interval_times_fraction_killed'] = float(args.pop('adjust_step_interval_times_fraction_killed', 1))
            movement_args['adjust_step_interval'] = int(round(movement_args['adjust_step_interval_times_fraction_killed']/(float(ns_args['n_cull'])/float(ns_args['n_walkers']))))
#       if movement_args['adjust_step_interval'] < 20:
#           print("WARNING: step size adjustment would be done too often, at every ", movement_args['adjust_step_interval'], " iteration")
#           print("WARNING: adjust_step_interval is increased to 20")
#           movement_args['adjust_step_interval'] = 20
        movement_args['full_auto_step_sizes'] = str_to_logical(args.pop('full_auto_step_sizes', "T"))

        movement_args['MC_adjust_step_factor'] = float(args.pop('MC_adjust_step_factor', 1.5))
        movement_args['MC_adjust_min_rate'] = float(args.pop('MC_adjust_min_rate', 0.25))
        movement_args['MC_adjust_max_rate'] = float(args.pop('MC_adjust_max_rate', 0.75))
        movement_args['GMC_adjust_min_rate'] = float(args.pop('GMC_adjust_min_rate', 0.25))
        movement_args['GMC_adjust_max_rate'] = float(args.pop('GMC_adjust_max_rate', 0.75))
        movement_args['GMC_dir_perturb_angle'] = float(args.pop('GMC_dir_perturb_angle', -1.0))
        movement_args['GMC_dir_perturb_angle_during'] = float(args.pop('GMC_dir_perturb_angle_during', 0.0))
        movement_args['MD_adjust_step_factor'] = float(args.pop('MD_adjust_step_factor', 1.1))
        movement_args['MD_adjust_min_rate'] = float(args.pop('MD_adjust_min_rate', 0.50))
        movement_args['MD_adjust_max_rate'] = float(args.pop('MD_adjust_max_rate', 0.95))

        movement_args['2D'] = str_to_logical(args.pop('2D', "F"))
        movement_args['keep_atoms_fixed'] = int(args.pop('keep_atoms_fixed', 0))
        movement_args['apply_Z_wall'] = str_to_logical(args.pop('apply_Z_wall', "F"))
        if movement_args['apply_Z_wall']:
            # TODO: RY - the wall_dist should be a parameter and should allow user to change its value
            movement_args['wall_dist'] = 10.00 # LIVIA - review this hard coded value 
        else:
            movement_args['wall_dist'] = 0.0

        if len(args) > 0:
            exit_error(str(args)+"\nUnknown arguments read in\n", 2)

        if rank == 0:
            print("ns_args ", pprint.pformat(ns_args))
            print("movement_args ", pprint.pformat(movement_args))

        # initialize in-situ analyzers
        try:
            ns_analyzers=[]
            for analyzer_str in ns_args['ns_run_analyzers'].split(";"):
                try:
                    (analyzer_name, ns_analyzer_interval) = analyzer_str.split()
                    ns_analyzer_interval = int(ns_analyzer_interval)
                    try:
                        analyzer_module = importlib.import_module(analyzer_name.strip())
                    except:
                        analyzer_module = importlib.import_module("ns_run_analyzers."+analyzer_name.strip())
                    ns_analyzers.append((analyzer_module.NSAnalyzer(comm), ns_analyzer_interval))
                    if rank == 0:
                        print("Got NSAnalyzer from", analyzer_name.strip(), "ns_analyzer_interval", ns_analyzer_interval)
                except:
                    if rank == 0:
                        print("Failed to get NSAnalyzer from", analyzer_name.strip())
        except:
            if rank == 0:
                print("no ns_run_analyzers set")
            ns_analyzers = None

        # initialise potential
        if do_calc_ASE:
            pot = importlib.import_module(ns_args['ASE_calc_module']).calc
        elif do_calc_internal or do_calc_fortran:
            pass
        elif do_calc_lammps:
            init_cmds = [s.strip() for s in ns_args['LAMMPS_init_cmds'].split(';')]
            header_cmds = [s.strip() for s in ns_args['LAMMPS_header'].split(';')]
            header_extra_cmds = [s.strip() for s in ns_args['LAMMPS_header_extra'].split(';')]
            if ns_args['LAMMPS_serial']:
                lammps_comm = None
            else:
                lammps_comm = calculator_comm
            if ns_args['debug'] >= 5:
                pot = LAMMPSlib(lmpcmds=init_cmds, atom_types=ns_args['LAMMPS_atom_types'], log_file='lammps.%d.log' % rank, keep_alive=True, lammps_name=ns_args['LAMMPS_name'],
                                lammps_header=header_cmds, lammps_header_extra=header_extra_cmds, comm=lammps_comm, read_molecular_info=ns_args['LAMMPS_molecular_info'])
            else:
                pot = LAMMPSlib(lmpcmds=init_cmds, atom_types=ns_args['LAMMPS_atom_types'], keep_alive=True, lammps_name=ns_args['LAMMPS_name'],
                                lammps_header=header_cmds, lammps_header_extra=header_extra_cmds, comm=lammps_comm, read_molecular_info=ns_args['LAMMPS_molecular_info'])
            if rank == 0:
                print("PRE START_LAMMPS")
                sys.stdout.flush()
            pot.start_lammps()  # so top level things like units will be set
            if rank == 0:
                print("POST START_LAMMPS")
                sys.stdout.flush()
            pot.first_propagate = True
        else:
            exit_error("Need some way of initializing calculator\n", 3)

        # figure out numbers of local walkers
        rank_of_walker = [0] * ns_args['n_walkers']
        if size <= 1:
            n_walkers = ns_args['n_walkers']
        else:
            n_walkers_per_task = ns_args['n_walkers']//size  # using // ensures division gives an integer value
            if n_walkers_per_task*size != ns_args['n_walkers']:
                exit_error("number of walkers %d not divisible by number of MPI processes %d\n" % (ns_args['n_walkers'], size), 5)
            last_walker = 0
            for i_rank in range(size):
                first_walker = last_walker
                last_walker = last_walker + n_walkers_per_task
                if last_walker > ns_args['n_walkers']:
                    last_walker = ns_args['n_walkers']
                if i_rank == rank:
                    n_walkers = last_walker-first_walker
                    my_first_walker = first_walker
                    my_last_walker = last_walker
                if last_walker > first_walker:
                    rank_of_walker[first_walker:last_walker] = [i_rank]*(last_walker-first_walker)

        # figure out number of configs that will be culled on each task
        n_cull = ns_args['n_cull']
        n_extra_walk_per_task = ns_args['n_extra_walk_per_task']
        max_n_cull_per_task = int(n_cull/size)
        if max_n_cull_per_task * size != n_cull:
            max_n_cull_per_task += 1

        # internal model, LJ eps=1, sigma=1, cutoff=3, with PBC cube l=pbc[0,0]
        internal_cutoff = 3.0
        Eshift = internal_cutoff**-12 - internal_cutoff**-6

        set_n_from_expected('n_model_calls')
        if rank == 0:
            print("Using n_model_calls = ", movement_args['n_model_calls'])

        # create list of species, and check for possible problems
        try:
            species_list = ns_args['start_species'].split(',')
        except:
            species_list = []
            if rank == 0:
                init_atoms = ase.io.read(ns_args['start_config_file'], parallel=False)
                atomic_numbers = init_atoms.get_atomic_numbers()
                for Z in set(atomic_numbers):
                    n_of_Z = sum(atomic_numbers == Z)
                    if 'masses' in init_atoms.arrays:
                        mass_of_Z = init_atoms.get_masses()[np.where(atomic_numbers == Z)[0][0]]
                        species_list.append("%d %d %f" % (Z, n_of_Z, mass_of_Z))
                    else:
                        species_list.append("%d %d" % (Z, n_of_Z))
            if comm is not None:
                species_list = comm.bcast(species_list, root=0)


        if do_calc_lammps:
            if not ns_args['LAMMPS_atom_types'] == 'TYPE_EQUALS_Z':
                used_chem_symbols = {ase.data.chemical_symbols[int(species.split()[0])] for species in species_list}
                if not used_chem_symbols == set(ns_args['LAMMPS_atom_types'].keys()):
                    exit_error("species in start_species must correspond to those in LAMMPS_atom_types\n", 1)
        mass_list = []
        Z_list = []
        warned_explicit_mass = False
        for species in species_list:
            species_fields = species.split()
            Z_list.append(int(species_fields[0]))
            if len(species_fields) == 3:
                if not warned_explicit_mass:
                    if rank == 0:
                        sys.stderr.write("WARNING: setting masses explicitly. Not recommended, do only if you're sure it's necessary\n")
                    warned_explicit_mass = True
                type_mass = float(species_fields[2])
                mass_list.append(type_mass)

        if len(mass_list) > 0:
            mass_list = np.array(mass_list)
            if np.any(mass_list != mass_list[0]) and not movement_args['atom_velo_rej_free_fully_randomize']:
                exit_error("ERROR: Masses are not all equal, and atom_velo_rej_free_fully_randomize is false. Refusing to produce incorrect results\n", 1)

        created_temp_restart_file = False
        if ns_args['restart_file'] == "AUTO":
            if rank == 0:
                print("DOING restart_file=AUTO")
                import glob
                sfx = ns_args['config_file_format']

                # try to match .0. or .ALL., but using sloppy regexp
                print("checking snapshots glob", glob.iglob('%ssnapshot.[0-9]*.[0A]*.%s' % (ns_args['out_file_prefix'], sfx)))
                try:
                    newest_snapshot = max(glob.iglob('%ssnapshot.[0-9]*.[0A]*.%s' % (ns_args['out_file_prefix'], sfx)), key=os.path.getmtime)
                    print("restarting from ", newest_snapshot)
                    if re.search('\.ALL\.%s$' % sfx, newest_snapshot) is not None:  # latest snapshot is a combined one for all nodes
                        print("snapshot is combined for all nodes")
                        ns_args['restart_file'] = newest_snapshot
                    else:
                        snapshot_root = re.sub('\.0\.%s' % sfx, "", newest_snapshot)
                        restart_file = snapshot_root+".ALL."+sfx
                        print("creating combined snapshot file", restart_file, "from", [tf for tf in glob.iglob('%s.[0-9]*.%s' % (snapshot_root, sfx))])
                        created_temp_restart_file = True
                        with open(restart_file, "w") as snapshot_out:
                            for snapshot_file in glob.iglob('%s.[0-9]*.%s' % (snapshot_root, sfx)):
                                with open(snapshot_file, "r") as snapshot_in:
                                    for line in snapshot_in:
                                        snapshot_out.write(line)
                        ns_args['restart_file'] = restart_file
                except:
                    print("no snapshot files found")
                    ns_args['restart_file'] = ''

            if comm is not None:
                ns_args['restart_file'] = comm.bcast(ns_args['restart_file'], root=0)
        sys.stdout.flush()

        # set up walkers
        walkers = []
        if ns_args['restart_file'] == '':  # start from scratch
            start_first_iter = 0
            # create initial config
            if rank == 0:
                if ns_args['start_config_file'] is not None:
                    init_atoms = ase.io.read(ns_args['start_config_file'], parallel=False)
                    if not 'masses' in init_atoms.arrays:
                        init_atoms.set_masses([1.0] * len(init_atoms))
                else:
                    # create atoms structs from a list of atomic numbers and numbers of atoms
                    # always create it slightly smaller than the max to avoid numerical instability with nearly identical volumes
                    if movement_args['2D']:
                        lc = 0.999*(ns_args['max_volume_per_atom']/ns_args['Z_cell_axis'])**(1.0/2.0)
                        init_atoms = ase.Atoms(cell=(lc, lc, ns_args['Z_cell_axis']), pbc=(1, 1, 1))
                    else:
                        lc = 0.999*ns_args['max_volume_per_atom']**(1.0/3.0)
                        init_atoms = ase.Atoms(cell=(lc, lc, lc), pbc=(1, 1, 1))

                    for species in species_list:
                        species_fields = species.split()
                        type_Z = int(species_fields[0])
                        type_n = int(species_fields[1])
                        if len(species_fields) == 2:
                            init_atoms += ase.Atoms([type_Z] * type_n, masses=[1.0] * type_n)
                        elif len(species_fields) == 3:
                            type_mass = float(species_fields[2])
                            init_atoms += ase.Atoms([type_Z] * type_n, masses=[type_mass] * type_n)
                        else:
                            exit_error("Each entry in start_species must include atomic number, multiplicity, and optionally mass", 5)

                    init_atoms.set_cell(init_atoms.get_cell()*float(len(init_atoms))**(1.0/3.0), scale_atoms=True)

                ase.io.write(sys.stdout, init_atoms, parallel=False, format=ns_args['config_file_format'])
                # ase.io.write(sys.stdout, init_atoms, format=ns_args['config_file_format'])
            else:  # rank != 0
                init_atoms = None

            # bcast atoms created on rank == 0
            if comm is not None:
                init_atoms = comm.bcast(init_atoms, root=0)

            # create extra data arrays if needed
            if ns_args['n_extra_data'] > 0:
                init_atoms.arrays['ns_extra_data'] = np.zeros( (len(init_atoms), ns_args['n_extra_data']) )

            # clone initial config into array of walkers
            for i_walker in range(n_walkers):
                walkers.append(init_atoms.copy())

            # set up data structures to track configs as they are cloned and evolve
            if ns_args['track_configs']:
                if comm is None:
                    config_ind = 0
                else:
                    config_ind = comm.rank*n_walkers
            for at in walkers:
                at.set_velocities(np.zeros((len(walkers[0]), 3)))
                if ns_args['track_configs']:
                    at.info['config_ind'] = config_ind
                    at.info['from_config_ind'] = -1
                    at.info['config_ind_time'] = -1
                    config_ind += 1
                if do_calc_ASE or do_calc_lammps:
                    at.calc = pot

            if ns_args['track_configs']:
                if comm is not None:
                    cur_config_ind = comm.size*n_walkers
                else:
                    cur_config_ind = n_walkers

            # V should have prob distrib p(V) = V^N.
            # Using transformation rule p(y) = p(x) |dx/dy|, with p(y) = y^N and p(x) = 1,
            #       one gets dx/dy = y^N
            #                x = y^{N+1}
            #                y = x^{1/(N+1)}
            if ns_args['start_energy_ceiling_per_atom'] is not None:
                ns_args['start_energy_ceiling'] = ns_args['start_energy_ceiling_per_atom'] * len(init_atoms)
            ns_args['start_energy_ceiling'] += movement_args['MC_cell_P']*ns_args['max_volume_per_atom']*len(init_atoms)
            # initial positions are just random, up to an energy ceiling
            for (i_at, at) in enumerate(walkers):
                # randomize cell if P is set, both volume (from appropriate distribution) and shape (from uniform distribution with min aspect ratio limit)

                if ns_args['random_initialise_cell']:
                    energy = float('nan')
                    if movement_args['MC_cell_P'] > 0.0:
                        if movement_args['2D']:
                            lc = (len(at)*ns_args['max_volume_per_atom']/ns_args['Z_cell_axis']*rng.float_uniform(0.0,1.0)**(1.0/float(len(at)+1)))**(1.0/2.0)
                            temp_cell = np.identity(3) * lc
                            temp_cell[2,2] = ns_args['Z_cell_axis']
                            at.set_cell( temp_cell )
                            do_cell_shape_walk(at, movement_args)

                        else:
                            lc = (len(at)*ns_args['max_volume_per_atom']*rng.float_uniform(0.0,1.0)**(1.0/float(len(at)+1)))**(1.0/3.0)
                            at.set_cell( np.identity(3) * lc )
                            do_cell_shape_walk(at, movement_args)

                if ns_args['random_initialise_pos']:
                    # random initial positions
                    energy = float('nan')
                    n_try = 0
                    while n_try < ns_args['random_init_max_n_tries'] and (math.isnan(energy) or energy > ns_args['start_energy_ceiling']):
                        at.set_scaled_positions( rng.float_uniform(0.0, 1.0, (len(at), 3) ) )
                        if movement_args['2D']:  # zero the Z coordiates in a 2D simulation
                            temp_at=at.get_positions()
                            temp_at[:,2]=0.0
                            at.set_positions(temp_at)
                        energy = eval_energy(at)
                        n_try += 1
                    if math.isnan(energy) or energy > ns_args['start_energy_ceiling']:
                        sys.stderr.write("WARNING: rank %d failed to generate initial config by random positions under max energy %f in %d tries\n" % (rank, ns_args['start_energy_ceiling'], ns_args['random_init_max_n_tries']))

                    # try FORTRAN config initializer
                    n_try = 0
                    if do_calc_fortran:
                        while n_try < ns_args['random_init_max_n_tries'] and (math.isnan(energy) or energy > ns_args['start_energy_ceiling']):
                            f_MC_MD.init_config(at, ns_args['start_energy_ceiling']-movement_args['MC_cell_P']*ns_args['max_volume_per_atom']*len(init_atoms))
                            energy = eval_energy(at)
                            n_try += 1
                        if math.isnan(energy) or energy > ns_args['start_energy_ceiling']:
                            sys.stderr.write("WARNING: rank %d failed to generate initial config by fortran config initializer under max energy %f in %d tries\n" % (rank, ns_args['start_energy_ceiling'], ns_args['random_init_max_n_tries']))

                    # try python config initializer
                    n_try = 0
                    while n_try < ns_args['random_init_max_n_tries'] and (math.isnan(energy) or energy > ns_args['start_energy_ceiling']):
                        energy = additive_init_config(at, ns_args['start_energy_ceiling'])
                        n_try += 1

                    # quit if failed to generate acceptable config
                    if math.isnan(energy) or energy > ns_args['start_energy_ceiling']:
                        exit_error("Rank %d failed to generate initial config by random, fortran, or (atom by atom addition) python initializer under max energy %f in %d tries each\n" % (rank, ns_args['start_energy_ceiling'], ns_args['random_init_max_n_tries']), 4)

                energy = eval_energy(at)
                at.info['ns_energy'] = rand_perturb_energy(energy, ns_args['random_energy_perturbation'])
                at.info['volume'] = at.get_volume()

            # Done initialising atomic positions. Now initialise momenta

            # set KEmax from P and Vmax
            if (movement_args['do_velocities']):
                if ns_args['KEmax_max_T'] > 0.0:
                    KEmax = 1.5*len(walkers[0])*ns_args['kB']*ns_args['KEmax_max_T']
                elif movement_args['MC_cell_P'] > 0.0:
                    KEmax = 1.5*movement_args['MC_cell_P']*len(walkers[0])*ns_args['max_volume_per_atom']
                else:
                    exit_error("do_velocities is set, but neither KEmax_max_T nor MC_cell_P are > 0, so no heuristic for setting KEmax",4)
                for at in walkers:
                    at.info['KEmax']=KEmax
            else:
                KEmax = -1.0

            # set initial velocities, rejection free
            if movement_args['do_velocities']:
                for at in walkers:
                    rej_free_perturb_velo(at, None, KEmax) # adds KE to at.info['ns_energy']

            # swap atomic numbers if doing semi-grand canonical ensemble
            ns_args['swap_atomic_numbers'] = (movement_args['n_semi_grand_steps'] > 0)

        else:  # set up walkers with a restart
            # swap atomic numbers if doing semi-grand canonical ensemble
            ns_args['swap_atomic_numbers'] = (
                    movement_args['n_semi_grand_steps'] > 0)
            at_list = ase.io.read(ns_args['restart_file'], index=":")  # LIVIA
            for r in range(size):
                if rank == r:
                    walkers = at_list[r*n_walkers:(r+1)*n_walkers]  # TODO: RBW – split walkers on different processes? maybe we need to set things up (energies?) before splitting?
                    print(rank, r, walkers)
            for at in walkers:
                if np.any(at.get_atomic_numbers()
                          != walkers[0].get_atomic_numbers()):
                    ns_args['swap_atomic_numbers'] = True

            # broadcast swap_atomic_numbers in case it was overridden to True
            # by presence of configurations with different atomic number lists
            if comm is not None:
                ns_args['swap_atomic_numbers'] = comm.bcast(
                    ns_args['swap_atomic_numbers'], root=0)

            if ns_args['track_configs']:
                if comm is not None:
                    cur_config_ind = comm.size*n_walkers
                else:
                    cur_config_ind = n_walkers
                if 'config_ind' in walkers[0].info:
                    max_config_ind = max(
                        [at.info['config_ind'] for at in walkers])
                    if comm is not None:
                        comm.allreduce(max_config_ind, op=MPI.MAX)
                    cur_config_ind = max_config_ind+1

            # Calculate ns_energy for manual surface configurations from restart
            if movement_args['keep_atoms_fixed'] > 0:
                # print("RBW: calc ns_energy for man surf configs from restart") # debug
                for (i_at, at) in enumerate(walkers):
                    if do_calc_ASE or do_calc_lammps:
                        at.set_calculator(pot)
                    energy = eval_energy(at)
                    at.info['ns_energy'] = rand_perturb_energy(
                        energy, ns_args['random_energy_perturbation'])
                    at.info['volume'] = at.get_volume()

            if movement_args['do_velocities']:
                KEmax = walkers[0].info['KEmax']
            else:
                KEmax = -1.0

            # Set init velocities for manual surface configurations from restart
            # RBW checked this and seed-averaged initial total energies are ~
            # the same, which means velocities *seem* to be set properly
            # LBP: only perturb if this is a "fake" restart to start a surface sim.
            #       proper restart does not need perturbation as that messes up KE
            if movement_args['keep_atoms_fixed'] > 0 and walkers[0].info['iter'] < 1:
                # print("RBW: set init vels for man surf configs from restart") # debug
                if movement_args['do_velocities']:
                    for at in walkers:
                        # TODO: RBW – make sure surface restart files contain KEmax calculated from number of free atoms, otherwise, their kinetic energy can be inefficiently high
                        energy = eval_energy(at)
                        rej_free_perturb_velo(at, None, KEmax)
                        # adds KE to at.info['ns_energy']
                        energy = eval_energy(at)

            for at in walkers:
                if ns_args['n_extra_data'] > 0 and (not 'ns_extra_data' in at.arrays or at.arrays['ns_extra_data'].size/len(at) != ns_args['n_extra_data']):
                    at.arrays['ns_extra_data'] = np.zeros( (len(at), ns_args['n_extra_data']) )
                if do_calc_ASE or do_calc_lammps:
                    at.calc = pot
                at.info['ns_energy'] = rand_perturb_energy(eval_energy(at), ns_args['random_energy_perturbation'])

                if 'iter' in at.info:
                    start_first_iter = at.info['iter'] + 1
                else:
                    print("ERROR: no iteration number information was found "
                          "in the restart file")
                    exit_error("no iteration number information was found in "
                               "the restart file\n", 5)

                key_found = False
                for key in at.info:
                    # check if 'volume=' info is present in the file used for
                    # restart
                    if key == 'volume':
                        movement_args['MC_cell_volume_per_atom_step_size'] = (
                                at.info['volume']/10.0/len(at))
                        key_found = True
                if not key_found:
                    print( "WARNING: no volume information was found in the restart file. If volume changes will be done, the starting stepsize will be the default")
                    
        sys.stdout.flush()

        # add GMC direction if needed
        movement_args['do_GMC'] = ((movement_args['atom_algorithm'] == 'MC') and (movement_args['MC_atom_Galilean']))
        if movement_args['do_GMC']:
            for at in walkers:
                if 'GMC_direction' not in at.arrays:
                    at.arrays['GMC_direction'] = np.zeros((len(at), 3))

        # scale initial MC_atom_step_size by max_vol^(1/3)
        max_lc = (ns_args['max_volume_per_atom']*len(walkers[0]))**(1.0/3.0)
        movement_args['MC_atom_step_size'] *= max_lc
        movement_args['MC_atom_step_size_max'] *= max_lc
        # scale MC_cell_shear_step_size by max_vol^1.0)
        movement_args['MC_cell_shear_step_size'] *= max_lc
        movement_args['MC_cell_shear_step_size_max'] *= max_lc

        if ns_analyzers is not None:
            for (ns_analyzer, ns_analyzer_interval) in ns_analyzers:
                if ns_analyzer_interval < 0:
                    ns_analyzer.analyze(walkers, -1, "initial_walk start")
        # do initial walks if needed
        if ns_args['initial_walk_N_walks'] > 0 and ns_args['restart_file'] == '':
            if rank == 0:
                print("doing initial_walks", ns_args['initial_walk_N_walks'])
                t0 = time.time()
            (Emax, Vmax, cull_rank, cull_ind) = max_energy(walkers, 1)
            # WARNING: this assumes that all walkers have same numbers of atoms
            Emax = Emax[0] + ns_args['initial_walk_Emax_offset_per_atom']*len(walkers[0])

            # do full walks, not shortened walks that account for NS algorithm
            # parallelization
            save_n_model_calls = movement_args['n_model_calls']
            movement_args['n_model_calls'] = movement_args['n_model_calls_expected']

            walk_stats_adjust = {}
            zero_stats(walk_stats_adjust, movement_args)
            for i_initial_walk in range(ns_args['initial_walk_N_walks']):
                if rank == 0:
                    print("initial walk start iter ", i_initial_walk, "time",
                          time.time()-t0)
                print_prefix = "%d initial_walk %d" % (rank, i_initial_walk)

                if i_initial_walk > 0 and (i_initial_walk-1) % ns_args['initial_walk_adjust_interval'] == 0:  # first adjust is after first walk
                    # could be done before first walk if full_auto_set_stepsize
                    # didn't need walk_stats
                    full_auto_set_stepsizes(walkers, walk_stats_adjust,
                                            movement_args, comm, Emax, -1, size)
                    walk_stats_adjust = {}
                    zero_stats(walk_stats_adjust, movement_args)

                for (i_at, at) in enumerate(walkers):
                    print_prefix = "%d initial_walk %d at %d" % (rank, i_initial_walk, i_at)
                    if ns_args['profile'] == rank:
                        import cProfile
                        pr = cProfile.Profile()
                        walk_stats = pr.runcall(walk_single_walker, at=at, movement_args=movement_args, Emax=Emax, KEmax=-1)
                        pr.dump_stats(ns_args['out_file_prefix']+'initial_walk.profile.stats')
                    else:
                        walk_stats = walk_single_walker(at, movement_args, Emax, -1)
                    accumulate_stats(walk_stats_adjust, walk_stats)

                if ns_analyzers is not None:
                    for (ns_analyzer, ns_analyzer_interval) in ns_analyzers:
                        if ns_analyzer_interval < 0 and (ns_analyzer_interval == -1 or (i_initial_walk+1) % (-ns_analyzer_interval) == 0):
                            ns_analyzer.analyze(walkers, -1, "initial_walk %d" % i_initial_walk)

                if ns_args['snapshot_interval'] > 0 and (i_initial_walk+1) % ns_args['snapshot_interval'] == 0:
                    save_snapshot(i_initial_walk-ns_args['initial_walk_N_walks'])

            # restore walk lengths for rest of NS run
            movement_args['n_model_calls'] = save_n_model_calls

            (KEmax, _, _, _) = max_energy(walkers, 1, kinetic_only=True)

            if ns_args['initial_walk_only']:
                if comm is not None:
                    MPI.Finalize()
                sys.exit(0)

        sys.stdout.flush()
        n_atoms = len(walkers[0])
        # do NS

        # open the file where the trajectory will be printed
        if ns_args['restart_file'] == '':  # start from scratch, so if this file exists, overwrite it
            traj_io = open(ns_args['out_file_prefix']+'traj.%d.%s' % (rank, ns_args['config_file_format']), "w")
            if ns_args['track_configs'] and ns_args['track_configs_write']:
                track_traj_io = open(ns_args['out_file_prefix']+'track_traj.%d.%s' % (rank, ns_args['config_file_format']), "w")
            else:
                track_traj_io = None
        else:  # restart, so the existing file should be appended
            # concatenate existing traj file to before restart
            print(rank, "truncating traj file to start_first_iter",
                  start_first_iter)
            mode = "r+"
            if movement_args['keep_atoms_fixed'] > 0:
                mode = "a+"
            with open(ns_args['out_file_prefix']+'traj.%d.%s' % (rank, ns_args['config_file_format']), mode) as f:
                prev_pos = None
                # loop this way with "while True" and "f.readline()" because directly looping over f does not 
                # set position reported by f.tell() to end of line
                # make sure there's somplace to truncate to if traj file has only one config and it's already too late
                prev_pos = 0
                prev_prev_pos = 0
                while True:
                    l = f.readline()
                    if not l:
                        break
                    m = re.search(r"\biter=(\d+)\b", l)
                    if m is not None:
                        cur_iter = int(m.group(1))
                        if cur_iter >= start_first_iter:
                            # rewind back to two lines back, before start of this config
                            f.seek(prev_prev_pos)
                            f.truncate()
                            break
                    prev_prev_pos = prev_pos
                    prev_pos = f.tell()
            traj_io = open(ns_args['out_file_prefix']+'traj.%d.%s' % (rank, ns_args['config_file_format']), "a")
            if ns_args['track_configs'] and ns_args['track_configs_write']:
                track_traj_io = open(ns_args['out_file_prefix']+'track_traj.%d.%s' % (rank, ns_args['config_file_format']), "a")
            else:
                track_traj_io = None

            # Read the existing traj file and look for the point where we restart from. Truncate the rest.
            # This part is not used because the ASE.io.read takes soooo long, that it makes a restart impossible.
            #traj_io = open(ns_args['out_file_prefix']+'traj.%d.%s' % (rank, ns_args['config_file_format']), "r+")
            #i = 0
            #while True:
            #    at=(ase.io.read(traj_io, format=ns_args['config_file_format'],index=i))
            #   print("ASE.io.read trajectory", rank, i, at.info['iter'])
            #    if at.info['iter'] >= start_first_iter:
            #         at=(ase.io.read(traj_io, format=ns_args['config_file_format'],index=i-1))
            #         traj_io.truncate()
            #         break
            #    i += 1

        sys.stdout.flush()
        if ns_args['E_dump_interval'] > 0 and rank == 0:
            E_dump_io = open(ns_args['out_file_prefix']+'E_dump', "w")
            E_dump_io.write("n_walkers %d\n" % ns_args['n_walkers'])
        else:
            E_dump_io = None

        # open the file where the energies will be printed
        if rank == 0:
            if ns_args['restart_file'] == '':  # start from scratch, so if this file exists, overwrite it
                energy_io = open(ns_args['out_file_prefix']+'energies', 'w')
            else:  # restart, so the existing file should be appended
                try:
                    energy_io = open(ns_args['out_file_prefix']+'energies', 'r+')
                    tmp_iter = 0
                    line = energy_io.readline()  # read the first line of nwalker,ncull..etc information
                    i = 0
                    while True:  # we do create an infinite loop here :(
                        line = energy_io.readline()            # read lines one by one
                        if not line:                           # something went wrong, exit the infinit loop
                            print("WARNING: end of .energies file reached without finding the iteration number", start_first_iter)
                            break
                        i = i + 1
                        if i % 10000 == 0:
                            print(rank, "reading .energies file line %d" % i)
                        if i % n_cull == 0:                    # if this is n_cull-th line, examine the stored iteration
                            tmp_split = line.split()
                            tmp_iter = int(tmp_split[0])       # tmp_iter contains the iteration number of the line as an integer number
                        if tmp_iter == start_first_iter - 1:   # if this is the iteration same as in the snapshot,
                            print(rank, "truncating energy file at line ", i)
                            energy_io.truncate()                #delete the rest of the file, as we are restarting from here
                            break
                except FileNotFoundError:
                    print("WARNING: got restart file, but no corresponding energies file, so creating new one from scratch")
                    energy_io = open(ns_args['out_file_prefix']+'energies', 'w')

        sys.stdout.flush()

        if ns_args['profile'] == rank:
            print("if")
            import cProfile
            pr = cProfile.Profile()
            final_iter = pr.runcall(do_ns_loop)
            pr.dump_stats(ns_args['out_file_prefix']+'profile.stats')
        else:
            print("else")
            final_iter = do_ns_loop()

        # cleanup post loop
        save_snapshot(final_iter)  # this is the final configuration

        for at in walkers:
            print(rank, ": final energy ", at.info['ns_energy'])

        if rank == 0:
            energy_io.close()
        traj_io.close()
        if track_traj_io is not None:
            track_traj_io.close()
        if E_dump_io is not None:
            E_dump_io.close()

        if comm is not None:
            MPI.Finalize()
        sys.exit(0)


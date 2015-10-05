import re, math, time, os
import numpy as np, ase, ase.io
import ns_rng
import stacktrace
from itertools import izip

def usage():
    """ Print help to the standard output about the usage of the code and input parameters. The current list of parameters is the following:

.. glossary::  

    ``max_volume_per_atom=float``
       | Maximum volume per atom allowed during the run.
       | default: 1.0e3

    ``start_species=int int [ float ] [, int int [ float ] ... ]``
       | MANDATORY
       | Atomic number; multiplicity; mass (amu). Info repeated for each species, separated by commas, mass is optional. Mutually exclusive with restart_*, one is required

    ``restart_file=path_to_file``
       | File for restart configs. Mutually exclusive with start_*, one is required. The file should contain the state of the walkers to continue from along with the restart iteration number. Normally such a file can be the concatenated snapshot files. 

    ``n_walkers=int``
       | MANDATORY
       | Total number of walkers, i.e. the size of the live set used. It has to be a multiple of the number of processors used.

    ``n_cull=int``
       | Number of walkers to kill at each NS iteration.
       | deafult: 1

    ``n_extra_walk_per_task=int``
       | default: 0

    ``n_iter_per_walker=int``
       | MANDATORY
       | Number of nested sampling iteration cycles performed per walker. Thus the total number of iterations will be ``n_walkers`` * ``n_iter_per_walker``

    ``min_Emax=float``
       | Termination condition based on Emax.
       | No default.

    ``out_file_prefix=str``
       | String used as prefix for the different output files. 
       | (None)

    ``energy_calculator= ( quip | lammps | internal | fortran)``
       | Energy calculator.
       | default: fortran

    ``n_extra_data=int``
       | Amount of extra data per atom to pass around.
       | default: 0

    ``KEmax_max_T=float`` 
       | Maximum temperature for estimating KEmax if *P* = 0, i.e. fixed V ensemble (will be multiplied by kB ~= 8.6e-5 eV/K)
       | default: 1.0e5

    ``start_energy_ceiling=float``
       | Maximum potential energy for initial configurations.  P*Vmax is added to this automatically in case of NpT runs.
       | default: 1.0e9

    ``n_steps_expected=int``
       | (0, one of these is required)

    ``n_steps=int``
       | (0, one of these is required)

    ``n_steps_unblocked_atom=int``
       | (0, one of these is required)

    ``n_steps_unblocked_cell_volume=int``
       | (0, one of these is required)

    ``n_steps_unblocked_cell_shear=int``
       | (0, one of these is required)

    ``n_steps_unblocked_cell_stretch=int``
       | (0, one of these is required)

    ``n_steps_unblocked_swap=int``
       | (0, one of these is required)

    ``atom_n_substeps=int`` 
       | (10, number of steps (MC sweeps or MD time steps in each segement))

    ``cell_n_substeps=int`` 
       | (1, number of MC sweeps in each segement)

    ``swap_n_substeps=int`` 
       | (0, number of atom swaps in each segement)

    ``velo_n_substeps=int`` 
       | (0, number of MC sweeps in each velocity MC segement)

    ``random_energy_perturbation=float`` 
       | (1.0e-12)

    ``atom_algorithm=[MC | MD]``
       | (MANDATORY)
       | Use either Monte Carlo or Molecular dynamics to explore.
    ``MC_atom_velocities=[T | F]`` 
       | (F, supported only for energy_calculator=fortran)
    ``MC_atom_velocities_pre_perturb=[T | F]``
       | (F, Perturb velocities (rejection free) before MC + velocities walk)
    ``MC_atom_step_size=float`` 
       | (0.1, in units of (max_volume_per_atom * N_atoms)^(1/3) 
    ``MC_atom_step_size_max=float`` 
       | (0.5, in units of (max_volume_per_atom * N_atoms)^(1/3) 
    ``MC_atom_uniform_rv=[T | F]`` 
       | default: F
    ``atom_velo_random_order_perturb=[T | F]`` 
       | Perturb velocities as part of main loop steps
       | default: F
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
     | default: 0.5 (ASE time units)
    ``MD_atom_energy_fuzz=float``
     |  (1.0e-2. Tolerance for rejecting non-energy conserving trajectories, as fraction of KE)
    ``MD_atom_reject_energy_violation=[ T | F ]``
     |  (F, use energy conservation violation (exceeding MD_atom_energy_fuzz * KE) to reject MD trajectories)
    ``atom_velo_rej_free_perturb_randomize=[T | F]``
     |  (F. If true, randomize velocities completely rather than actually perturbing.
    ``atom_velo_rej_free_perturb_angle=float``
     |  (0.3. Max angle in radians for random rotations.)
    ``MC_atom_velo_step_size=float``
     |  (0.1)
    ``MC_atom_velo_step_size_max=float``
     |  (1.0)
    ``MC_atom_velo_walk_rej_free=[T | F]``
     |  (T. If true, use rejection free algorithm for MC_atom_walk
    ``MC_cell_P=float``
     | Pressure value to be used. (Note: the unit of pressure depends on the energy calculator and the potential model used)
     | default: 0.0 
    ``MC_cell_volume_per_atom_step_size=float``
     | Initial volume stepsize for volume change.
     | Default: 5% of the maximum allowed volume
    ``MC_cell_volume_per_atom_step_size_max=float``
     |  (5000.0)
    ``MC_cell_volume_per_atom_prob=float``
     |  (1.0)
    ``MC_cell_stretch_step_size=float``
     |  (0.01)
    ``MC_cell_stretch_step_size_max=float``
     |  (0.05)
    ``MC_cell_stretch_prob=float``
     |  (1.0)
    ``MC_cell_shear_step_size=float``
     |  (0.01)
    ``MC_cell_shear_step_size_max=float``
     |  (0.05)
    ``MC_cell_shear_prob=float``
     |  (1.0)
    ``MC_cell_min_aspect_ratio=float``
     |  (0.9)
    ``cell_shape_equil_steps=int``
     |  (100)
    ``adjust_step_interval_per_walker=float`` 
     | Multipled by number of walkers to get actual interval in iterations, negative for only using last iteration, 0 for no adjust.
     | default: 0.25
    ``MC_adjust_step_factor=float``
     |  default: 1.5
    ``MC_adjust_min_rate=float``
     |  default: 0.25
    ``MC_adjust_max_rate=float``
     |  default: 0.75
    ``MD_adjust_step_factor=float``
     |  default: 1.5
    ``MD_adjust_min_rate=float``
     |  default: 0.95
    ``MD_adjust_max_rate=float``
     |  default: 1.00
    ``QUIP_pot_args=str``
     |  (MANDATORY if energy_calculator=quip)
    ``QUIP_pot_params_file=str``
     |  (MANDATORY if energy_calculator=quip)
    ``FORTRAN_model=str``
     |  (MANDATORY if energy_calculator=fortran)
    ``LAMMPS_init_cmds=str``
     |  (MANDATORY if energy_calculator=lammps)
    ``LAMMPS_name=str``
     |  ('', arch name for lammps shared object file)
    ``LAMMPS_header=str``
     |  (lammpslib.py value default, override lammpslib.py header commands for energy_calculator=lammps)
    ``LAMMPS_header_extra=str``
     |  ('', extra lammpslib.py header commands for energy_calculator=lammps)
    ``config_file_format=str``
     |  (extxyz)
    ``rng=( numpy | internal | rngstream )``
     |  (numpy)
    ``profile=rank_to_profile``
     |  (-1)
    ``2D=[ T | F ]``
     |  (F, unsupported)
    ``debug=debug_level``
     |  (0, <= 0 for no debugging tests/prints)
    ``snapshot_interval=int``
     |  Iteration interval at which a snapshot is created: every process prints out its current walkers in extended xyz format. If it is set <=0, no snapshots will be printed except the final positions at the end of the nested sampling run. Note that when new snapshots are printed, the previous set is deleted. The snapshot files are convenient source to see how the sampling progresses, but these are also the basis to restart a sampling! When using restart, the walkers will be read from these files.
     |  default: 1000
    ``traj_interval=int``
     |  Iteration interval at which the currently culled configuration(s) is/are printed to the trajectory output, in the set format. If it is set <=0, no trajectory files will be printed at all. Useful option for larger runs as the trajectory files can become huge. 
     |  default: 1 
    ``random_seed=seed_shift``
     |  (-1, < 0 for seed from /dev/urandom)
    ``no_extra_walks_at_all=[ T | F ]``
     | (F)
    ``track_configs=[ T | F ]``
     | Track configrations across all walks/clones
     | (F)

    """
    sys.stderr.write("Usage: %s [ -no_mpi ] < input\n" % sys.argv[0])
    sys.stderr.write("input:\n")
    sys.stderr.write("max_volume_per_atom=float (1e3)\n")
    sys.stderr.write("start_species=int int [ float ] [, int int [ float ] ... ] (MANDATORY, atomic_number multiplicity mass (amu). Info repeated for each species, separated by commas, mass is optional. Mutually exclusive with restart_*, one is required\n")
    sys.stderr.write("restart_file=path_to_file (file for restart configs. Mutually exclusive with start_*, one is required)\n")
    sys.stderr.write("n_walkers=int (MANDATORY)\n")
    sys.stderr.write("n_cull=int (1, number of walkers to kill at each NS iteration)\n")
    sys.stderr.write("n_extra_walk_per_task=int (0)\n")
    sys.stderr.write("n_iter_per_walker=int (MANDATORY)\n")
    sys.stderr.write("min_Emax=float (None.  Termination condition based on Emax)\n")
    sys.stderr.write("out_file_prefix=str (None)\n")
    sys.stderr.write("energy_calculator= ( quip | lammps | internal | fortran) (fortran)\n")
    sys.stderr.write("n_extra_data=int (0, amount of extra data per atom to pass around)\n")
    sys.stderr.write("\n")
    sys.stderr.write("KEmax_max_T=float (1e5, maximum temperature for estimating KEmax if P == 0, i.e. fixed V ensemble (will be multiplied by kB ~= 8.6e-5 eV/K))\n")
    sys.stderr.write("start_energy_ceiling=float (1.0e9, max potential energy for initial configs.  P*Vmax is added to this automatically)\n")
    sys.stderr.write("\n")
    sys.stderr.write("n_steps_expected=int (0, one of these is required)\n")
    sys.stderr.write("n_steps=int (0, one of these is required)\n")
    sys.stderr.write("n_steps_unblocked_atom=int (0, one of these is required)\n")
    sys.stderr.write("n_steps_unblocked_cell_volume=int (0, one of these is required)\n")
    sys.stderr.write("n_steps_unblocked_cell_shear=int (0, one of these is required)\n")
    sys.stderr.write("n_steps_unblocked_cell_stretch=int (0, one of these is required)\n")
    sys.stderr.write("n_steps_unblocked_swap=int (0, one of these is required)\n")
    sys.stderr.write("\n")
    sys.stderr.write("atom_n_substeps=int (10, number of steps (MC sweeps or MD time steps in each segement))\n")
    sys.stderr.write("cell_n_substeps=int (1, number of MC sweeps in each segement)\n")
    sys.stderr.write("swap_n_substeps=int (0, number of atom swaps in each segement)\n")
    sys.stderr.write("velo_n_substeps=int (0, number of MC sweeps in each velocity MC segement)\n")
    sys.stderr.write("\n")
    sys.stderr.write("random_energy_perturbation=float (1.0e-12)\n")
    sys.stderr.write("atom_algorithm=[MC | MD] (MANDATORY)\n")
    sys.stderr.write("\n")
    sys.stderr.write("MC_atom_velocities=[T | F] (F, supported only for energy_calculator=fortran)\n")
    sys.stderr.write("MC_atom_velocities_pre_perturb=[T | F] (F, Perturb velocities (rejection free) before MC + velocities walk)\n")
    sys.stderr.write("MC_atom_step_size=float (0.1, in units of (max_volume_per_atom * N_atoms)^(1/3) )\n")
    sys.stderr.write("MC_atom_step_size_max=float (0.5, in units of (max_volume_per_atom * N_atoms)^(1/3) )\n")
    sys.stderr.write("MC_atom_uniform_rv=[T | F] (F)\n")
    sys.stderr.write("\n")
    sys.stderr.write("atom_velo_random_order_perturb=[T | F] (F. Perturb velocities as part of main loop steps\n")
    sys.stderr.write("\n")
    sys.stderr.write("MD_atom_velo_pre_perturb=[T | F] (F. Perturb velocities before MD trajectory\n")
    sys.stderr.write("MD_atom_velo_post_perturb=[T | F] (T. Perturb velocities after MD trajectory\n")
    sys.stderr.write("MD_atom_velo_flip_accept=[T | F] (F)\n")
    sys.stderr.write("MD_atom_timestep=float (0.1, in ASE time units)\n")
    sys.stderr.write("MD_atom_timestep_max=float (0.5, in ASE time units)\n")
    sys.stderr.write("MD_atom_energy_fuzz=float (1.0e-2. Tolerance for rejecting non-energy conserving trajectories, as fraction of KE)\n")
    sys.stderr.write("MD_atom_reject_energy_violation=[ T | F ] (F, use energy conservation violation (exceeding MD_atom_energy_fuzz * KE) to reject MD trajectories)\n")
    sys.stderr.write("\n")
    sys.stderr.write("atom_velo_rej_free_perturb_randomize=[T | F] (F. If true, randomize velocities completely rather than actually perturbing.\n")
    sys.stderr.write("atom_velo_rej_free_perturb_angle=float (0.3. Max angle in radians for random rotations.)\n")
    sys.stderr.write("MC_atom_velo_step_size=float (0.1)\n")
    sys.stderr.write("MC_atom_velo_step_size_max=float (1.0)\n")
    sys.stderr.write("MC_atom_velo_walk_rej_free=[T | F] (T. If true, use rejection free algorithm for MC_atom_walk\n")
    sys.stderr.write("\n")
    sys.stderr.write("\n")
    sys.stderr.write("MC_cell_P=float (0.0)\n")
    sys.stderr.write("MC_cell_volume_per_atom_step_size=float (100.0)\n")
    sys.stderr.write("MC_cell_volume_per_atom_step_size_max=float (5000.0)\n")
    sys.stderr.write("MC_cell_volume_per_atom_prob=float (1.0)\n")
    sys.stderr.write("MC_cell_stretch_step_size=float (0.01)\n")
    sys.stderr.write("MC_cell_stretch_step_size_max=float (0.05)\n")
    sys.stderr.write("MC_cell_stretch_prob=float (1.0)\n")
    sys.stderr.write("MC_cell_shear_step_size=float (0.01)\n")
    sys.stderr.write("MC_cell_shear_step_size_max=float (0.05)\n")
    sys.stderr.write("MC_cell_shear_prob=float (1.0)\n")
    sys.stderr.write("MC_cell_min_aspect_ratio=float (0.9)\n")
    sys.stderr.write("cell_shape_equil_steps=int (100)\n")
    sys.stderr.write("\n")
    sys.stderr.write("adjust_step_interval_per_walker=float (0.25, multipled by number of walkers to get actual interval in iterations, negative for only using last iteration, 0 for no adjust)\n")
    sys.stderr.write("MC_adjust_step_factor=float (1.5)\n")
    sys.stderr.write("MC_adjust_min_rate=float (0.25)\n")
    sys.stderr.write("MC_adjust_max_rate=float (0.75)\n")
    sys.stderr.write("MD_adjust_step_factor=float (1.5)\n")
    sys.stderr.write("MD_adjust_min_rate=float (0.95)\n")
    sys.stderr.write("MD_adjust_max_rate=float (1.00)\n")
    sys.stderr.write("\n")
    sys.stderr.write("QUIP_pot_args=str (MANDATORY if energy_calculator=quip)\n")
    sys.stderr.write("QUIP_pot_params_file=str (MANDATORY if energy_calculator=quip)\n")
    sys.stderr.write("FORTRAN_model=str (MANDATORY if energy_calculator=fortran)\n")
    sys.stderr.write("LAMMPS_init_cmds=str (MANDATORY if energy_calculator=lammps)\n")
    sys.stderr.write("LAMMPS_name=str ('', arch name for lammps shared object file)\n")
    sys.stderr.write("LAMMPS_header=str (lammpslib.py value default, override lammpslib.py header commands for energy_calculator=lammps)\n")
    sys.stderr.write("LAMMPS_header_extra=str ('', extra lammpslib.py header commands for energy_calculator=lammps)\n")
    sys.stderr.write("LAMMPS_atom_types=symbol int [, symbol int ] ... ('', mapping from atomic symbols to type numbers for LAMMPS ASE interface)\n")
    sys.stderr.write("\n")
    sys.stderr.write("config_file_format=str (extxyz)\n") # julia
    sys.stderr.write("rng=( numpy | internal | rngstream ) (numpy)\n") # julia
    sys.stderr.write("profile=rank_to_profile (-1)\n")
    sys.stderr.write("2D=[ T | F ] (F, unsupported)\n")
    sys.stderr.write("debug=debug_level (0, <= 0 for no debugging tests/prints)\n")
    sys.stderr.write("snapshot_interval=int (1000, <=0 for no snapshots except final positions)\n")
    sys.stderr.write("traj_interval=int (1, <=0 for no trajectory)\n")
    sys.stderr.write("random_seed=seed_shift (-1, < 0 for seed from /dev/urandom)\n")
    sys.stderr.write("no_extra_walks_at_all=[ T | F ] (F)\n")
    sys.stderr.write("track_configs=[ T | F ] (F)\n")

def exit_error(message, stat):
    sys.stderr.write(message)
    try:
	comm.Abort(stat)
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

def eval_energy(at, do_PE=True, do_KE=True, do_PV=True):
    """Calls appropriate functions to calcaulet the potential energy, kinetic energy and the p*V term. 
    """
    # potential
    if do_PE:
	if do_calc_quip or do_calc_lammps:
	    if do_calc_lammps:
		#NB maybe not do this every time?
		at.wrap()
	    energy = at.get_potential_energy()
	elif do_calc_internal:
	    energy = energy_internal(at)
	elif do_calc_fortran:
	    energy = f_MC_MD.eval_energy(at)
	else:
	    sys.stderr.write("No way to eval_energy()\n", 5)
    else:
	energy = 0.0

    # kinetic
    if do_KE and at.has('momenta') and at.has('masses'):
	energy += at.get_kinetic_energy()

    if do_PV:
	energy += movement_args['MC_cell_P']*at.get_volume()

    return energy

def propagate_NVE_quippy(at, dt, n_steps):
    old_velo = at.get_velocities()
    if old_velo is not None:
	if not hasattr(at, 'velo'):
	    at.add_property('velo', 0.0, n_cols=3)
	at.velo[:,:] = old_velo.transpose()/(ase.units.Ang/ase.units.fs)
    ds=quippy.DynamicalSystem(at)

    # dt is being converted from ASE units to fs
    if ns_args['debug'] >= 10:
	ds.run(pot, dt=dt/ase.units.fs, n_steps=n_steps, summary_interval=1, write_interval=0, save_interval=0)
    else:
	ds.run(pot, dt=dt/ase.units.fs, n_steps=n_steps, summary_interval=0, write_interval=0, save_interval=0)

    at.set_velocities(at.velo.transpose()*(ase.units.Ang/ase.units.fs))

def propagate_NVE_lammps(at, dt, n_steps):
    if pot.first_MD:
	pot.lmp.command('fix 1 all nve')
	pot.first_MD=False

    # NB maybe not do this every time?
    at.wrap()

    # dt being passed in ASE units
    pot.propagate(at, properties=['energy','forces'],system_changes=['positions'], n_steps=n_steps, dt=dt)

def velo_rv_mag(n):
    if movement_args['2D']:
	unit_rv[:,2] = 0.0
	nDOF = 2.0
    else:
	nDOF = 3.0
    # In 3D rv_mag should have prob distrib p(r) = r^(3N-1).
    # Using transformation rule p(y) = p(x) |dx/dy|, with p(y) = y^{3N-1} and p(x) = 1,
    #       one gets dx/dy = y^{3N-1}
    #                x = y^{3N}
    #                y = x^{1/3N}
    return rng.float_uniform(0.0,1.0)**(1.0/(nDOF*n))

def velo_unit_rv(n):
    unit_rv = rng.normal( 1.0, (n, 3) )
    unit_rv /= np.linalg.norm(unit_rv)
    return unit_rv

def gen_random_velo(at, KEmax, unit_rv=None):
    if unit_rv is None:
	unit_rv = velo_unit_rv(len(at))
    rv_mag = velo_rv_mag(len(at))

    # from Baldock thesis Eq. 11.10
    #     p^{**} = r \mathbf{S} \hat{\mathbf{r}}
    # and 11.11
    #     S_{ij} = \delta_{ij} (2 m_i [ E_{lim} - U(q))])^{1/2}
    # p_i = r (2 m_i)^{1/2} (E-U)^{1/2} \hat{r}_i
    # v_i = p_i / m_i = r (2/m)^{1/2} (E-U)^{1/2} \hat{r}_i

    masses = at.get_masses()
    velocities = rv_mag * np.sqrt(2.0/np.array([masses,]*3).transpose()) * np.sqrt(KEmax) * unit_rv

    return velocities

def pairwise(iterable):
    a = iter(iterable)
    return izip(a, a)

def rej_free_perturb_velo(at, Emax, KEmax, rotate=True):
#DOC
#DOC \vspace*{\baselineskip}
#DOC perturb\_velo
#DOC \begin{itemize}

    if not at.has('momenta') or not at.has('masses'):
	return

    KEmax_use = None
    if Emax is not None:
	initial_KE = eval_energy(at, do_PE=False, do_PV=False)
	KEmax_use = Emax - (at.info['ns_energy'] - initial_KE)
    if KEmax is not None:
	if Emax is None:
	    KEmax_use = KEmax
	else:
	    KEmax_use = min(KEmax_use, KEmax)
    if KEmax_use is None:
	exit_error("rej_free_perturb_velo() called with Emax and KEmax both None\n", 9)

    orig_KE = at.get_kinetic_energy( )
    #DOC \item if atom\_velo\_randomize, pick random velocities consistent with Emax
    if movement_args['atom_velo_rej_free_perturb_randomize']:
	# randomize completely
	at.set_velocities(gen_random_velo(at, KEmax_use))
    #DOC \item else perturb velocities
    #DOC \begin{itemize}
    else: # perturb
	velo = at.get_velocities()
	velo_mag = np.linalg.norm(velo)
	#DOC \item if current velocity=0, can't rescale, so pick random velocities consistent with Emax
	if velo_mag == 0.0:
	    at.set_velocities(gen_random_velo(at, KEmax_use))
	#DOC \item else, pick new random magnitude consistent with Emax, random rotation of current direction with angle uniform in $\pm$ atom\_velo\_perturb\_angle
	else:
	    # pick new random magnitude - count on dimensionality to make change small
	    # WARNING: check this for variable masses
	    masses_2D = at.get_masses().reshape( (len(at),1) )
	    momenta = gen_random_velo(at, KEmax_use, velo/velo_mag) * masses_2D

	    if rotate:
		# apply random rotations
		indices = [ (int(i/3), i%3) for i in range(3*len(at)) ]
		rng.shuffle_in_place(indices)
		for ((ind_1_a,ind_1_c), (ind_2_a,ind_2_c)) in pairwise(indices):
		    ang = rng.float_uniform(-movement_args['atom_velo_rej_free_perturb_angle'],movement_args['atom_velo_rej_free_perturb_angle'])
		    c_ang = np.cos(ang)
		    s_ang = np.sin(ang)
		    v_1 = momenta[ind_1_a,ind_1_c] * c_ang + momenta[ind_2_a,ind_2_c] * s_ang
		    v_2 = -momenta[ind_1_a,ind_1_c] * s_ang + momenta[ind_2_a,ind_2_c] * c_ang
		    momenta[ind_1_a,ind_1_c] = v_1
		    momenta[ind_2_a,ind_2_c] = v_2

	    at.set_velocities(momenta / masses_2D)

    new_KE = at.get_kinetic_energy()
    # rej_free_perturb_velo expects at.info['ns_energy'] to be set correctly initially
    at.info['ns_energy'] += new_KE-orig_KE
    #DOC \end{itemize}
#DOC \end{itemize}

def do_MC_atom_velo_walk(at, movement_args, Emax, KEmax):
    if movement_args['MC_atom_velo_walk_rej_free']:
	rej_free_perturb_velo(at, Emax, KEmax)
	return {}

    n_steps = movement_args['velo_n_substeps']
    step_size = movement_args['MC_atom_velo_step_size']

    initial_KE = eval_energy(at, do_PE=False, do_PV=False)
    KEmax_use = Emax - (at.info['ns_energy'] - initial_KE)
    if KEmax is not None and KEmax < KEmax_use:
	KEmax_use = KEmax

    if do_calc_fortran:
	(n_accept, final_KE) = f_MC_MD.MC_atom_walk_velo(at, n_steps, step_size, KEmax_use)
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
	    for i_at in at_list:
		d_KE = 0.5*masses[i_at]*(np.sum((velocities[i_at,:]+d_vel[i_at,:])**2) - np.sum(velocities[i_at,:]**2))
		if KE + d_KE < KEmax_use:
		    velocities[i_at,:] += d_vel[i_at,:]
		    KE += d_KE
		    n_accept += 1
	at.set_velocities(velocities)
	at.info['ns_energy'] += KE-initial_KE

    return {'MC_atom_velo' : (n_steps*len(at), n_accept)}

def do_MD_atom_walk(at, movement_args, Emax, KEmax):
#DOC
#DOC \vspace*{\baselineskip}
#DOC do\_MD\_atom\_walk
#DOC \begin{itemize}

    """ perform MD walk on the configuration """
    orig_E = at.info['ns_energy']
    if orig_E >= Emax:
	print print_prefix, ": WARNING: orig_E =",orig_E," >= Emax =",Emax

    #DOC \item if MD\_atom\_velo\_pre\_perturb, call perturb\_velo() for magnitude and rotation
    if movement_args['MD_atom_velo_pre_perturb']:
	do_MC_atom_velo_walk(at, movement_args, Emax, KEmax)

    pre_MD_pos = at.get_positions()
    pre_MD_velo = at.get_velocities()
    if ns_args['n_extra_data'] > 0:
	pre_MD_extra_data = at.arrays['ns_extra_data'].copy()

    pre_MD_E = at.info['ns_energy']

    #DOC \item propagate in time atom\_n\_substeps time steps of length MD\_atom\_timestep
    if do_calc_quip:
	propagate_NVE_quippy(at, dt=movement_args['MD_atom_timestep'], n_steps=movement_args['atom_n_substeps'])
	final_E = eval_energy(at)
    elif do_calc_lammps:
	propagate_NVE_lammps(at, dt=movement_args['MD_atom_timestep'], n_steps=movement_args['atom_n_substeps'])
	final_E = eval_energy(at)
    elif do_calc_fortran:
	final_E = f_MC_MD.MD_atom_NVE_walk(at, n_steps=movement_args['atom_n_substeps'], timestep=movement_args['MD_atom_timestep'], debug=ns_args['debug'])
	final_E += eval_energy(at,do_PE=False, do_KE=False)
    else:
	exit_error("Need some non-quippy, non-fortran, non-lammps way of doing MD\n",3)

    reject_fuzz = False
    final_KE = eval_energy(at, do_PE=False, do_PV=False)
    #DOC \item optionally accept/reject entire move on E deviating by less than energy fuzz times kinetic energy
    if abs(final_E-pre_MD_E) > movement_args['MD_atom_energy_fuzz']*final_KE:
	if movement_args['MD_atom_reject_energy_violation']:
	    reject_fuzz = True
	else:
	    print print_prefix, ": WARNING: MD energy deviation > fuzz*final_KE. Pre-MD, post-MD, difference, final_KE ", pre_MD_E, final_E, final_E-pre_MD_E, final_KE

    #DOC \item accept/reject entire move on E $<$ Emax
    reject_Emax = (final_E >= Emax)
    reject_KEmax = (KEmax is not None and final_KE >= KEmax)

    #DOC \item if reject
    #DOC \begin{itemize}
    if reject_fuzz or reject_Emax or reject_KEmax: # reject
	#DOC \item set positions, velocities, energy back to value before perturbation (maybe should be after?)
	# print print_prefix, ": WARNING: reject MD traj Emax ", Emax, " initial E ", orig_E, " velo perturbed E ", pre_MD_E, " final E ",final_E, " KEmax ", KEmax, " KE ", final_KE
	at.set_positions(pre_MD_pos)
	if movement_args['MD_atom_velo_flip_accept']:
	    at.set_velocities(pre_MD_velo)
	else:
	    at.set_velocities(-pre_MD_velo)
	if ns_args['n_extra_data'] > 0:
	    at.arrays['ns_extra_data'][...] = pre_MD_extra_data
	at.info['ns_energy'] = pre_MD_E
	n_accept = 0
    #DOC \end{itemize}
    #DOC \item else
    else: # accept
    #DOC \begin{itemize}
	#DOC \item flip velocities
	# remember to reverse velocities on acceptance to preserve detailed balance, since velocity is (sometimes) being perturbed, not completely randomized
	if movement_args['MD_atom_velo_flip_accept']:
	    at.set_velocities(-at.get_velocities()) # is there a faster way of doing this in ASE?  Can you do at.velocities?
	at.info['ns_energy'] = final_E
	n_accept = 1
    #DOC \end{itemize}

    #DOC \item if MD\_atom\_velo\_post\_perturb, call perturb\_velo() for magnitude and rotation
    if movement_args['MD_atom_velo_post_perturb']:
	do_MC_atom_velo_walk(at, movement_args, Emax, KEmax)

    return {'MD_atom' : (1, n_accept) }
#DOC \end{itemize}

def do_MC_atom_walk(at, movement_args, Emax, KEmax):
#DOC
#DOC \vspace*{\baselineskip}
#DOC do\_MC\_atom\_walk
#DOC \begin{itemize}

    n_steps = movement_args['atom_n_substeps']
    step_size = movement_args['MC_atom_step_size']
    step_size_velo = movement_args['MC_atom_velo_step_size']
    n_accept=0
    n_accept_velo = None

    #DOC \item if MC\_atom\_velocities and MC\_atom\_velocities\_pre\_perturb, perturb velocities, magnitude and and rotation
    if movement_args['MC_atom_velocities'] and movement_args['MC_atom_velocities_pre_perturb']:
	do_MC_atom_velo_walk(at, movement_args, Emax, KEmax)

    #DOC \item loop atom\_n\_substeps times
    #DOC \begin{itemize}
    if do_calc_fortran:
	if movement_args['MC_atom_velocities']:
	    (n_accept, n_accept_velo, final_E) = f_MC_MD.MC_atom_walk(at, n_steps, step_size, Emax-eval_energy(at, do_PE=False, do_KE=False), KEmax, step_size_velo)
	    at.info['ns_energy'] = final_E + eval_energy(at, do_PE=False, do_KE=False)
	else:
	    (n_accept, final_E) = f_MC_MD.MC_atom_walk(at, n_steps, step_size, Emax-eval_energy(at, do_PE=False))
	    at.info['ns_energy'] = final_E + eval_energy(at, do_PE=False, do_KE=True)
    else:
	if movement_args['MC_atom_velocities']:
	    exit_error("MC_atom_velocities only supported for FORTRAN calculator\n", 8)
	dz=0.0
	for i_MC_step in range(n_steps):
	    #DOC \item loop over atoms in random order
	    #DOC \begin{itemize}
	    at_list=list(range(len(at)))
	    rng.shuffle_in_place(at_list)
	    for i_at in at_list:
		#DOC \item propose single atom move
		if movement_args['MC_atom_uniform_rv']:
		    dx = rng.float_uniform(-step_size,step_size)
		    dy = rng.float_uniform(-step_size,step_size)
		    if not movement_args['2D']:
			dz = rng.float_uniform(-step_size,step_size)
		else:
		    dx = rng.normal(step_size)
		    dy = rng.normal(step_size)
		    if not movement_args['2D']:
			dz = rng.normal(step_size)
		orig_energy = at.info['ns_energy']
		orig_pos = at.get_positions()
		if ns_args['n_extra_data'] > 0:
		    orig_extra_data = at.arrays['ns_extra_data'].copy()
		new_pos = orig_pos.copy()
		new_pos[i_at,:] += (dx, dy,dz)
		at.set_positions(new_pos)
		#DOC \item accept/reject on E $<$ Emax
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
	    #DOC \end{itemize}
     #DOC \end{itemize}

    out = {}
    if n_accept_velo is not None:
	out['MC_atom_velo'] = (n_steps*len(at), n_accept_velo)
    out['MC_atom'] = (n_steps*len(at), n_accept)
    return out
#DOC \end{itemize}


def propose_volume_step(at, step_size):
    dV = rng.normal(step_size*len(at))
    orig_V = at.get_volume()
    new_V = orig_V+dV
    if new_V < 0: # negative number cannot be raised to fractional power, so this is only to avoid fatal error during the run
	new_V=abs(new_V)
	print "Warning, the step_size for volume change might be too big, resulted negativ new volume", step_size, dV, orig_V+dV
    #print "TRANSFORM", new_V, orig_V, dV, step_size, len(at)
    transform = np.identity(3)*(new_V/orig_V)**(1.0/3.0)
    p_accept = min(1.0, (new_V/orig_V)**len(at))
    return (p_accept, transform)

def propose_shear_step(at, step_size):
    # pick random vector
    rnd_vec_ind = rng.int_uniform(0, 3)
    # turn other two into orthonormal pair
    other_vec_ind = range(3)
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
    for i in range(3):
	vi = cell[i,:]
	vnorm_hat = np.cross(cell[(i+1)%3,:],cell[(i+2)%3,:])
	vnorm_hat /= np.sqrt(np.dot(vnorm_hat,vnorm_hat))
	min_aspect_ratio = min(min_aspect_ratio, abs(np.dot(vnorm_hat,vi)))
    return min_aspect_ratio/(vol**(1.0/3.0))

def do_cell_step(at, Emax, p_accept, transform):
    if p_accept < 1.0 and rng.float_uniform(0.0,1.0) > p_accept:
	return False

    # save old config and apply transformation
    orig_cell = at.get_cell()
    new_cell = np.dot(orig_cell,transform)
    new_vol = abs(np.dot(new_cell[0,:],np.cross(new_cell[1,:],new_cell[2,:])))
    # check size and shape constraints
    if new_vol > ns_args['max_volume_per_atom']*len(at):
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
    new_energy = eval_energy(at)
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
    items = possible_moves.items()

    for i in range(movement_args['cell_shape_equil_steps']):
	rng.shuffle_in_place(items)
	for key, propose_step_func in items:
	    (p_accept, transform) = propose_step_func(at, movement_args[key+"_step_size"])
	    do_cell_step(at, None, p_accept, transform)

def do_MC_swap_walk(at, movement_args, Emax, KEmax):
#DOC
#DOC \vspace*{\baselineskip}
#DOC do\_MC\_swap\_walk
#DOC \begin{itemize}
    if movement_args['swap_n_substeps'] <= 0:
	return {}

    Z = at.get_atomic_numbers()

    #DOC \item return if all atomic numbers are identical
    if (Z[:] == Z[0]).all(): 
        # don't try to swap when all atoms are the same
        return {}

    #DOC \item loop swap\_n\_substeps times:
    #DOC \begin{itemize}
    for i in range(movement_args['swap_n_substeps']):
	#DOC \item pick two atoms with distinct atomic numbers
	i1 = rng.int_uniform(0,len(at))
	i2 = rng.int_uniform(0,len(at))
	while Z[i1] == Z[i2]:
	    i2 = rng.int_uniform(0,len(at))
	p_1_orig = at.positions[i1,:].copy()
	p_2_orig = at.positions[i2,:].copy()
	at.positions[i1,:] = p_2_orig
	at.positions[i2,:] = p_1_orig
	if ns_args['n_extra_data'] > 0:
	    extra_data_1_orig = at.arrays['ns_extra_data'][i1,...].copy()
	    extra_data_2_orig = at.arrays['ns_extra_data'][i2,...].copy()
	    at.arrays['ns_extra_data'][i1,...] = extra_data_2_orig
	    at.arrays['ns_extra_data'][i2,...] = extra_data_1_orig

	#DOC \item accept swap if energy < Emax
	new_energy = eval_energy(at)

	if new_energy < Emax: # accept
	    at.info['ns_energy'] = new_energy
	else:
	    t = at.positions[i1,:]
	    at.positions[i1,:] = p_1_orig
	    at.positions[i2,:] = p_2_orig
	    if ns_args['n_extra_data'] > 0:
		at.arrays['ns_extra_data'][i1,...] = extra_data_1_orig
		at.arrays['ns_extra_data'][i2,...] = extra_data_2_orig

    #DOC \end{itemize}
    return {}
#DOC \end{itemize}

def do_MC_cell_volume_step(at, movement_args, Emax, KEmax):
    (p_accept, transform) = propose_volume_step(at, movement_args['MC_cell_volume_per_atom_step_size'])
    if do_cell_step(at, Emax, p_accept, transform):
	return {'MC_cell_volume_per_atom' : (1, 1) }
    else:
	return {'MC_cell_volume_per_atom' : (1, 0) }

def do_MC_cell_shear_step(at, movement_args, Emax, KEmax):
    (p_accept, transform) = propose_shear_step(at, movement_args['MC_cell_shear_step_size'])
    if do_cell_step(at, Emax, p_accept, transform):
	return {'MC_cell_shear' : (1, 1) }
    else:
	return {'MC_cell_shear' : (1, 0) }

def do_MC_cell_stretch_step(at, movement_args, Emax, KEmax):
    (p_accept, transform) = propose_stretch_step(at, movement_args['MC_cell_stretch_step_size'])
    if do_cell_step(at, Emax, p_accept, transform):
	return {'MC_cell_stretch' : (1, 1) }
    else:
	return {'MC_cell_stretch' : (1, 0) }


def do_MC_cell_walk(at, movement_args, Emax, KEmax):
#DOC
#DOC \vspace*{\baselineskip}
#DOC do\_MC\_cell\_walk
#DOC \begin{itemize}
    if movement_args['cell_n_substeps'] <= 0:
	return {}

    if movement_args['MC_cell_P'] <= 0.0 and rank == 0:
	sys.stderr.write("WARNING: doing MC_cell_walk but pressure == %f <= 0.0\n" % movement_args['MC_cell_P'])

    possible_moves = {
       'MC_cell_volume_per_atom': do_MC_cell_volume_step,
       'MC_cell_shear': do_MC_cell_shear_step,
       'MC_cell_stretch': do_MC_cell_stretch_step }

    items = possible_moves.items()
    out = {}
    #DOC \item loop cell\_n\_substeps times:
    #DOC \begin{itemize}
    for i in range(movement_args['cell_n_substeps']):
	rng.shuffle_in_place(items)
	#DOC \item do in random order
	#DOC \begin{itemize}
	for key, do_step_func in items:
	    #DOC \item volume move if rng $<$ attempt probablity
	    #DOC \item shear move if rng $<$ attempt probablity
	    #DOC \item stretch move if rng $<$ attempt probablity
	    step_rv = rng.float_uniform(0.0, 1.0)
	    if step_rv < movement_args[key+"_prob"]:
		accumulate_stats(out, do_step_func(at, movement_args, Emax, KEmax))
	#DOC \end{itemize}
    #DOC \end{itemize}

    return out
#DOC \end{itemize}


def do_atom_walk(at, movement_args, Emax, KEmax):
    if movement_args['atom_algorithm'] == 'MC':
	out = do_MC_atom_walk(at, movement_args, Emax, KEmax)
    elif movement_args['atom_algorithm'] == 'MD':
	out = do_MD_atom_walk(at, movement_args, Emax, KEmax)
    else:
	exit_error("do_atom_walk got unknown 'atom_algorithm' = '%s'\n" % movement_args['atom_algorithm'], 5)

    return out

def rand_perturb_energy(energy, perturbation, Emax=None):
    if Emax is None:
	if abs(energy)<= 1.0:
	    energy += rng.float_uniform(-1.0,0.0)*perturbation
	else:
	    energy *= (1.0+rng.float_uniform(-1.0,0.0)*perturbation)
    else:
	if abs(energy) <= 1.0:
	    pert = rng.float_uniform(-1.0,0.0)*perturbation
	    n_tries = 0
	    while energy+pert >= Emax and n_tries < 100:
		pert = rng.float_uniform(-1.0,0.0)*perturbation
		n_tries += 1
	    if energy+pert >= Emax:
		print print_prefix, "WARNING: failed to do random energy perturbation below Emax ", energy, Emax
	    energy += pert
	else:
	    pert = 1.0 + rng.float_uniform(-1.0,0.0)*perturbation
	    n_tries = 0
	    while energy*pert >= Emax and n_tries < 100:
		pert = 1.0 + rng.float_uniform(-1.0,0.0)*perturbation
		n_tries += 1
	    if energy*pert >= Emax:
		print print_prefix, "WARNING: failed to do random energy perturbation below Emax ", energy, Emax
	    energy *= pert

    return energy

def walk_single_walker(at, movement_args, Emax, KEmax):
    """Do random walk on a single atoms object."""
#DOC
#DOC \vspace*{\baselineskip}
#DOC walk\_single\_walker
#DOC \begin{itemize}

    #DEBUG print "walk_single_walker start w/o PE ", eval_energy(at, do_PE=False), " w/ PE ", eval_energy(at), "Emax ", Emax #DEBUG
    possible_moves = [do_atom_walk, do_MC_cell_walk, do_MC_swap_walk]
    if movement_args['do_velocities'] and movement_args['atom_velo_random_order_perturb']:
	possible_moves.append(do_MC_atom_velo_walk)
    out = {}

    #DOC \item loop n\_steps times
    #DOC \begin{itemize}
    for i_step in range(movement_args['n_steps']):
	#DOC \item do in random order
	#DOC \begin{itemize}
	rng.shuffle_in_place(possible_moves)
	for move in possible_moves:
	    #DOC \item do\_(MC $|$ MD)\_atom\_walk
	    #DOC \item do\_MC\_cell\_walk
	    #DOC \item do\_MC\_swap\_walk
	    accumulate_stats(out, move(at, movement_args, Emax, KEmax))
	#DOC \end{itemize}
    #DOC \end{itemize}

    #DOC \item create list of possible unblocked moves, each type x repeated n\_steps\_unblocked\_x times (x=do\_atom\_walk, do\_MC\_cell\_*\_step, do\_MC\_swap\_walk)
    possible_moves = ( [ do_atom_walk ] * movement_args['n_steps_unblocked_atom'] +
		       [ do_MC_cell_volume_step ] * movement_args['n_steps_unblocked_cell_volume'] +
		       [ do_MC_cell_shear_step ] * movement_args['n_steps_unblocked_cell_shear'] +
		       [ do_MC_cell_stretch_step ] * movement_args['n_steps_unblocked_cell_stretch'] +
		       [ do_MC_swap_walk ] * movement_args['n_steps_unblocked_swap'] )
    #DOC \begin{itemize}
    #DOC \item do in random order
    rng.shuffle_in_place(possible_moves)
    for move in possible_moves:
	accumulate_stats(out, move(at, movement_args, Emax, KEmax))
    #DOC \end{itemize}

    #DOC \item perturb final energy by random\_energy\_perturbation
    # perturb final energy
    at.info['ns_energy'] = rand_perturb_energy(at.info['ns_energy'],ns_args['random_energy_perturbation'],Emax)

    #DEBUG print "walk_single_walker end ", eval_energy(at, do_PE=False), eval_energy(at) #DEBUG

    return out
#DOC \end{itemize}



def max_energy(walkers, n):
    """Collect the current energies of the walkers from all the processes and chooses the right number of highest energies to be culled"""
    # do local max
    energies_loc = np.array([ at.info['ns_energy'] for at in walkers])
    if comm is not None:
	energies = np.zeros( (comm.size*len(energies_loc)) )
	# comm.barrier() #BARRIER
	comm.Allgather( [ energies_loc, MPI.DOUBLE ], [ energies, MPI.DOUBLE ] )
	energies = energies.flatten()
    else:
	energies = energies_loc
   
    # n is n_cull
    Emax_ind = energies.argsort()[-1:-n-1:-1]
    Emax = energies[Emax_ind]
    # WARNING: assumes that each node has equal number of walkers
    rank_of_max = np.floor(Emax_ind/len(walkers)).astype(int)
    ind_of_max = np.mod(Emax_ind,len(walkers))

    return (Emax, rank_of_max, ind_of_max)

def median_PV(walkers):
    # do local max
    PVs_loc = np.array([ eval_energy(at, do_PE=False, do_KE=False, do_PV=True) for at in walkers])

    if comm is not None:
	PVs = np.zeros( (comm.size*len(PVs_loc)) )
	# comm.barrier() #BARRIER
	comm.Allgather( [ PVs_loc, MPI.DOUBLE ], [ PVs, MPI.DOUBLE ] )
	PVs = PVs.flatten()
    else:
	PVs = PVs_loc

    if len(PVs) % 2 == 0:
	PV_median = (PVs[int(len(PVs)/2)-1] + PVs[int(len(PVs)/2)])/2.0
    else:
	PV_median = PVs[int(len(PVs)/2)]

    return PV_median

def adjust_step_sizes(walk_stats, movement_args, comm):
    """
    Adjust step size to keep the acceptance ratio at the desired level.
    """
    for key in walk_stats:
	(n_try, n_accept) = walk_stats[key]
	if comm is not None:
	    n_try_s = np.array( [n_try], dtype = np.int)
	    n_accept_s = np.array( [n_accept], dtype = np.int)
	    n_try_g = np.zeros( (1), dtype=np.int)
	    n_accept_g = np.zeros( (1), dtype=np.int)
	    # comm.barrier() #BARRIER
	    comm.Allreduce([n_try_s, MPI.INT], [n_try_g, MPI.INT], MPI.SUM)
	    comm.Allreduce([n_accept_s, MPI.INT], [n_accept_g, MPI.INT], MPI.SUM)
	    n_try = n_try_g[0]
	    n_accept = n_accept_g[0]

	if n_try > 0:
	    rate = float(n_accept)/float(n_try)

	    if comm is None or comm.rank == 0:
		print print_prefix, "accept rate for %s = %f (%d)" % (key, rate, n_try)

	    if key.find("MC") == 0:
		min_rate = movement_args['MC_adjust_min_rate']
		max_rate = movement_args['MC_adjust_max_rate']
		suffix="step_size"
	    elif key.find("MD") == 0:
		min_rate = movement_args['MD_adjust_min_rate']
		max_rate = movement_args['MD_adjust_max_rate']
		suffix="timestep"
	    else:
		exit_error("adjust_step_size got key '%s', neither MC nor MD\n" % key, 5)

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
		print print_prefix, "adjust_step_sizes adjusted %s %s to %f" % (key, dir, movement_args[key+"_"+suffix])

def zero_stats(d, movement_args):
    for key in movement_args:
	m = re.search('(.*)_(step_size|timestep)$', key)
	if m is not None:
	    d[m.group(1)] = (0, 0)

def accumulate_stats(d_cumul, d):
    for key in d:
	if key in d_cumul:
	    d_cumul[key] = tuple([i1+i2 for i1,i2 in zip(d[key],d_cumul[key])])
	else:
	    d_cumul[key] = d[key]


# figure out n_steps to walk on each iteration to get correct expected number
def set_n_steps(prop):
    if movement_args[prop+'_expected'] > 0:
	if movement_args[prop] > 0:
	    exit_error("Got both "+prop+" and "+prop+"_expected, conflict\n", 5)

	if rank == 0:
	    print "Calculating %s from %s_expected=%d" % (prop, prop, movement_args[prop+'_expected'])

	if max_n_cull_per_task*size == n_cull and n_extra_walk_per_task == 0: # no extra walkers
	    movement_args[prop] = movement_args[prop+'_expected']
	    if rank == 0:
		print "No extra walkers (n_cull mod n_tasks == 0), trivial, so average n_walks at kill is 1, and %s=%s_expected" % (prop, prop)
	else:
	    # f_c = n_c/n_t [ fraction of total that are culled (and walked once)]
	    # f_e = n_e/(n_t-n_c) [fraction of ones that aren't culled that are also walked ]

	    # n(1) <- f_c n_t + (1-f_c)*(1-f_e) n(1)
	    # n(l >= 2) <- (1-f_c)(1-f_e) n(l) + (1-f_c) f_e n(l-1)

	    # f = (1-f_c)f_e / (f_c+f_e-f_c f_e)

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
		print "Calculated average n_walks at kill = ",n_walks, " from f_cull ",f_cull," f_extra (due to otherwise idle processors doing walks and explicitly requested extra walks) ",f_extra
                if n_walks < 1:
                    print "WARNING: got average n_walks < 1, but will always do at least 1 walk, so effective %s_expected will be higher than requested" % prop
                print "Setting %s = ceiling(%s_expected/n_walks)" % (prop, prop)
	    movement_args[prop] = int(math.ceil(movement_args[prop+'_expected']/n_walks))

    else:
	if movement_args[prop] > 0 and rank == 0:
	    print "WARNING: using absolute number of "+prop

    if rank == 0:
        print "Final value of %s=%d" % (prop, movement_args[prop])

def additive_init_config(at, Emax):
    if do_calc_lammps:
	exit_error("got energy %f ceiling %f, but python additive_init_config doesn't work with LAMMPS, since it varies list of atoms\n" % (energy, ns_args['start_energy_ceiling']), 10)
    pos = at.get_positions()
    for i_at in range(1,len(at)):
	at_new = at[0:i_at+1]
	if do_calc_quip or do_calc_lammps:
	    at_new.set_calculator(at.get_calculator())
	at_new.set_positions(pos[0:i_at+1,:])
	success = False
	for i_try in range(10):
	    pos[i_at,:] = np.dot(at_new.get_cell(), rng.float_uniform(0.0, 1.0, (3) ))
	    at_new.set_positions(pos[0:i_at+1,:])
	    energy = eval_energy(at_new)
	    if energy < Emax:
		success = True
		break
	if not success:
	    exit_error("Failed 10 times to insert atom %d with Emax %f" % (i_at, Emax), 7)
    at.set_positions(pos)
    return energy

def save_snapshot(id):
    """
    Save the current walker configurations' as a snapshot in the file ``out_file_prefix.iter.rank.config_file_format``
    """
    #QUIP_IO if have_quippy:
	#QUIP_IO snapshot_io = quippy.CInOutput(ns_args['out_file_prefix']+'snapshot.%s.%d.extxyz' % (id,rank), action=quippy.OUTPUT)
    #QUIP_IO else:
	#QUIP_IO try:
	    #QUIP_IO snapshot_file=ns_args['out_file_prefix']+'snapshot.'+('%d' % id)+'.%05d.'+('%04d' % rank)+'.extxyz'
	#QUIP_IO except:
	    #QUIP_IO snapshot_file=ns_args['out_file_prefix']+'snapshot.'+id+'.%05d.'+('%04d' % rank)+'.extxyz'
    try:
	snapshot_io = open(ns_args['out_file_prefix']+'snapshot.%s.%d.%s' % (id,rank, ns_args['config_file_format']), "w")
    except:
	snapshot_io = open(ns_args['out_file_prefix']+'snapshot.%d.%d.%s' % (id,rank, ns_args['config_file_format']), "w")

    for at in walkers:
	#QUIP_IO if have_quippy:
	    #QUIP_IO at.write(snapshot_io)
	#QUIP_IO else:
	    #QUIP_IO ase.io.write(snapshot_file % i_at, ase.Atoms(at))
        at.info['iter']=id
	ase.io.write(snapshot_io, at, format=ns_args['config_file_format'])

    snapshot_io.close()

def clean_prev_snapshot(prev_snapshot_iter):
    if prev_snapshot_iter is not None:
	snapshot_file=ns_args['out_file_prefix']+'snapshot.%d.%d.%s' % (prev_snapshot_iter, rank, ns_args['config_file_format'])
	try:
	    os.remove(snapshot_file)
	except:
	    print print_prefix, ": WARNING: Failed to delete '%s'" % snapshot_file

def do_ns_loop():
    """ 
    This is the main nested sampling loop, doing the iterations.
    """
    global print_prefix
    global cur_config_ind

    if rank == 0:
	nD = 3
	if movement_args['2D']:
	    nD = 2
	if ns_args['restart_file'] == '':
	    if movement_args['do_velocities']:
		nExtraDOF = 0
	    else:
		nExtraDOF = n_atoms*nD
	    energy_io.write("%d %d %d\n" % (ns_args['n_walkers'], ns_args['n_cull'], nExtraDOF) )

    ## print print_prefix, ": random state ", np.random.get_state()
    ## if rank == 0:
	## print print_prefix, ": common random state ", common_random_state

    if ns_args['debug'] >= 10 and size <= 1:
	for at in walkers:
	    at.info['n_walks'] = 0

    for at in walkers:
        at.info['KEmax']=KEmax 
	if movement_args['MC_cell_P'] > 0:
	    print rank, ": initial enthalpy ", at.info['ns_energy'], " PE ", eval_energy(at, do_KE=False, do_PV=False), " KE ", eval_energy(at, do_PE=False, do_PV=False)
	else:
	    print rank, ": initial energy ", at.info['ns_energy']

    walk_stats_cumul={}
    zero_stats(walk_stats_cumul, movement_args)

    initial_time = time.time()
    prev_time = initial_time

    Emax_of_step = None
    Emax_save = []
    i_ns_step_save = []
    walker_list = []

    verbose=False

    # to avoid errors of unassigned values, if in case of a restart the final number of iter is the same as the satring, stop.
    if start_first_iter == ns_args['n_iter']:
	print "WARNING: Increase the n_iter_per_walker variable in the input if you want NS cycles to be performed."
        exit_error("satring iteration and the total number of required iterations are the same,hence no NS cycles will be performed\n",11)

    # actual iteration cycle starts here
    for i_ns_step in range(start_first_iter, ns_args['n_iter']):
	print_prefix="%d %d" % (rank, i_ns_step)

        if ns_args['debug'] >= 4 and ns_args['track_configs']:
            for at in walkers:
                print print_prefix, "INFO: 10 config_ind ", at.info['config_ind'], " from ", at.info['from_config_ind'], " at ", at.info['config_ind_time']

	if movement_args['adjust_step_interval'] < 0:
	    zero_stats(walk_stats_cumul, movement_args)

	if ns_args['debug'] >= 20:
	    print print_prefix, "%30s" % ": LOOP_TE START 00 ",i_ns_step, [ eval_energy(at) for at in walkers ]
	    print print_prefix, "%30s" % ": LOOP_PE START 01 ",i_ns_step, [ eval_energy(at, do_KE=False) for at in walkers ]
	    print print_prefix, "%30s" % ": LOOP_X START 02 ",i_ns_step, [ at.positions[0,0] for at in walkers ]

	# get list of highest energy configs
	(Emax, cull_rank, cull_ind) = max_energy(walkers, n_cull)
	Emax_next = Emax[-1]
	if rank == 0 and Emax_of_step is not None and Emax[0] > Emax_of_step:
	    print print_prefix, ": WARNING: energy above Emax ", Emax_of_step, " bad energies: ", Emax[np.where(Emax > Emax_of_step)], cull_rank[np.where(Emax > Emax_of_step)], cull_ind[np.where(Emax > Emax_of_step)]
	    # comm.barrier()
	    # exit_error("Energy above Emax\n", 5)

	if rank == 0 and (i_ns_step > start_first_iter and Emax_next >= Emax_of_step):
	    print "WARNING: Emax not decreasing ",Emax_of_step, Emax_next
	Emax_of_step=Emax_next

	if ns_args['min_Emax'] is not None and Emax_of_step < ns_args['min_Emax']:
	    if rank == 0:
                # if the termination was set by a minimum energy, and it is reached, stop.
		print "Leaving loop because Emax=",Emax_of_step," < min_Emax =",ns_args['min_Emax']
	    break

	if rank == 0:
	    cur_time=time.time()
	    if (cur_time > prev_time+60 or i_ns_step == 0 or i_ns_step == ns_args['n_iter'] or (ns_args['n_iter'] > 0 and i_ns_step % max(int(ns_args['n_iter']/1000),1) == 0)):
		print i_ns_step, "Emax_of_step ", Emax_of_step, " loop time ", cur_time-prev_time
		prev_time = cur_time

	cull_list=[None] * size
	for r in range(size):
	    entries_for_this_rank = np.where(cull_rank == r)[0]
	    cull_list[r] = cull_ind[entries_for_this_rank]
            if rank == 0 and ns_args['debug'] >= 4 and len(cull_ind[entries_for_this_rank]) > 0:
                print print_prefix, "INFO: 20 cull ", cull_ind[entries_for_this_rank], " on ",r


	# record Emax walkers energies
	if rank == 0:
            for E in Emax:
		energy_io.write("%d %.60f\n" % (i_ns_step, E))
	    energy_io.flush()

            ## Save the energies and corresponding iteration numbers in a list then print them out only when printing a snapshot
            #Emax_save.extend(Emax)
            #i_ns_step_save.extend(n_cull*[i_ns_step])
            ## if it is time to print (i.e. at the same iteration when a snapshot is written, or at every iter if no snapshots - for smooth restarts)
	    #if ns_args['snapshot_interval'] < 0 or i_ns_step % ns_args['snapshot_interval'] == ns_args['snapshot_interval']-1:
	    #    for istep,E in zip(i_ns_step_save,Emax_save):
	    #        energy_io.write("%d %.60f\n" % (istep, E))
	    #    energy_io.flush()
            #    #empty the save lists, so they are ready for the next bunch of saved energies
            #    Emax_save[:]=[]
            #    i_ns_step_save[:]=[]

	# record Emax walkers configurations
	if cull_list[rank] is not None:
	    for i in cull_list[rank]:
		if ns_args['debug'] >= 10 and size <= 1:
		    print print_prefix, "walker killed at age ",walkers[i].info['n_walks']
		walkers[i].info['volume'] = walkers[i].get_volume()
		walkers[i].info['ns_P'] = movement_args['MC_cell_P']
		walkers[i].info['iter'] = i_ns_step
		if walkers[i].has('masses') and walkers[i].has('momenta'):
		    walkers[i].info['ns_KE'] = walkers[i].get_kinetic_energy()
		#QUIP_IO if have_quippy:
		    #QUIP_IO walkers[i].write(traj_io)
		#QUIP_IO else:
		    #QUIP_IO ase.io.write(traj_file % i_ns_step, ase.Atoms(walkers[i]))

                # store culled config in list to be written (when snapshot_interval has passed) every traj_interval steps
	        if ns_args['traj_interval'] > 0 and i_ns_step % ns_args['traj_interval'] == 0:
		    walker_list.append(walkers[i].copy())

                # if tracking all configs, save this one that has been culled
                if track_traj_io is not None:
                    at = walkers[i].copy()
                    at.info['culled'] = True
                    ase.io.write(track_traj_io, at, format=ns_args['config_file_format'])

	# print the recorded Emax walkers configurations to output file
        if ns_args['snapshot_interval'] < 0 or i_ns_step % ns_args['snapshot_interval'] == ns_args['snapshot_interval']-1:
	    for at in walker_list:
	        ase.io.write(traj_io, at, format=ns_args['config_file_format'])
	    traj_io.flush()
	    walker_list=[]

	# calculate how many will be culled on each rank
	n_cull_of_rank = np.array([ sum(cull_rank == r) for r in range(size) ])

	# label configs to be culled
	status = np.empty( (size, n_walkers), np.object_)
	status[:,:] = ''
	for r in range(size):
	    status[r,cull_ind[np.where(cull_rank == r)[0]]] = 'c_t'

	if ns_args['debug'] >= 10:
	    initial_PE_loc = [ eval_energy(at, do_KE=False) for at in walkers ]
	    initial_PE = np.array(comm.allgather(initial_PE_loc)).flatten()
	    initial_changed = initial_PE[np.where(status.flatten() == 'c_t')]
	    initial_unchanged = initial_PE[np.where(status.flatten() == '')]

	if ns_args['debug'] >= 30:
	    for r in range(len(status)):
		print print_prefix, ": initial status ", r, [ s for s in status[r,:] ]

	# find load balance by cloning on top of excess maxima
	recv_ind=[]
	recv_rank=[]
	send_ind=[]
	send_rank=[]
	cull_inds_to_remove=[]

	if n_cull > 1: # send/recv for fixing load balance
	    # CHECK FOR RANDOMNESS ISSUES AND WHICH NODES ARE USED FOR CLONES
	    for r in range(size):
		# maybe remote_r should be chosen completely randomly, rather than close to task of extra culled configs
		for dr in np.array(zip(np.array(range(1,size)), -np.array(range(1,size)))).flatten():
		    if n_cull_of_rank[r] <= max_n_cull_per_task: # not too many that need to be culled on this rank
			break
		    # this rank has too many to cull, must receive replacement from another node
		    remote_r = (r+dr) % size
		    if n_cull_of_rank[remote_r] < max_n_cull_per_task: # send from r+dr to r
			n_transfer = min(n_cull_of_rank[r]-max_n_cull_per_task, max_n_cull_per_task-n_cull_of_rank[remote_r])
			recv_rank.extend([r]*n_transfer)
			send_rank.extend([remote_r]*n_transfer)
			local_ind = np.where(status[r,:] == 'c_t')[0][0:n_transfer]
			recv_ind.extend(local_ind)
			remote_ind = np.where(status[remote_r,:] == '')[0][0:n_transfer]
			send_ind.extend(remote_ind)
			status[r,local_ind] = 'c_s'
			status[remote_r,remote_ind] = 'c_t_a'
			n_cull_of_rank[r] -= n_transfer
			n_cull_of_rank[remote_r] += n_transfer


	# save local random state, and switch to common one
	rng.switch_to_common()

	# select clones
	for r in range(size):
	    list_clone_target = np.where(status[r,:] == 'c_t')[0]
	    # assign clones
	    n_remaining_clones = len(list_clone_target)
	    while n_remaining_clones > 0:
		remote_r = rng.int_uniform(0,size)
		n_avail_remote = sum(status[remote_r,:] == '')
		if n_avail_remote > 0: # something is available on remote_r
		    # send from random avail walker on remote_r to clone_target on r
		    n_transfer = min(n_remaining_clones, n_avail_remote)

		    # set ranks
		    send_rank.extend([remote_r]*n_transfer)
		    recv_rank.extend([r]*n_transfer)

		    # set indices
		    r_is = []
		    for ii in range(n_transfer):
			r_i = rng.int_uniform(0, n_walkers)
			while status[remote_r,r_i] != '':
			    r_i = rng.int_uniform(0, n_walkers)
			# now r_i should be something with status ''
			status[remote_r,r_i] = 'c_s'
			r_is.append(r_i)
		    send_ind.extend(r_is)

		    status[r,list_clone_target[0:n_transfer]] = 'c_t_a'
		    recv_ind.extend(list_clone_target[0:n_transfer])

		    if n_transfer < len(list_clone_target):
			list_clone_target = list_clone_target[n_transfer:]
		    n_remaining_clones -= n_transfer

	if ns_args['debug'] >= 20:
	    print print_prefix, "%30s" % ": LOOP_TE POST_LOC_CLONE 15 ",i_ns_step, [ eval_energy(at) for at in walkers ]
	    print print_prefix, "%30s" % ": LOOP_PE POST_LOC_CLONE 16 ",i_ns_step, [ eval_energy(at, do_KE=False) for at in walkers ]
	    print print_prefix, "%30s" % ": LOOP_X POST_LOC_CLONE 17 ",i_ns_step, [ at.positions[0,0] for at in walkers ]

	# make into numpy arrays so that mathematical operations will work
	send_rank = np.array(send_rank)
	send_ind = np.array(send_ind)
	recv_rank = np.array(recv_rank)
	recv_ind = np.array(recv_ind)

	if ns_args['debug'] >= 10:
	    if rank == 0:
		for i in range(len(send_rank)):
		    print print_prefix, "send from ",send_rank[i],send_ind[i]," to ",recv_rank[i], recv_ind[i]

	# save new common state, and restore to local state
	rng.switch_to_local()

	if n_cull == 1:
	    if send_rank[0] == recv_rank[0] and send_rank[0] == rank: # local copy
		walkers[recv_ind[0]].set_positions(walkers[send_ind[0]].get_positions())
		walkers[recv_ind[0]].set_cell(walkers[send_ind[0]].get_cell())
		if movement_args['do_velocities']:
		    walkers[recv_ind[0]].set_velocities(walkers[send_ind[0]].get_velocities())
		if ns_args['n_extra_data'] > 0:
		    walkers[recv_ind[0]].arrays['ns_extra_data'][...] = walkers[send_ind[0]].arrays['ns_extra_data']
		if ns_args['track_configs']:
                    walkers[recv_ind[0]].info['config_ind'] = walkers[send_ind[0]].info['config_ind']
                    walkers[recv_ind[0]].info['from_config_ind'] = walkers[send_ind[0]].info['from_config_ind']
                    walkers[recv_ind[0]].info['config_ind_time'] = walkers[send_ind[0]].info['config_ind_time']
		walkers[recv_ind[0]].info['ns_energy'] = eval_energy(walkers[recv_ind[0]])
		if ns_args['debug'] >= 10 and size <= 1:
		    walkers[recv_ind[0]].info['n_walks'] = 0
	    else: # need send/recv
		n_send = 3*(n_atoms + 3)
		if movement_args['do_velocities']:
		    n_send += 3*n_atoms
		if ns_args['n_extra_data'] > 0:
		    n_send += ns_args['n_extra_data']*n_atoms
                if ns_args['track_configs']:
                    n_send += 3
		buf = np.zeros ( n_send )
		if send_rank[0] == rank: # only one config is sent/received
		    buf_o = 0
		    buf[buf_o:buf_o+3*n_atoms] = walkers[send_ind[0]].get_positions().reshape( (3*n_atoms) ); buf_o += 3*n_atoms
		    buf[buf_o:buf_o+3*3] = walkers[send_ind[0]].get_cell().reshape( (3*3) ); buf_o += 3*3
		    if movement_args['do_velocities']:
			buf[buf_o:buf_o+3*n_atoms] = walkers[send_ind[0]].get_velocities().reshape( (3*n_atoms) ); buf_o += 3*n_atoms
		    if ns_args['n_extra_data'] > 0:
			buf[buf_o:buf_o+ns_args['n_extra_data']*n_atoms] = walkers[send_ind[0]].arrays['ns_extra_data'].reshape( (ns_args['n_extra_data']*n_atoms) ); buf_o += ns_args['n_extra_data']*n_atoms
                    if ns_args['track_configs']:
                        buf[buf_o] = walkers[send_ind[0]].info['config_ind']; buf_o += 1
                        buf[buf_o] = walkers[send_ind[0]].info['from_config_ind']; buf_o += 1
                        buf[buf_o] = walkers[send_ind[0]].info['config_ind_time']; buf_o += 1
		    comm.Send([buf,  MPI.DOUBLE], dest=recv_rank[0], tag=100)
		elif recv_rank[0] == rank:
		    comm.Recv([buf, MPI.DOUBLE], source=send_rank[0], tag=100)
		    buf_o = 0
		    walkers[recv_ind[0]].set_positions(buf[buf_o:buf_o+3*n_atoms].reshape( (n_atoms, 3) )); buf_o += 3*n_atoms
		    walkers[recv_ind[0]].set_cell(buf[buf_o:buf_o+3*3].reshape( (3, 3) )); buf_o += 3*3
		    if movement_args['do_velocities']:
			walkers[recv_ind[0]].set_velocities(buf[buf_o:buf_o+3*n_atoms].reshape( (n_atoms, 3) )); buf_o += 3*n_atoms
		    if ns_args['n_extra_data'] > 0:
			walkers[recv_ind[0]].arrays['ns_extra_data'][...] = buf[buf_o:buf_o+3*n_atoms].reshape( walkers[recv_ind[0]].arrays['ns_extra_data'].shape ); buf_o += ns_args['n_extra_data']*n_atoms
		    if ns_args['track_configs']:
                        walkers[recv_ind[0]].info['config_ind'] = int(buf[buf_o]); buf_o += 1
                        walkers[recv_ind[0]].info['from_config_ind'] = int(buf[buf_o]); buf_o += 1
                        walkers[recv_ind[0]].info['config_ind_time'] = int(buf[buf_o]); buf_o += 1
		    walkers[recv_ind[0]].info['ns_energy'] = eval_energy(walkers[recv_ind[0]])

	else: # complicated construction of sending/receiving buffers
	    # figure out how much is sent per config
	    n_data_per_config = 1+3*(n_atoms + 3)
	    if movement_args['do_velocities']:
		n_data_per_config += 3*n_atoms
	    if ns_args['n_extra_data'] > 0:
		n_data_per_config += ns_args['n_extra_data']*n_atoms
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
		if ns_args['n_extra_data'] > 0:
		    send_data[data_o:data_o+ns_args['n_extra_data']*n_atoms] = walkers[i_send].arrays['ns_extra_data'].reshape( (ns_args['n_extra_data']*n_atoms) ); data_o += ns_args['n_extra_data']*n_atoms
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
		if ns_args['n_extra_data'] > 0:
		    walkers[i_recv].arrays['ns_extra_data'][...] = recv_data[data_o:data_o+ns_args['n_extra_data']*n_atoms].reshape( walkers[i_recv].arrays['ns_extra_data'].shape ); data_o += ns_args['n_extra_data']*n_atoms
		if ns_args['track_configs']:
                    walkers[i_recv].info['config_ind'] = int(recv_data[data_o]); data_o += 1
                    walkers[i_recv].info['from_config_ind'] = int(recv_data[data_o]); data_o += 1
                    walkers[i_recv].info['config_ind_time'] = int(recv_data[data_o]); data_o += 1
		recv_displ_t[r_send] = data_o

	if ns_args['debug'] >= 20:
	    print print_prefix, "%30s" % ": LOOP_TE POST_CLONE 20 ",i_ns_step, [ eval_energy(at) for at in walkers ]
	    print print_prefix, "%30s" % ": LOOP_PE POST_CLONE 21 ",i_ns_step, [ eval_energy(at, do_KE=False) for at in walkers ]
	    print print_prefix, "%30s" % ": LOOP_X POST_CLONE 22 ",i_ns_step, [ at.positions[0,0] for at in walkers ]

        if ns_args['track_configs']:
            # loop over _all_ clone targets and increment cur_config_ind, setting appropriate configs' new config_ind as needed 
            for r in range(size):
                clone_walk_ind = np.where(status[r,:] == 'c_t_a')[0]
                for i_at in clone_walk_ind:
                    if r == rank:
                        walkers[i_at].info['from_config_ind'] = walkers[i_at].info['config_ind']
                        walkers[i_at].info['config_ind'] = cur_config_ind
                        walkers[i_at].info['config_ind_time'] = i_ns_step
                    cur_config_ind += 1
	# move cloned walkers

	# walk clone targets
	if ns_args['debug'] >= 4:
	    for i in np.where(status[rank,:] == 'c_s')[0]:
		print print_prefix, "INFO: 30 clone source ", rank, i
	clone_walk_ind = np.where(status[rank,:] == 'c_t_a')[0]
	for i_at in clone_walk_ind:
	    if ns_args['debug'] >= 4:
		print print_prefix, "INFO: 40 WALK clone_target ", rank, i_at
	    walk_stats = walk_single_walker(walkers[i_at], movement_args, Emax_of_step, KEmax)
            # if tracking all configs, save this one that has been walked
            if track_traj_io is not None:
                walkers[i_at].info['iter'] = i_ns_step
                ase.io.write(track_traj_io, walkers[i_at], format=ns_args['config_file_format'])
	    #print "WALK on rank ", rank, "at iteration ", i_ns_step, " walker ", i_at
	    if ns_args['debug'] >= 10 and size <= 1:
		walkers[i_at].info['n_walks'] += movement_args['n_steps']
	    accumulate_stats(walk_stats_cumul, walk_stats)

	if ns_args['debug'] >= 20:
	    print print_prefix, "%30s" % ": LOOP_TE POST_CLONE_WALK 25 ",i_ns_step, [ eval_energy(at) for at in walkers ]
	    print print_prefix, "%30s" % ": LOOP_PE POST_CLONE_WALK 26 ",i_ns_step, [ eval_energy(at, do_KE=False) for at in walkers ]

	# check that everything that should have been changed has, and things that shouldn't have, haven't
	if ns_args['debug'] >= 10:
            final_PE_loc = [ eval_energy(at, do_KE=False) for at in walkers ]
            if comm is not None:
                final_PE = np.array(comm.allgather(final_PE_loc)).flatten()
            else:
                final_PE = final_PE_loc
	    if rank == 0:
		final_status = status.flatten()
		for e in initial_unchanged:
		    if e not in final_PE:
			print "initial_PE ", initial_PE
			print "final_PE ", final_PE
			print "final_status ", final_status
			print "WARNING: energy that should have been unchanged ", e," missing from final energies"
		for e in initial_changed:
		    if e in final_PE:
			print "initial_PE ", initial_PE
			print "final_PE ", final_PE
			print "final_status ", final_status
			print "WARNING: energy that should have been changed ", e," still there in final energies"


	# walk extras
	if not ns_args['no_extra_walks_at_all']:
	    for ii in range(max_n_cull_per_task - len(clone_walk_ind)+n_extra_walk_per_task):
		r_i = rng.int_uniform(0, n_walkers)
		# WARNING: this may select walkers for extra walks multiple times, yet never re-walk ones that were walked as clone targets
		while status[rank,r_i] != '' and status[rank,r_i] != 'c_s':
		    r_i = rng.int_uniform(0, n_walkers)
		if ns_args['debug'] >= 4:
		    print print_prefix, "INFO: 50 WALK extra ",rank, r_i
		walk_stats = walk_single_walker(walkers[r_i], movement_args, Emax_of_step, KEmax)
                # if tracking all configs, save this one that has been walked
                if track_traj_io is not None:
                    walkers[i_at].info['iter'] = i_ns_step
                    ase.io.write(track_traj_io, walkers[i_at], format=ns_args['config_file_format'])
	        #print "WALK EXTRA on rank ", rank, "at iteration ", i_ns_step, " walker ", r_i
		if ns_args['debug'] >= 10 and size <= 1:
		    walkers[r_i].info['n_walks'] += movement_args['n_steps']
		accumulate_stats(walk_stats_cumul, walk_stats)

	if movement_args['adjust_step_interval'] != 0 and i_ns_step % abs(movement_args['adjust_step_interval']) == abs(movement_args['adjust_step_interval'])-1:
	    adjust_step_sizes(walk_stats_cumul, movement_args, comm)
	    zero_stats(walk_stats_cumul, movement_args)

	if ns_args['debug'] >= 20:
	    print print_prefix, "%30s" % ": LOOP_TE END 30 ",i_ns_step, [ eval_energy(at) for at in walkers ]
	    print print_prefix, "%30s" % ": LOOP_PE END 31 ",i_ns_step, [ eval_energy(at,do_KE=False) for at in walkers ]

	if ns_args['debug'] >= 30:
	    for r in range(len(status)):
		print print_prefix, ": final status ", r, [ s for s in status[r,:] ]

	if ns_args['snapshot_interval'] > 0 and i_ns_step % ns_args['snapshot_interval'] == ns_args['snapshot_interval']-1:
	    global prev_snapshot_iter
	    save_snapshot(i_ns_step)
	    clean_prev_snapshot(prev_snapshot_iter)
	    prev_snapshot_iter = i_ns_step
    cur_time = time.time()
    if rank == 0:
	print "LOOP TIME total ",cur_time-initial_time, " per iter ", (cur_time-initial_time)/(i_ns_step+1)

def main():
	""" Main function """
        global movement_args
        global ns_args, start_first_iter
        global max_n_cull_per_task
        global size, rank, comm, rng, np, sys
        global n_cull, n_walkers, n_walkers_per_task
        global n_extra_walk_per_task, prev_snapshot_iter
        global do_calc_quip, do_calc_lammps, do_calc_internal, do_calc_fortran
        global energy_io, traj_io, walkers
        global n_atoms, KEmax, pot
        global MPI, quippy, f_MC_MD
        global track_traj_io, cur_config_ind

	import sys

	stacktrace.listen()

	print_prefix=""
	if len(sys.argv) != 1 and len(sys.argv) != 2:
	    usage()
	    sys.exit(1)

	use_mpi=True
	if len(sys.argv) == 2:
	    if sys.argv[1] == "-no_mpi":
		use_mpi=False
	    else:
		usage()
		sys.exit(1)

	try:
	    import quippy
	    have_quippy=True
	except:
	    have_quippy=False

	# initialize mpi
	comm = None
	rank = 0
	size = 1
	if use_mpi:
	    try:
		from mpi4py import MPI
	    except:
		sys.stderr.write("Failed to import mpi4py\n")
		sys.exit(10)
	    comm = MPI.COMM_WORLD

	if use_mpi:
	    try:
		rank = comm.Get_rank()
		size = comm.Get_size()
	    except:
		exit_error("Failed to get rank or size\n", 10)

	if comm is not None:
	    print "comm ", comm, " size ", size, " rank ", rank

	# read inputs on root, then bcast
	if rank == 0:
	    lines=sys.stdin.readlines()
	    if len(lines) == 0:
		try:
		    infile=open("ns_inputs","r")
		except:
		    exit_error("Failed to read ns_inputs file\n", 1)
		lines = infile.readlines()
	    args={}
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
	    args = comm.bcast(args,root=0)

	# parse args
	ns_args={}

	# convert from strings to actual args
	try:
	    ns_args['n_walkers'] = int(args.pop('n_walkers'))
	except:
	    exit_error("need number of walkers n_walkers\n",1)

	ns_args['n_cull'] = int(args.pop('n_cull', 1))

	try:
	    ns_args['n_iter_per_walker'] = int(args.pop('n_iter_per_walker'))
	except:
	    exit_error("need number of iterations n_iter_per_walker\n",1)
	ns_args['n_iter'] = ns_args['n_iter_per_walker']*ns_args['n_walkers']

	try:
	    ns_args['min_Emax'] = float(args.pop('min_Emax'))
	except:
	    ns_args['min_Emax'] = None

	ns_args['restart_file'] = args.pop('restart_file', '')
	if ns_args['restart_file'] == '':
	    try:
		ns_args['start_species'] = args.pop('start_species')
	    except:
		exit_error("need species for initial configs start_species\n",1)

	ns_args['max_volume_per_atom'] = float(args.pop('max_volume_per_atom', 1.0e3))

	ns_args['out_file_prefix'] = args.pop('out_file_prefix', '')
	if ns_args['out_file_prefix'] != '':
	    ns_args['out_file_prefix'] += '.'
	ns_args['profile'] = int(args.pop('profile', -1))
	ns_args['debug'] = int(args.pop('debug', -1))
	ns_args['snapshot_interval'] = int(args.pop('snapshot_interval', 1000))
	ns_args['traj_interval'] = int(args.pop('traj_interval', 1))
	ns_args['delta_random_seed'] = int(args.pop('delta_random_seed', -1))
	ns_args['n_extra_walk_per_task'] = int(args.pop('n_extra_walk_per_task', 0))
	ns_args['random_energy_perturbation'] = float(args.pop('random_energy_perturbation', 1.0e-12))
	ns_args['n_extra_data'] = int(args.pop('n_extra_data', 0))

	ns_args['start_energy_ceiling'] = float(args.pop('start_energy_ceiling', 1.0e9))
	ns_args['KEmax_max_T'] = float(args.pop('KEmax_max_T', 1.0e5))
	kB = 8.6173324e-5

	# parse energy_calculator
	ns_args['energy_calculator'] = args.pop('energy_calculator', 'fortran')
	do_calc_quip = False
	do_calc_lammps = False
	do_calc_internal = False
	do_calc_fortran = False
	if ns_args['energy_calculator'] == 'quip':
	    if not have_quippy:
		exit_error("Got energy_calculator=quip but not quippy module\n", 3)
	    do_calc_quip=True
	    try:
		ns_args['QUIP_pot_args'] = args.pop('QUIP_pot_args')
	    except:
		exit_error("need QUIP potential args QUIP_pot_args\n",1)
	    try:
		ns_args['QUIP_pot_params_file'] = args.pop('QUIP_pot_params_file')
	    except:
		exit_error("need QUIP potential params file QUIP_pot_params_file\n",1)
	elif ns_args['energy_calculator'] == 'lammps':
	    try:
		from lammpslib import LAMMPSlib
	    except:
		exit_error("energy_calculator=lammps and failed to import lammpslib module\n", 1)
	    do_calc_lammps=True
	    try:
		ns_args['LAMMPS_init_cmds'] = args.pop('LAMMPS_init_cmds')
	    except:
		exit_error("need LAMMPS initialization commands LAMMPS_init_cmds\n",1)
	    ns_args['LAMMPS_name'] = args.pop('LAMMPS_name', '')
	    ns_args['LAMMPS_header'] = args.pop('LAMMPS_header', 'units metal; atom_style atomic; atom_modify map array sort 0 0')
	    ns_args['LAMMPS_header_extra'] = args.pop('LAMMPS_header_extra', '')
            ns_args['LAMMPS_atom_types'] = None
	    LAMMPS_atom_types = args.pop('LAMMPS_atom_types', '')
            if len(LAMMPS_atom_types) > 0:
                ns_args['LAMMPS_atom_types'] = {}
                for type_pair in [s.strip() for s in LAMMPS_atom_types.split(',')]:
                    f = type_pair.split()
                    ns_args['LAMMPS_atom_types'][f[0]] = f[1]
	elif ns_args['energy_calculator'] == 'internal':
	    do_calc_internal=True
	elif ns_args['energy_calculator'] == 'fortran':
	    import fortranMCMDpy
	    do_calc_fortran=True
	    try:
		ns_args['FORTRAN_model'] = args.pop('FORTRAN_model')
	    except:
		exit_error("need FORTRAN model FORTRAN_model\n",1)
	    f_MC_MD = fortranMCMDpy.fortran_MC_MD(ns_args['FORTRAN_model'])
	    f_MC_MD.init_model()
	else:
	    exit_error("energy_calculator=%s unknown\n" % ns_args['energy_calculator'], 3)

	ns_args['no_extra_walks_at_all'] = str_to_logical(args.pop('no_extra_walks_at_all', "F"))

	ns_args['track_configs'] = str_to_logical(args.pop('track_configs', "F"))

	ns_args['config_file_format'] = args.pop('config_file_format', 'extxyz')

	ns_args['rng'] = args.pop('rng', 'numpy')

	if ns_args['rng'] == 'numpy':
	    rng = ns_rng.NsRngNumpy(ns_args['delta_random_seed'],comm)
	# elif ns_args['rng'] == 'julia':
	#    import julia
	#    j = julia.Julia()
	#    rng = ns_rng.NsRngJulia(j)
	elif ns_args['rng'] == 'rngstream':
	    import rngstream
	    rng = ns_rng.NsRngStream(ns_args['delta_random_seed'],comm)
	elif ns_args['rng'] == 'internal':
	    rng = ns_rng.NsRngInternal(ns_args['delta_random_seed'],comm)
	else:
	    exit_error("rng=%s unknown\n" % ns_args['rng'], 3)

	if do_calc_fortran:
	    l_seed = f_MC_MD.seed_size()
	    seed = np.array( [0] * l_seed , dtype=np.int32)
	    for i in range(l_seed):
		seed[i] = rng.int_uniform(1,sys.maxint)
	    f_MC_MD.set_seed(seed)

	movement_args={}

	movement_args['n_steps_expected'] = int(args.pop('n_steps_expected', 0))
	movement_args['n_steps'] = int(args.pop('n_steps', 0))
	movement_args['n_steps_unblocked_atom'] = int(args.pop('n_steps_unblocked_atom', 0))
	movement_args['n_steps_unblocked_cell_volume'] = int(args.pop('n_steps_unblocked_cell_volume', 0))
	movement_args['n_steps_unblocked_cell_shear'] = int(args.pop('n_steps_unblocked_cell_shear', 0))
	movement_args['n_steps_unblocked_cell_stretch'] = int(args.pop('n_steps_unblocked_cell_stretch', 0))
	movement_args['n_steps_unblocked_swap'] = int(args.pop('n_steps_unblocked_swap', 0))
	if (movement_args['n_steps_expected'] <= 0 and
	    movement_args['n_steps'] <= 0 and
	    movement_args['n_steps_unblocked_atom'] <= 0 and
	    movement_args['n_steps_unblocked_cell_volume'] <= 0 and
	    movement_args['n_steps_unblocked_cell_shear'] <= 0 and
	    movement_args['n_steps_unblocked_cell_stretch'] <= 0 and
	    movement_args['n_steps_unblocked_swap'] <= 0):
	    exit_error("Got all of n_steps_* == 0\n", 3)

	movement_args['atom_n_substeps'] = int(args.pop('atom_n_substeps', 10))
	movement_args['cell_n_substeps'] = int(args.pop('cell_n_substeps', 1))
	movement_args['swap_n_substeps'] = int(args.pop('swap_n_substeps', 0))
	movement_args['velo_n_substeps'] = int(args.pop('velo_n_substeps', 8))
	try:
	    movement_args['atom_algorithm'] = args.pop('atom_algorithm')
	except:
	    exit_error("Failed to read algorithm for atom motion atom_algorithm", 1)

	if movement_args['atom_algorithm'] != 'MC' and movement_args['atom_algorithm'] != 'MD':
	    exit_error("Got unknown atom_algorithm '%s'\n" % movement_args['atom_algorith,'], 3)

	movement_args['MC_atom_velocities'] = str_to_logical(args.pop('MC_atom_velocities', "F"))
	movement_args['MC_atom_velocities_pre_perturb'] = str_to_logical(args.pop('MC_atom_velocities_pre_perturb', "F"))
	movement_args['MC_atom_step_size'] = float(args.pop('MC_atom_step_size', 0.1))
	movement_args['MC_atom_step_size_max'] = float(args.pop('MC_atom_step_size_max', 0.5))
	movement_args['MC_atom_velo_step_size'] = float(args.pop('MC_atom_velo_step_size', 0.1))
	movement_args['MC_atom_velo_step_size_max'] = float(args.pop('MC_atom_velo_step_size_max', 1.0))
	movement_args['MC_atom_uniform_rv'] = str_to_logical(args.pop('MC_atom_uniform_rv', "F"))
	movement_args['do_velocities'] = (movement_args['atom_algorithm'] == 'MD' or movement_args['MC_atom_velocities'])

	movement_args['atom_velo_random_order_perturb'] = str_to_logical(args.pop('atom_velo_random_order_perturb', "F"))

	movement_args['MD_atom_velo_pre_perturb'] = str_to_logical(args.pop('MD_atom_velo_pre_perturb', "F"))
	movement_args['MD_atom_velo_post_perturb'] = str_to_logical(args.pop('MD_atom_velo_post_perturb', "T"))
	movement_args['MD_atom_velo_flip_accept'] = str_to_logical(args.pop('MD_atom_velo_flip_accept', "F"))
	movement_args['atom_velo_rej_free_perturb_randomize'] = str_to_logical(args.pop('atom_velo_rej_free_perturb_randomize', "F"))
	movement_args['atom_velo_rej_free_perturb_angle'] = float(args.pop('atom_velo_rej_free_perturb_angle', 0.3))
	movement_args['MC_atom_velo_walk_rej_free'] = str_to_logical(args.pop('MC_atom_velo_walk_rej_free', "T"))

	movement_args['MD_atom_timestep'] = float(args.pop('MD_atom_timestep', 0.1))
	movement_args['MD_atom_timestep_max'] = float(args.pop('MD_atom_timestep_max', 0.2))
	movement_args['MD_atom_energy_fuzz'] = float(args.pop('MD_atom_energy_fuzz', 1.0e-2))
	movement_args['MD_atom_reject_energy_violation'] = str_to_logical(args.pop('MD_atom_reject_energy_violation', "F"))

	movement_args['MC_cell_P'] = float(args.pop('MC_cell_P', 0.0))

	default_value = ns_args['max_volume_per_atom']/20.0 # 5% of maximum allowed volume per atom
	movement_args['MC_cell_volume_per_atom_step_size'] = float(args.pop('MC_cell_volume_per_atom_step_size', default_value))
	movement_args['MC_cell_volume_per_atom_step_size_max'] = float(args.pop('MC_cell_volume_per_atom_step_size_max', 5000.0))
	movement_args['MC_cell_volume_per_atom_prob'] = float(args.pop('MC_cell_volume_per_atom_prob', 1.0))
	movement_args['MC_cell_stretch_step_size'] = float(args.pop('MC_cell_stretch_step_size', 0.01))
	movement_args['MC_cell_stretch_step_size_max'] = float(args.pop('MC_cell_stretch_step_size_max', 0.05))
	movement_args['MC_cell_stretch_prob'] = float(args.pop('MC_cell_stretch_prob', 1.0))
	movement_args['MC_cell_shear_step_size'] = float(args.pop('MC_cell_shear_step_size', 0.01))
	movement_args['MC_cell_shear_step_size_max'] = float(args.pop('MC_cell_shear_step_size_max', 0.05))
	movement_args['MC_cell_shear_prob'] = float(args.pop('MC_cell_shear_prob', 1.0))

	movement_args['MC_cell_min_aspect_ratio'] = float(args.pop('MC_cell_min_aspect_ratio', 0.9))
	movement_args['cell_shape_equil_steps'] = int(args.pop('cell_shape_equil_steps', 100))

	movement_args['adjust_step_interval_per_walker'] = float(args.pop('adjust_step_interval_per_walker', 0.25))
	movement_args['adjust_step_interval'] = int(movement_args['adjust_step_interval_per_walker']*ns_args['n_walkers'])
	if movement_args['adjust_step_interval'] < 20:
	    print "WARNING: step size adjustment would be done too often, at every ", movement_args['adjust_step_interval'], " iteration"
	    print "WARNING: adjust_step_interval is increased to 20"
	    movement_args['adjust_step_interval'] = 20
	movement_args['MC_adjust_step_factor'] = float(args.pop('MC_adjust_step_factor', 1.5))
	movement_args['MC_adjust_min_rate'] = float(args.pop('MC_adjust_min_rate', 0.25))
	movement_args['MC_adjust_max_rate'] = float(args.pop('MC_adjust_max_rate', 0.75))
	movement_args['MD_adjust_step_factor'] = float(args.pop('MD_adjust_step_factor', 1.5))
	movement_args['MD_adjust_min_rate'] = float(args.pop('MD_adjust_min_rate', 0.95))
	movement_args['MD_adjust_max_rate'] = float(args.pop('MD_adjust_max_rate', 1.00))

	movement_args['2D'] = str_to_logical(args.pop('2D', "F"))

	if 'QUIP_pot_params_file' in ns_args:
	    if not have_quippy:
		exit_error("Got QUIP_pot_params but no quippy module\n", 3)
	    try:
		if rank == 0:
		    ns_args['QUIP_pot_params'] = open(ns_args['QUIP_pot_params_file'],"r").read()
		else:
		    ns_args['QUIP_pot_params'] = None
		if comm is not None:
		    ns_args['QUIP_pot_params'] = comm.bcast(ns_args['QUIP_pot_params'], root=0)
	    except:
		exit_error("Failed to read params file '%s'\n" % ns_args['QUIP_pot_params_file'], 1)

	if len(args) > 0:
	    exit_error(str(args)+"\nUnknown arguments read in\n", 2)

	if rank == 0:
	    print "ns_args ",ns_args
	    print "movement_args ",movement_args

	# initialise potential
	if do_calc_quip:
	    pot = quippy.Potential(ns_args['QUIP_pot_args'], param_str=ns_args['QUIP_pot_params'], calculation_always_required=True, cutoff_skin=1.0)
	elif do_calc_internal or do_calc_fortran:
	    pass
	elif do_calc_lammps:
	    init_cmds = [s.strip() for s in ns_args['LAMMPS_init_cmds'].split(';')]
	    header_cmds = [s.strip() for s in ns_args['LAMMPS_header'].split(';')]
	    header_extra_cmds = [s.strip() for s in ns_args['LAMMPS_header_extra'].split(';')]
	    if ns_args['debug'] >= 5:
		pot = LAMMPSlib(lmpcmds=init_cmds, atom_types=ns_args['LAMMPS_atom_types'], log_file='lammps.%d.log' % rank, keep_alive=True, lammps_name=ns_args['LAMMPS_name'],
				lammps_header=header_cmds, lammps_header_extra=header_extra_cmds)
	    else:
		pot = LAMMPSlib(lmpcmds=init_cmds, atom_types=ns_args['LAMMPS_atom_types'], keep_alive=True, lammps_name=ns_args['LAMMPS_name'],
				lammps_header=header_cmds, lammps_header_extra=header_extra_cmds)
	    print "PRE START_LAMMPS"
	    pot.start_lammps() # so that top level things like units will be set
	    print "POST START_LAMMPS"
	    pot.first_MD=True
	else:
	    exit_error("Need some way of initializing calculator\n",3)

	# figure out numbers of local walkers
	rank_of_walker = [0]*ns_args['n_walkers']
	if size <= 1:
	    n_walkers = ns_args['n_walkers']
	else:
	    n_walkers_per_task = ns_args['n_walkers']/size
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

	# internal model, LJ eps=1, sigma=1, cutoff=3,  with PBC cube l = pbc[0,0]
	internal_cutoff = 3.0
	Eshift = internal_cutoff**-12 - internal_cutoff**-6

	set_n_steps('n_steps')
	if  rank == 0:
	    print "Using n_steps = ", movement_args['n_steps']

	walkers=[]
	if ns_args['restart_file'] == '': # start from scratch
            start_first_iter = 0
	    # create initial config
	    if rank == 0:
		# create atoms structs from a list of atomic numbers and numbers of atoms
		lc = ns_args['max_volume_per_atom']**(1.0/3.0)
		init_atoms = ase.Atoms(cell=(lc, lc, lc), pbc=(1,1,1))
		species = ns_args['start_species'].split(',')
		for specie in species:
		    specie_fields = specie.split()
		    type_Z = int(specie_fields[0])
		    type_n = int(specie_fields[1])
		    if len(specie_fields) == 2:
			init_atoms += ase.Atoms([type_Z] * type_n)
		    elif len(specie_fields) == 3:
			type_mass = float(specie_fields[2])
			init_atoms += ase.Atoms([type_Z] * type_n, masses=[type_mass] * type_n)
		    else:
			exit_error("Each entry is start_species must include atomic number, multiplicity, and optionally mass", 5)
		init_atoms.set_cell(init_atoms.get_cell()*float(len(init_atoms))**(1.0/3.0), scale_atoms=True)
		if have_quippy:
		    init_atoms = quippy.Atoms(init_atoms)
		ase.io.write(sys.stdout, init_atoms, format=ns_args['config_file_format'])
	    else:
		init_atoms = None

	    if comm is not None:
		init_atoms = comm.bcast(init_atoms, root=0)
	    if do_calc_quip:
		init_atoms.set_cutoff(pot.cutoff(), cutoff_skin=1.0)

	    # make sure masses are set if MD is going to be used
	    if movement_args['do_velocities'] and not init_atoms.has('masses'):
		if have_quippy:
		    if not hasattr(init_atoms,'mass'):
			init_atoms.add_property('mass', 0.0)
			init_atoms.mass[:] = quippy.ElementMass[init_atoms.Z]
		    init_atoms.set_masses(init_atoms.mass/quippy.MASSCONVERT)
		else:
		    exit_error("MD set, but masses weren't specified in start_species, and quippy is not available for automatically setting masses\n", 3)
	    # create extra data arrays if needed
	    if ns_args['n_extra_data'] > 0:
		init_atoms.arrays['ns_extra_data'] = np.zeros( (len(init_atoms), ns_args['n_extra_data']) )

	    # clone initial config into array of walkers
	    for i_walker in range(n_walkers):
		walkers.append(init_atoms.copy())

            if ns_args['track_configs']:
                if comm is None:
                    config_ind = 0
                else:
                    config_ind = comm.rank*n_walkers
	    for at in walkers:
		at.set_velocities(np.zeros( (len(walkers[0]), 3) ))
                if ns_args['track_configs']:
                    at.info['config_ind'] = config_ind
                    at.info['from_config_ind'] = -1
                    at.info['config_ind_time'] = -1
                    config_ind += 1
		if do_calc_quip or do_calc_lammps:
		    at.set_calculator(pot)

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
	    ns_args['start_energy_ceiling'] += movement_args['MC_cell_P']*ns_args['max_volume_per_atom']*len(init_atoms)
	    # initial positions are just random, up to an energy ceiling
	    for at in walkers:
		# randomize cell if P is set, both volume (from appropriate distribution) and shape (from uniform distribution with min aspect ratio limit)
		if movement_args['MC_cell_P'] > 0.0:
		    lc = (len(at)*ns_args['max_volume_per_atom']*rng.float_uniform(0.0,1.0)**(1.0/float(len(at)+1)))**(1.0/3.0)
		    at.set_cell( np.identity(3) * lc )
		    do_cell_shape_walk(at, movement_args)

		# random initial positions
		energy = float('nan')
		n_try = 0
		while n_try < 10 and math.isnan(energy) or energy > ns_args['start_energy_ceiling']:
		    at.set_scaled_positions( rng.float_uniform(0.0, 1.0, (len(at), 3) ) )
		    energy = eval_energy(at)
		    n_try += 1
		if math.isnan(energy) or energy > ns_args['start_energy_ceiling']:
		    sys.stderr.write("WARNING: rank %d failed to generate initial config by random positions under max energy %f in 10 tries\n" % (rank, ns_args['start_energy_ceiling']), 4)

		# try FORTRAN config initializer
		n_try = 0
		if do_calc_fortran:
		    while n_try < 10 and math.isnan(energy) or energy > ns_args['start_energy_ceiling']:
			f_MC_MD.init_config(at, ns_args['start_energy_ceiling']-movement_args['MC_cell_P']*ns_args['max_volume_per_atom']*len(init_atoms))
			energy = eval_energy(at)
			n_try += 1
		    if math.isnan(energy) or energy > ns_args['start_energy_ceiling']:
			sys.stderr.write("WARNING: rank %d failed to generate initial config by fortran config initializer under max energy %f in 10 tries\n" % (rank, ns_args['start_energy_ceiling']), 4)

		# try python config initializer
		n_try = 0
		while n_try < 10 and math.isnan(energy) or energy > ns_args['start_energy_ceiling']:
		    energy = additive_init_config(at, ns_args['start_energy_ceiling'])
		    n_try += 1

		# quit if failed to generate acceptable config
		if math.isnan(energy) or energy > ns_args['start_energy_ceiling']:
		    exit_error("Rank %d failed to generate initial config by random, fortran, or (atom by atom addition) python initializer under max energy %f in 10 tries each\n" % (rank, ns_args['start_energy_ceiling']), 4)

		at.info['ns_energy'] = rand_perturb_energy(energy, ns_args['random_energy_perturbation'])
		at.info['volume'] = at.get_volume()

		# set KEmax from P and Vmax
		if movement_args['do_velocities']:
		    if movement_args['MC_cell_P'] > 0.0:
			KEmax = median_PV(walkers)
		    else:
			KEmax = kB*ns_args['KEmax_max_T']
                    at.info['KEmax']=KEmax
                 
		else:
		    KEmax = None

		# set initial velocities, rejection free
		if movement_args['do_velocities']:
		    rej_free_perturb_velo(at, None, KEmax)

	else: # doing a restart
	    if rank == 0: # read on head task and send to other tasks
		i_at = 0
		for r in range(size):
		    at_list=[]
		    for i_walker in range(n_walkers):
			print "read walker ",i_at, " from ",ns_args['restart_file']
			if comm is None or r == 0:
			    walkers.append(ase.io.read(ns_args['restart_file'], i_at))
			else:
			    at_list.append(ase.io.read(ns_args['restart_file'], i_at))
			i_at += 1
		    if r > 0:
			comm.send(at_list, dest=r, tag=1)
	    else: # receive from head task
		walkers = comm.recv(source=0, tag=1)

	    for at in walkers:
		if ns_args['n_extra_data'] > 0 and (not 'ns_extra_data' in at.arrays or at.arrays['ns_extra_data'].size/len(at) != ns_args['n_extra_data']):
		    at.arrays['ns_extra_data'] = np.zeros( (len(at), ns_args['n_extra_data']) )
		if do_calc_quip or do_calc_lammps:
		    at.set_calculator(pot)
		at.info['ns_energy'] = rand_perturb_energy(eval_energy(at), ns_args['random_energy_perturbation'])
                KEmax = at.info['KEmax']

		key_found = False
		for key in at.info: # check if 'iter=' info is present in the file used for restart
		    if key == 'iter':
                        start_first_iter = at.info['iter']+1
			key_found = True
		if not key_found:
		    print "WARNING: no iteration number information was found in the restart file"
                    exit_error("no iteration number information was found in the restart file\n",5)

		key_found = False
		for key in at.info: # check if 'volume=' info is present in the file used for restart
		    if key == 'volume':
		        movement_args['MC_cell_volume_per_atom_step_size'] = at.info['volume']/10.0/at.get_number_of_atoms()
			key_found = True
		if not key_found:
		    print "WARNING: no volume information was found in the restart file. If volume changes will be done, the starting stepsize will be the default"
	            

	    if have_quippy:
                walkers = [quippy.Atoms(at) for at in walkers]

	# scale MC_atom_step_size by max_vol^(1/3)
	max_lc = (ns_args['max_volume_per_atom']*len(walkers[0]))**(1.0/3.0)
	movement_args['MC_atom_step_size'] *= max_lc
	movement_args['MC_atom_step_size_max'] *= max_lc

	n_atoms = len(walkers[0])
	prev_snapshot_iter = None
	# do NS
	#QUIP_IO if have_quippy:
	    #QUIP_IO traj_io = quippy.CInOutput(ns_args['out_file_prefix']+'traj.%d.extxyz' % rank, action=quippy.OUTPUT)
	#QUIP_IO else:
	    #QUIP_IO traj_file = ns_args['out_file_prefix']+'traj.%08d.'+('%04d' % rank)+'.extxyz'

        # open the file where the trajectory will be printed
	if ns_args['restart_file'] == '': # start from scratch, so if this file exists, overwrite it 
            traj_io = open(ns_args['out_file_prefix']+'traj.%d.%s' % (rank, ns_args['config_file_format']), "w")
            if ns_args['track_configs']:
                track_traj_io = open(ns_args['out_file_prefix']+'track_traj.%d.%s' % (rank, ns_args['config_file_format']), "w")
            else:
                track_traj_io = None
        else: # restart, so the existing file should be appended
            traj_io = open(ns_args['out_file_prefix']+'traj.%d.%s' % (rank, ns_args['config_file_format']), "a")
            if ns_args['track_configs']:
                track_traj_io = open(ns_args['out_file_prefix']+'track_traj.%d.%s' % (rank, ns_args['config_file_format']), "a")
            else:
                track_traj_io = None

            # Read the existing traj file and look for the point where we restart from. Truncate the rest. 
	    # This part is not used because the ASE.io.read takes soooo long, that it makes a restart impossible.
            #traj_io = open(ns_args['out_file_prefix']+'traj.%d.%s' % (rank, ns_args['config_file_format']), "r+")
	    #i = 0
	    #while True:
            #    at=(ase.io.read(traj_io, format=ns_args['config_file_format'],index=i))
	    #	print "ASE.io.read trajectory", rank, i, at.info['iter']
            #    if at.info['iter'] >= start_first_iter:
	    #         at=(ase.io.read(traj_io, format=ns_args['config_file_format'],index=i-1))
	    #         traj_io.truncate()
            #         break
	    #    i += 1

        # open the file where the energies will be printed
	if rank == 0:
	    if ns_args['restart_file'] == '': # start from scratch, so if this file exists, overwrite it 
	        energy_io = open(ns_args['out_file_prefix']+'energies', 'w')
            else: # restart, so the existing file should be appended
	        energy_io = open(ns_args['out_file_prefix']+'energies', 'r+')
		tmp_iter = 0
                line = energy_io.readline() # read the first line of nwalker,ncull..etc information
		i = 0
		while True: # we do create an infinite loop here :(
		    line=energy_io.readline()              # read lines one by one
	  	    if not line:                           # something went wrong, exit the infinit loop
                        print "WARNING: end of .energies file reached without finding the iteration number", start_first_iter
			break
		    i = i+1
		    if i%n_cull==0:                        # if this is n_cull-th line, examine the stored iteration
                        tmp_split = line.split()
                        tmp_iter = int(tmp_split[0])       # tmp_iter contains the iteration number of the line as an integer number
                    if tmp_iter == start_first_iter-1:     # if this is the iteration same as in the snapshot, 
                        energy_io.truncate()                #delete the rest of the file, as we are restarting from here
                        break

	if ns_args['profile'] == rank:
	    import cProfile
	    cProfile.run('do_ns_loop()',ns_args['out_file_prefix']+'profile.stats')
	else:
	    do_ns_loop()

	# cleanup post loop
	save_snapshot(ns_args['n_iter']-1) # this is the final configuration
	#clean_prev_snapshot(prev_snapshot_iter) # !!!! This line is commented out as sometimes it deleted the final snapshot - check the prev_snashot_iter vale in case of normal termination!!!!!!!!!!!!

	for at in walkers:
	    print rank, ": final energy ", at.info['ns_energy']

	if rank == 0:
	    energy_io.close()
	traj_io.close()
        if track_traj_io is not None:
            track_traj_io.close()

	if comm is not None:
	    MPI.Finalize()
        sys.exit(0)

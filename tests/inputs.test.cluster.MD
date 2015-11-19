max_volume_per_atom=108.0
start_species=1 6 1.0
start_energy_ceiling=40.0

n_walkers=1024
n_cull=1
n_iter_times_fraction_killed=120
out_file_prefix=test.cluster.MD

################################################################################

atom_algorithm=MD

n_model_calls_expected=640
n_atom_steps=1
atom_traj_len=10
n_cell_volume_steps=0
n_cell_shear_steps=0
n_cell_stretch_steps=0

MD_atom_timestep=0.01
MD_atom_timestep_max=0.025

FORTRAN_model=example_LJ_min_image_model.so
energy_calculator=fortran

rng=numpy

debug=0

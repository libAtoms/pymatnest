max_volume_per_atom=108.0
start_species=1 6
start_energy_ceiling=1.0e9

n_walkers=1024
n_cull=1
n_iter_times_fraction_killed=100
out_file_prefix=test.cluster.MC

################################################################################

atom_algorithm=MC
n_model_calls_expected=256
n_atom_steps=1
atom_traj_len=8
n_cell_volume_steps=0
n_cell_shear_steps=0
n_cell_stretch_steps=0


MC_atom_step_size=5.0
MC_atom_step_size_max=5.0

FORTRAN_model=example_LJ_min_image_model.so
energy_calculator=fortran

rng=numpy

debug=0

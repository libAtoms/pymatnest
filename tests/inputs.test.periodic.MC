max_volume_per_atom=2.5e6

start_species=1 64
start_energy_ceiling=1.0e9

n_walkers=640
n_cull=1
n_iter_times_fraction_killed=3000
out_file_prefix=test.periodic.MC

################################################################################

atom_algorithm=MC
n_model_calls_expected=640

n_atom_steps=1
atom_traj_len=8
n_cell_volume_steps=4
n_cell_shear_steps=4
n_cell_stretch_steps=4


MC_atom_step_size=5.0
MC_atom_step_size_max=5.0

################################################################################
MC_cell_P=2.37e-4

MC_cell_volume_per_atom_step_size=100.0
MC_cell_volume_per_atom_step_size_max=5000.0

MC_cell_shear_step_size=0.01
MC_cell_shear_step_size_max=0.1

MC_cell_stretch_step_size=0.01
MC_cell_stretch_step_size_max=0.1

MC_cell_min_aspect_ratio=0.9
################################################################################

FORTRAN_model=example_LJ_model.so
energy_calculator=fortran

rng=numpy

debug=0

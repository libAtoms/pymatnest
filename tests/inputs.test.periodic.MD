max_volume_per_atom=2.5e6

start_species=1 64 1.0
start_energy_ceiling=100.0

n_walkers=640
n_cull=1
n_iter_times_fraction_killed=5000
out_file_prefix=test.periodic.MD

################################################################################

atom_algorithm=MD
n_model_calls_expected=1408


n_atom_steps=1
atom_traj_len=32
n_cell_volume_steps=4
n_cell_stretch_steps=4
n_cell_shear_steps=4

MD_atom_timestep=0.01
MD_atom_timestep_max=0.025

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

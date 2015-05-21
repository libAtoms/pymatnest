start_volume_per_atom=108.13333

start_species=1 6 1.0
start_energy_ceiling=1.0e9

n_walkers=128
n_cull=1
n_iter_per_walker=10
out_file_prefix=cluster.MC

atom_algorithm=MC

atom_n_steps_expected=32
atom_n_substeps=8


adjust_step_interval=-10

MC_atom_step_size=5.0
MC_atom_step_size_max=5.0

FORTRAN_model=example_LJ_min_image_model.so
energy_calculator=fortran

rng=numpy

debug=0

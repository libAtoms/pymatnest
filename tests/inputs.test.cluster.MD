start_volume_per_atom=108.0
start_species=1 6 1.0
start_energy_ceiling=40.0

n_walkers=1024
n_cull=1
n_iter_per_walker=120
out_file_prefix=test.cluster.MD

atom_algorithm=MD

atom_n_steps_expected=64
atom_n_substeps=10
cell_n_substeps=0

MD_atom_timestep=0.01
MD_atom_timestep_max=0.025
MD_atom_velo_pre_perturb=T

FORTRAN_model=example_LJ_min_image_model.so
energy_calculator=fortran

rng=numpy

debug=0

start_volume_per_atom=108.133333
start_species=1 6 1.0
start_energy_ceiling=40.0

n_walkers=128
n_cull=1
n_iter_per_walker=10
out_file_prefix=cluster.MD

atom_algorithm=MD

atom_n_steps_expected=50
atom_n_substeps=10

adjust_step_interval=-10

MD_atom_timestep=0.01
MD_atom_timestep_max=0.01
MD_atom_velo_perturb=T
MD_atom_velo_randomize=F
MD_atom_velo_perturb_angle=0.1

FORTRAN_model=example_LJ_min_image_model.so
energy_calculator=fortran

rng=numpy

debug=0

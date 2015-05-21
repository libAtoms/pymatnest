start_volume_per_atom=1116

start_species=1 64 1.0
start_energy_ceiling=320.0

n_walkers=128
n_cull=1
n_iter_per_walker=1700
out_file_prefix=test.periodic.MD

################################################################################
atom_n_steps_expected=128
atom_n_substeps=8
cell_n_substeps=4

atom_algorithm=MD

adjust_step_interval=10

MD_atom_timestep=0.01
MD_atom_timestep_max=0.025
MD_atom_velo_perturb=T
MD_atom_velo_randomize=F
MD_atom_velo_perturb_angle=0.1

################################################################################
MC_cell_P=0.00675

MC_cell_volume_step_size=0.01
MC_cell_volume_step_size_max=0.1

MC_cell_shear_step_size=0.01
MC_cell_shear_step_size_max=0.1

MC_cell_stretch_step_size=0.01
MC_cell_stretch_step_size_max=0.1

MC_cell_min_aspect_ratio=0.85
MC_cell_max_volume_per_atom=1139.0
################################################################################

FORTRAN_model=example_LJ_model.so
energy_calculator=fortran

rng=numpy

debug=0

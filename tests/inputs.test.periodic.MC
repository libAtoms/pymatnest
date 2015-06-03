start_volume_per_atom=1116

start_species=1 64
start_energy_ceiling=1.0e9

n_walkers=128
n_cull=1
n_iter_per_walker=1200
out_file_prefix=test.periodic.MC

################################################################################

atom_algorithm=MC
atom_n_steps_expected=32
atom_n_substeps=8
cell_n_substeps=1


MC_atom_step_size=5.0
MC_atom_step_size_max=5.0

################################################################################
MC_cell_P=0.00675

MC_cell_volume_per_atom_step_size=100.0
MC_cell_volume_per_atom_step_size_max=5000.0

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

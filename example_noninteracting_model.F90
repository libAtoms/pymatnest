! API documented in example_LJ_model.F90
subroutine ll_init_model()
   return
end subroutine ll_init_model

subroutine ll_init_config()
   return
end subroutine ll_init_config

double precision function ll_eval_energy(N, pos, n_extra_data, extra_data, cell)
use mat_mod
implicit none
   integer :: N
   double precision :: pos(3,N), cell(3,3)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data, N)

   ll_eval_energy = 0.0

end function ll_eval_energy

integer function ll_move_atom_1(N, pos, n_extra_data, extra_data, cell, d_i, d_pos, dEmax, dE)
use mat_mod
implicit none
   integer :: N
   double precision :: pos(3,N), cell(3,3)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data, N)
   integer :: d_i
   double precision :: d_pos(3)
   double precision :: dEmax, dE

   ll_move_atom_1 = 1
   dE = 0.0

end function ll_move_atom_1

function ll_eval_forces(N, pos, n_extra_data, extra_data, cell, forces) result(energy)
use mat_mod
implicit none
   integer :: N
   double precision :: pos(3,N), cell(3,3), forces(3,N)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data, N)
   double precision :: energy ! result

   integer :: i, j
   double precision :: dr(3), dr_mag, dr_l(3), dr_l0(3), pos_l(3,N)
   double precision :: cell_inv(3,3)
   integer :: dj1, dj2, dj3

   double precision :: E_offset  = 1.0/3.0**12 - 1.0/3.0**6, E_term

   energy = 0.0
   forces = 0.0

end function ll_eval_forces

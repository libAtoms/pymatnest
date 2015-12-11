! API documented in example_LJ_model.F90

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ll_init_model()
use example_LJ_params_mod
implicit none
   epsilon(1,1) = 1.0
   epsilon(2,2) = 1.0
   epsilon(1,2) = 1.5
   epsilon(2,1) = epsilon(1,2)

   sigma(1,1) = 3.0
   sigma(2,2) = 3.3
   sigma(1,2) = (sigma(1,1)+sigma(2,2))/2.0
   sigma(2,1) = sigma(1,2)

   cutoff = 3.0*sigma
   cutoff_sq = cutoff*cutoff

   E_offset = (sigma/cutoff)**12 - (sigma/cutoff)**6

end subroutine ll_init_model

subroutine ll_init_config()
   return
end subroutine ll_init_config

double precision function ll_eval_energy(N, Z, pos, n_extra_data, extra_data, cell)
use example_mat_mod
use example_LJ_params_mod
implicit none
   integer :: N
   integer :: Z(N)
   double precision :: pos(3,N), cell(3,3)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data, N)

   integer :: i, j, Z_i, Z_j
   double precision :: dr(3), dr_mag, dr_l(3)

   double precision :: cell_inv(3,3)

   call matrix3x3_inverse(cell, cell_inv)

   ll_eval_energy = 0.0
   do i=1, N
   Z_i = Z(i)
   do j=i+1, N
      Z_j = Z(j)
      dr = pos(:,i)-pos(:,j)
      dr_l = matmul(cell_inv, dr)
      dr_l = dr_l - floor(dr_l+0.5)
      dr = matmul(cell, dr_l)
      dr_mag = sqrt(sum(dr*dr))
      if (dr_mag < cutoff(Z_i,Z_j)) then
	 ll_eval_energy = ll_eval_energy + epsilon(Z_i,Z_j)*((sigma(Z_i,Z_j)/dr_mag)**12 - &
                                                             (sigma(Z_i,Z_j)/dr_mag)**6 - E_offset(Z_i,Z_j))
      endif
   end do
   end do

end function ll_eval_energy

integer function ll_move_atom_1(N, Z, pos, n_extra_data, extra_data, cell, d_i, d_pos, dEmax, dE)
use example_mat_mod
use example_LJ_params_mod
implicit none
   integer :: N
   integer :: Z(N)
   double precision :: pos(3,N), cell(3,3)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data, N)
   integer :: d_i
   double precision :: d_pos(3)
   double precision :: dEmax, dE

   integer :: i, j, Z_i, Z_j
   double precision :: dr0(3), dr1(3), dr0_l(3), dr1_l(3), dr0_mag, dr1_mag
   double precision :: cell_inv(3,3)

   call matrix3x3_inverse(cell, cell_inv)

   dE = 0.0
   i=d_i
   Z_i = Z(i)
   do j=1,N
      if (j == i) cycle
      Z_j = Z(j)

      dr0 = pos(:,i) - pos(:,j)
      dr1 = (pos(:,i)+d_pos(:)) - pos(:,j)

      dr0_l = matmul(cell_inv, dr0)
      dr0_l = dr0_l - floor(dr0_l+0.5)
      dr0 = matmul(cell, dr0_l)
      dr0_mag = sqrt(sum(dr0*dr0))
      dr1_l = matmul(cell_inv, dr1)
      dr1_l = dr1_l - floor(dr1_l+0.5)
      dr1 = matmul(cell, dr1_l)
      dr1_mag = sqrt(sum(dr1*dr1))

      if (dr0_mag < cutoff(Z_i,Z_j)) then
	 dE = dE - epsilon(Z_i,Z_j)*((sigma(Z_i,Z_j)/dr0_mag)**12 - (sigma(Z_i,Z_j)/dr0_mag)**6 - E_offset(Z_i,Z_j))
      endif
      if (dr1_mag < cutoff(Z_i,Z_j)) then
	 dE = dE + epsilon(Z_i,Z_j)*((sigma(Z_i,Z_j)/dr1_mag)**12 - (sigma(Z_i,Z_j)/dr1_mag)**6 - E_offset(Z_i,Z_j))
      endif

   end do

   if (dE < dEmax) then ! accept
      pos(:,i) = pos(:,i) + d_pos(:)
      ll_move_atom_1 = 1
   else ! reject
      dE = 0.0
      ll_move_atom_1 = 0
   endif

end function ll_move_atom_1

function ll_eval_forces(N, Z, pos, n_extra_data, extra_data, cell, forces) result(energy)
use example_mat_mod
use example_LJ_params_mod
implicit none
   integer :: N
   integer :: Z(N)
   double precision :: pos(3,N), cell(3,3), forces(3,N)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data, N)
   double precision :: energy ! result

   integer :: i, j, Z_i, Z_j
   double precision :: dr(3), dr_mag, dr_l(3)
   double precision :: cell_inv(3,3)

   call matrix3x3_inverse(cell, cell_inv)

   energy = 0.0
   forces = 0.0
   do i=1, N
   Z_i = Z(i)
   do j=i+1, N
      Z_j = Z(j)
      dr = pos(:,i)-pos(:,j)
      dr_l = matmul(cell_inv, dr)
      dr_l = dr_l - floor(dr_l+0.5)
      dr = matmul(cell, dr_l)
      dr_mag = sqrt(sum(dr*dr))
      if (dr_mag < cutoff(Z_i,Z_j)) then
	 energy = energy + epsilon(Z_i,Z_j)*((sigma(Z_i,Z_j)/dr_mag)**12 - (sigma(Z_i,Z_j)/dr_mag**6) - E_offset(Z_i,Z_j))
	 forces(:,i) = forces(:,i) - epsilon(Z_i,Z_j)*(-12.0*sigma(Z_i,Z_j)**12/dr_mag**13 + 6.0*sigma(Z_i,Z_j)**6/dr_mag**7)*(dr/dr_mag)
	 forces(:,j) = forces(:,j) + epsilon(Z_i,Z_j)*(-12.0*sigma(Z_i,Z_j)**12/dr_mag**13 + 6.0*sigma(Z_i,Z_j)**6/dr_mag**7)*(dr/dr_mag)
      endif
   end do
   end do

end function ll_eval_forces

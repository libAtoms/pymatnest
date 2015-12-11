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

   integer :: i, j
   double precision :: dr(3), dr_mag, dr_l(3)
   double precision :: E_offset  = 1.0/3.0**12 - 1.0/3.0**6

   double precision :: cell_inv(3,3)

   call matrix3x3_inverse(cell, cell_inv)

   ll_eval_energy = 0.0
   do i=1, N
   do j=i+1, N
      dr = pos(:,i)-pos(:,j)
      dr_l = matmul(cell_inv, dr)
      dr_l = dr_l - floor(dr_l+0.5)
      dr = matmul(cell, dr_l)
      dr_mag = sqrt(sum(dr*dr))
      if (dr_mag < 3.0) then
	 ll_eval_energy = ll_eval_energy + ((1.0/dr_mag**12 - 1.0/dr_mag**6) - E_offset)
      endif
   end do
   end do

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

   double precision :: E_offset  = 1.0/3.0**12 - 1.0/3.0**6

   integer :: i, j
   double precision :: dr0(3), dr1(3), dr0_l(3), dr1_l(3), dr0_mag, dr1_mag
   double precision :: cell_inv(3,3)

   call matrix3x3_inverse(cell, cell_inv)

   dE = 0.0
   i=d_i
   do j=1,N
      if (j == i) cycle

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

      if (dr0_mag < 3.0) then
	 dE = dE - ((1.0/dr0_mag**12 - 1.0/dr0_mag**6) - E_offset)
      endif
      if (dr1_mag < 3.0) then
	 dE = dE + ((1.0/dr1_mag**12 - 1.0/dr1_mag**6) - E_offset)
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

function ll_eval_forces(N, pos, n_extra_data, extra_data, cell, forces) result(energy)
use mat_mod
implicit none
   integer :: N
   double precision :: pos(3,N), cell(3,3), forces(3,N)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data, N)
   double precision :: energy ! result

   integer :: i, j
   double precision :: dr(3), dr_mag, dr_l(3)
   double precision :: cell_inv(3,3)

   double precision :: E_offset  = 1.0/3.0**12 - 1.0/3.0**6

   call matrix3x3_inverse(cell, cell_inv)

   energy = 0.0
   forces = 0.0
   do i=1, N
   do j=i+1, N
      dr = pos(:,i)-pos(:,j)
      dr_l = matmul(cell_inv, dr)
      dr_l = dr_l - floor(dr_l+0.5)
      dr = matmul(cell, dr_l)
      dr_mag = sqrt(sum(dr*dr))
      if (dr_mag < 3.0) then
	 energy = energy + ((1.0/dr_mag**12 - 1.0/dr_mag**6) - E_offset)
	 forces(:,i) = forces(:,i) - (-12.0/dr_mag**13 + 6.0/dr_mag**7)*(dr/dr_mag)
	 forces(:,j) = forces(:,j) + (-12.0/dr_mag**13 + 6.0/dr_mag**7)*(dr/dr_mag)
      endif
   end do
   end do

end function ll_eval_forces

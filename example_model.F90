module mat_mod
implicit none
private

public :: matrix3x3_inverse

interface operator(.cross.)
  module procedure cross_product
end interface

contains
  !% Return the scalar triple product $\mathbf{x} \cdot \mathbf{y} \times \mathbf{z}$
  !% of the 3-vectors 'x', 'y' and 'z'.
  pure function scalar_triple_product(x,y,z)  ! [x,y,z]

    double precision, dimension(3), intent(in) :: x,y,z
    double precision                           :: scalar_triple_product

    scalar_triple_product = + x(1) * ( y(2)*z(3) - y(3)*z(2) )    &
         - x(2) * ( y(1)*z(3) - y(3)*z(1) )    &
         + x(3) * ( y(1)*z(2) - y(2)*z(1) )

  end function scalar_triple_product


  pure function cross_product(x,y) ! x ^ y

    double precision, dimension(3), intent(in):: x,y
    double precision, dimension(3)            :: cross_product

    cross_product(1) = + ( x(2)*y(3) - x(3)*y(2) )
    cross_product(2) = - ( x(1)*y(3) - x(3)*y(1) )
    cross_product(3) = + ( x(1)*y(2) - x(2)*y(1) )

  end function cross_product


   ! Matrix3x3_Inverse
   !
   !% Calculate $3\times3$ matrix inverse of 'lattice' and store result in 'g'.
   !% Avoids overhead of calling \textsc{lapack} for simple case of $3\times3$ matrix.
   subroutine matrix3x3_inverse(matrix, g)
     double precision, intent(in)  :: matrix(3,3)
     double precision, intent(out) :: g(3,3)

     double precision :: stp

     stp = scalar_triple_product(matrix(:,1), matrix(:,2), matrix(:,3))

     g(1,:) = (matrix(:,2) .cross. matrix(:,3))/stp
     g(2,:) = (matrix(:,3) .cross. matrix(:,1))/stp
     g(3,:) = (matrix(:,1) .cross. matrix(:,2))/stp

   end subroutine matrix3x3_inverse
end module mat_mod


double precision function ll_eval_energy(N, pos, cell)
use mat_mod
implicit none
   integer :: N
   double precision :: pos(3,N), cell(3,3)

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

double precision function ll_eval_energy_1(N, pos, cell, d_i, d_pos)
use mat_mod
implicit none
   integer :: N
   double precision :: pos(3,N), cell(3,3)
   integer :: d_i
   double precision :: d_pos(3)

   double precision :: E_offset  = 1.0/3.0**12 - 1.0/3.0**6

   integer :: i, j
   double precision :: dr0(3), dr1(3), dr0_l(3), dr1_l(3), dr0_mag, dr1_mag
   double precision :: cell_inv(3,3)

   call matrix3x3_inverse(cell, cell_inv)

   ll_eval_energy_1 = 0.0
   i=d_i
   do j=1,N
      if (j == i) cycle

      dr0 = pos(:,i)-pos(:,j)
      dr1 = d_pos(:)-pos(:,j)

      dr0_l = matmul(cell_inv, dr0)
      dr0_l = dr0_l - floor(dr0_l+0.5)
      dr0 = matmul(cell, dr0_l)
      dr0_mag = sqrt(sum(dr0*dr0))
      dr1_l = matmul(cell_inv, dr1)
      dr1_l = dr1_l - floor(dr1_l+0.5)
      dr1 = matmul(cell, dr1_l)
      dr1_mag = sqrt(sum(dr1*dr1))

      if (dr0_mag < 3.0) then
	 ll_eval_energy_1 = ll_eval_energy_1 - ((1.0/dr0_mag**12 - 1.0/dr0_mag**6) - E_offset)
      endif
      if (dr1_mag < 3.0) then
	 ll_eval_energy_1 = ll_eval_energy_1 + ((1.0/dr1_mag**12 - 1.0/dr1_mag**6) - E_offset)
      endif

   end do


end function ll_eval_energy_1

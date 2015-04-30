! publically accessible things required for interface to pymatnest
!
! subroutine ll_init_model() 
!    initializes potential
!
! subroutine ll_init_config(N, pos, cell, Emax) 
!    initializes a configuration with energy < Emax
!    config will be tested for failure after return
!
! double precision function ll_eval_energy(N, pos, cell)
!    integer :: N ! number of atoms
!    double precision :: pos(3,N), cell(3,3) ! positions, cell vectors
!    returns energy
!
! double precision function ll_eval_denergy_1(N, pos, cell, d_i, d_pos)
!    integer :: N ! number of atoms
!    double precision :: pos(3,N), cell(3,3) ! positions, cell vectors
!    integer :: d_i ! index of atom to be perturbed, 1-based (called from fortran_MC())
!    double precision :: d_pos(3) ! displacement of perturbed atom
!    returns energy change
!
! double precision function ll_eval_forces(N, pos, cell, forces)
!    integer :: N ! number of atoms
!    double precision :: pos(3,N), cell(3,3), forces(3,N) ! positions, cell vectors, forces
!    returns energy

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ll_init_model()
   return
end subroutine ll_init_model

subroutine ll_init_config()
   return
end subroutine ll_init_config

double precision function ll_eval_energy(N, pos, cell)
use mat_mod
implicit none
   integer :: N
   double precision :: pos(3,N), cell(3,3)

   integer :: i, j
   double precision :: dr(3), dr_mag, dr_l(3), dr_l0(3), pos_l(3,N)
   double precision :: E_offset  = 1.0/3.0**12 - 1.0/3.0**6

   double precision :: cell_inv(3,3), E_term
   integer :: dj1, dj2, dj3

   call matrix3x3_inverse(cell, cell_inv)
   ! into lattice coodinates 
   pos_l = matmul(cell_inv, pos)

   ll_eval_energy = 0.0
   do i=1, N
   do j=i, N
      dr_l0 = pos_l(:,i)-pos_l(:,j)
      dr_l0 = dr_l0 - floor(dr_l0+0.5)
      do dj1=-1,1
      dr_l(1) = dr_l0(1) + real(dj1, 8)
      do dj2=-1,1
      dr_l(2) = dr_l0(2) + real(dj2, 8)
      do dj3=-1,1
      dr_l(3) = dr_l0(3) + real(dj3, 8)
	 if (i == j .and. dj1 == 0 .and. dj2 == 0 .and. dj3 == 0) cycle

	 dr = matmul(cell, dr_l)
	 dr_mag = sqrt(sum(dr*dr))
	 if (dr_mag < 3.0) then
	    E_term = ((1.0/dr_mag**12 - 1.0/dr_mag**6) - E_offset)
	    if (i == j) E_term = E_term * 0.5
	    ll_eval_energy = ll_eval_energy + E_term
	 endif
      end do
      end do
      end do
   end do
   end do

end function ll_eval_energy

double precision function ll_eval_denergy_1(N, pos, cell, d_i, d_pos)
use mat_mod
implicit none
   integer :: N
   double precision :: pos(3,N), cell(3,3)
   integer :: d_i
   double precision :: d_pos(3)

   double precision :: E_offset  = 1.0/3.0**12 - 1.0/3.0**6

   integer :: i, j
   double precision :: dr(3), drp(3), dr_l(3), drp_l(3), dr_l0(3), drp_l0(3), dr_mag, drp_mag, pos_l(3,N)

   double precision :: cell_inv(3,3) 
   integer :: dj1, dj2, dj3

   call matrix3x3_inverse(cell, cell_inv)
   ! into lattice coodinates 
   pos_l = matmul(cell_inv, pos)

   ll_eval_denergy_1 = 0.0
   i=d_i
   do j=1,N
      if (j == i) cycle

      dr_l0 = pos_l(:,i) - pos_l(:,j)
      dr_l0 = dr_l0 - floor(dr_l0+0.5)
      do dj1=-1,1
      dr_l(1) = dr_l0(1) + real(dj1, 8)
      do dj2=-1,1
      dr_l(2) = dr_l0(2) + real(dj2, 8)
      do dj3=-1,1
      dr_l(3) = dr_l0(3) + real(dj3, 8)

	 dr = matmul(cell, dr_l)
	 drp = dr + d_pos
	 dr_mag = sqrt(sum(dr*dr))
	 drp_mag = sqrt(sum(drp*drp))

	 if (dr_mag < 3.0) then
	    ll_eval_denergy_1 = ll_eval_denergy_1 -  ((1.0/dr_mag**12 - 1.0/dr_mag**6) - E_offset)
	 endif
	 if (drp_mag < 3.0) then
	    ll_eval_denergy_1 = ll_eval_denergy_1 + ((1.0/drp_mag**12 - 1.0/drp_mag**6) - E_offset)
	 endif

      end do
      end do
      end do
   end do

end function ll_eval_denergy_1

function ll_eval_forces(N, pos, cell, forces) result(energy)
use mat_mod
implicit none
   integer :: N
   double precision :: pos(3,N), cell(3,3), forces(3,N)
   double precision :: energy ! result

   integer :: i, j
   double precision :: dr(3), dr_mag, dr_l(3), dr_l0(3), pos_l(3,N)
   double precision :: cell_inv(3,3)
   integer :: dj1, dj2, dj3

   double precision :: E_offset  = 1.0/3.0**12 - 1.0/3.0**6, E_term

   call matrix3x3_inverse(cell, cell_inv)
   pos_l = matmul(cell_inv, pos)

   energy = 0.0
   forces = 0.0
   do i=1, N
   do j=i, N

      dr_l0 = pos_l(:,i) - pos_l(:,j)
      dr_l0 = dr_l0 - floor(dr_l0+0.5)
      do dj1=-1,1
      dr_l(1) = dr_l0(1) + real(dj1, 8)
      do dj2=-1,1
      dr_l(2) = dr_l0(2) + real(dj2, 8)
      do dj3=-1,1
      dr_l(3) = dr_l0(3) + real(dj3, 8)
      if (i == j .and. dj1 == 0 .and. dj2 == 0 .and. dj3 == 0) cycle

	 dr = matmul(cell, dr_l)
	 dr_mag = sqrt(sum(dr*dr))
	 if (dr_mag < 3.0) then
	    E_term = ((1.0/dr_mag**12 - 1.0/dr_mag**6) - E_offset)
	    if (i == j) E_term = E_term * 0.5
	    energy = energy + E_term
	    if (i /= j) then
	       forces(:,i) = forces(:,i) - (-12.0/dr_mag**13 + 6.0/dr_mag**7)*(dr/dr_mag)
	       forces(:,j) = forces(:,j) + (-12.0/dr_mag**13 + 6.0/dr_mag**7)*(dr/dr_mag)
	    endif
	 endif

      end do
      end do
      end do
   end do
   end do

end function ll_eval_forces

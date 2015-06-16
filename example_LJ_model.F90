! publically accessible things required for interface to pymatnest
!
! subroutine ll_init_model() 
!
!    initializes potential
!
! subroutine ll_init_config(N, pos, cell, Emax) 
!    integer :: N ! number of atoms
!    double precision :: pos(3,N), cell(3,3) ! positions, cell vectors
!    double precision :: Emax ! maximum energy for config acceptance
!
!    initializes a configuration with energy < Emax
!    config will be tested for failure after return
!
! double precision function ll_eval_energy(N, pos, n_extra_data, extra_data, cell)
!    integer :: N ! number of atoms
!    double precision :: pos(3,N), cell(3,3) ! positions, cell vectors
!    integer :: n_extra_data ! width of extra data array
!    double precision :: extra_data(n_extra_data, N) ! extra data on output
!
!    evaluates energy of a config, sets extra_data, returns energy
!
! integer function ll_move_atom_1(N, pos, n_extra_data, extra_data, cell, d_i, d_pos, dEmax, dE)
!    integer :: N ! number of atoms
!    double precision :: pos(3,N), cell(3,3) ! positions, cell vectors, on output updated (pos only) to be consistent with acceptance/rejection
!    integer :: n_extra_data ! width of extra data array
!    double precision :: extra_data(n_extra_data, N) ! extra data on input, on output updated to be consistent with acceptance/rejection
!    integer :: d_i ! index of atom to be perturbed, 1-based (called from fortran_MC())
!    double precision :: d_pos(3) ! displacement of perturbed atom
!    double precision :: dEmax ! maximum change in energy for move acceptance
!    double precision :: dE ! on output actual change in energy, 0.0 if move is rejected
!
!    moves an atom if dE < dEmax
!    if move is accepted, updates pos, extra_data, sets dE
!    if move is rejected, nothing is updated, dE set to 0.0
!    returns 1 for accept, 0 for reject
!
! double precision function ll_eval_forces(N, pos, n_extra_data, extra_data, cell, forces)
!    integer :: N ! number of atoms
!    double precision :: pos(3,N), cell(3,3), forces(3,N) ! positions, cell vectors, forces
!    integer :: n_extra_data ! width of extra data array
!    double precision :: extra_data(n_extra_data, N) ! extra data on output
!
!    evaluates forces, sets extra_data
!    returns energy

module LJ_params_mod
   double precision, parameter :: epsilon = 0.4
   double precision, parameter :: sigma = 3.0
   double precision, parameter :: cutoff = 9.0, cutoff_sq = cutoff*cutoff
   double precision, parameter :: E_offset  = (sigma/cutoff)**12 - (sigma/cutoff)**6
end module LJ_params_mod

module mat_mod

implicit none
private

public :: matrix3x3_inverse

public operator(.cross.)
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

double precision function ll_eval_energy(N, pos, n_extra_data, extra_data, cell)
use mat_mod
use LJ_params_mod
implicit none
   integer :: N
   double precision :: pos(3,N), cell(3,3)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data, N)

   integer :: i, j
   double precision :: dr(3), dr_mag, dr_mag_sq, dr_l(3), dr_l0(3), pos_l(3,N)

   double precision :: cell_inv(3,3), E_term
   integer :: dj1, dj2, dj3

   integer :: n_images
   double precision cell_height(3), v_norm_hat(3)

   call matrix3x3_inverse(cell, cell_inv)
   ! into lattice coodinates 
   pos_l = matmul(cell_inv, pos)

   if (n_extra_data == 1) extra_data = 0.0

   do i=1, 3
      v_norm_hat = cell(:,mod(i,3)+1) .cross. cell(:,mod(i+1,3)+1)
      v_norm_hat = v_norm_hat / sqrt(sum(v_norm_hat**2))
      cell_height(i) = abs(sum(v_norm_hat*cell(:,i)))
   end do
   n_images = ceiling(cutoff/minval(cell_height))

   ll_eval_energy = 0.0
   do i=1, N
   do j=i, N
      dr_l0 = pos_l(:,i)-pos_l(:,j)
      dr_l0 = dr_l0 - floor(dr_l0+0.5)
      do dj1=-n_images,n_images
      dr_l(1) = dr_l0(1) + real(dj1, 8)
      do dj2=-n_images,n_images
      dr_l(2) = dr_l0(2) + real(dj2, 8)
      do dj3=-n_images,n_images
      dr_l(3) = dr_l0(3) + real(dj3, 8)
	 if (i == j .and. dj1 == 0 .and. dj2 == 0 .and. dj3 == 0) cycle

	 dr(1) = sum(cell(1,:)*dr_l)
	 dr(2) = sum(cell(2,:)*dr_l)
	 dr(3) = sum(cell(3,:)*dr_l)
	 dr_mag_sq = sum(dr*dr)
	 if (dr_mag_sq < cutoff_sq) then
	    dr_mag = sqrt(dr_mag_sq)
	    E_term = epsilon*(((sigma/dr_mag)**12 - (sigma/dr_mag)**6) - E_offset)
	    if (i == j) E_term = E_term * 0.5
	    ll_eval_energy = ll_eval_energy + E_term

	    if (n_extra_data == 1) then
	       extra_data(1,i) = extra_data(1,i) + 0.5*E_term
	       extra_data(1,j) = extra_data(1,j) + 0.5*E_term
	    endif

	 endif
      end do
      end do
      end do
   end do
   end do

end function ll_eval_energy

integer function ll_move_atom_1(N, pos, n_extra_data, extra_data, cell, d_i, d_pos, dEmax, dE)
use mat_mod
use LJ_params_mod
implicit none
   integer :: N
   double precision :: pos(3,N), cell(3,3)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data, N)
   integer :: d_i
   double precision :: d_pos(3)
   double precision :: dEmax, dE

   integer :: i, j
   double precision :: dr(3), drp(3), dr_l(3), drp_l(3), dr_l0(3), drp_l0(3), dr_mag, drp_mag, &
		       dr_mag_sq, drp_mag_sq, pos_l(3,N), d_pos_l(3)

   double precision :: cell_inv(3,3) 
   integer :: dj1, dj2, dj3

   double precision, allocatable, save :: new_extra_data(:,:)

   integer n_images
   double precision cell_height(3), v_norm_hat(3)

   do i=1, 3
      v_norm_hat = cell(:,mod(i,3)+1) .cross. cell(:,mod(i+1,3)+1)
      v_norm_hat = v_norm_hat / sqrt(sum(v_norm_hat**2))
      cell_height(i) = abs(sum(v_norm_hat*cell(:,i)))
   end do
   n_images = ceiling(cutoff/minval(cell_height))

   call matrix3x3_inverse(cell, cell_inv)
   ! into lattice coodinates 
   do i=1, N
       pos_l(1,i) = sum(cell_inv(1,:)*pos(:,i))
       pos_l(2,i) = sum(cell_inv(2,:)*pos(:,i))
       pos_l(3,i) = sum(cell_inv(3,:)*pos(:,i))
   end do
   d_pos_l(1) = sum(cell_inv(1,:)*d_pos(:))
   d_pos_l(2) = sum(cell_inv(2,:)*d_pos(:))
   d_pos_l(3) = sum(cell_inv(3,:)*d_pos(:))


   if (n_extra_data == 1 .and. allocated(new_extra_data)) then
      if (any(shape(new_extra_data) /= shape(extra_data))) then
	 deallocate(new_extra_data)
      endif
   endif
   if (n_extra_data == 1 .and. .not. allocated(new_extra_data)) then
      allocate(new_extra_data(n_extra_data, N))
   endif

   if (n_extra_data == 1) new_extra_data = extra_data

   dE = 0.0
   i=d_i
   do j=1,N
      if (j == i) cycle

      dr_l0 = pos_l(:,i) - pos_l(:,j)
      dr_l0 = dr_l0 - floor(dr_l0+0.5)

      drp_l0 = pos_l(:,i)+d_pos_l(:) - pos_l(:,j)
      drp_l0 = drp_l0 - floor(drp_l0+0.5)

      do dj1=-n_images,n_images
      dr_l(1) = dr_l0(1) + real(dj1, 8)
      drp_l(1) = drp_l0(1) + real(dj1, 8)
      do dj2=-n_images,n_images
      dr_l(2) = dr_l0(2) + real(dj2, 8)
      drp_l(2) = drp_l0(2) + real(dj2, 8)
      do dj3=-n_images,n_images
      dr_l(3) = dr_l0(3) + real(dj3, 8)
      drp_l(3) = drp_l0(3) + real(dj3, 8)

         dr(1) = sum(cell(1,:)*dr_l)
         dr(2) = sum(cell(2,:)*dr_l)
         dr(3) = sum(cell(3,:)*dr_l)
         drp(1) = sum(cell(1,:)*drp_l)
         drp(2) = sum(cell(2,:)*drp_l)
         drp(3) = sum(cell(3,:)*drp_l)
	 dr_mag_sq = sum(dr*dr)
	 drp_mag_sq = sum(drp*drp)

	 if (dr_mag_sq < cutoff_sq) then
	    dr_mag = sqrt(dr_mag_sq)
	    dE = dE -  epsilon*(((sigma/dr_mag)**12 - (sigma/dr_mag)**6) - E_offset)
	    if (n_extra_data == 1) then
	       new_extra_data(1,i) = new_extra_data(1,i) - 0.5*((1.0/dr_mag**12 - 1.0/dr_mag**6) - E_offset)
	       new_extra_data(1,j) = new_extra_data(1,j) - 0.5*((1.0/dr_mag**12 - 1.0/dr_mag**6) - E_offset)
	    endif
	 endif
	 if (drp_mag_sq < cutoff_sq) then
	    drp_mag = sqrt(drp_mag_sq)
	    dE = dE + epsilon*(((sigma/drp_mag)**12 - (sigma/drp_mag)**6) - E_offset)
	    if (n_extra_data == 1) then
	       new_extra_data(1,i) = new_extra_data(1,i) + 0.5*epsilon*(((sigma/drp_mag)**12 - (sigma/drp_mag)**6) - E_offset)
	       new_extra_data(1,j) = new_extra_data(1,j) + 0.5*epsilon*(((sigma/drp_mag)**12 - (sigma/drp_mag)**6) - E_offset)
	    endif
	 endif

      end do
      end do
      end do
   end do

   if (dE < dEmax) then ! accept
      pos(:,i) = pos(:,i) + d_pos(:)
      if (n_extra_data == 1) extra_data = new_extra_data
      ll_move_atom_1 = 1
   else ! reject
      dE = 0.0
      ll_move_atom_1 = 0
   endif

end function ll_move_atom_1

function ll_eval_forces(N, pos, n_extra_data, extra_data, cell, forces) result(energy)
use mat_mod
use LJ_params_mod
implicit none
   integer :: N
   double precision :: pos(3,N), cell(3,3), forces(3,N)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data, N)
   double precision :: energy ! result

   integer :: i, j
   double precision :: dr(3), dr_mag, dr_mag_sq, dr_l(3), dr_l0(3), pos_l(3,N)
   double precision :: cell_inv(3,3), E_term
   integer :: dj1, dj2, dj3

   integer n_images
   double precision cell_height(3), v_norm_hat(3)

   do i=1, 3
      v_norm_hat = cell(:,mod(i,3)+1) .cross. cell(:,mod(i+1,3)+1)
      v_norm_hat = v_norm_hat / sqrt(sum(v_norm_hat**2))
      cell_height(i) = abs(sum(v_norm_hat*cell(:,i)))
   end do
   n_images = ceiling(cutoff/minval(cell_height))

   call matrix3x3_inverse(cell, cell_inv)
   do i=1, N
       pos_l(1,i) = sum(cell_inv(1,:)*pos(:,i))
       pos_l(2,i) = sum(cell_inv(2,:)*pos(:,i))
       pos_l(3,i) = sum(cell_inv(3,:)*pos(:,i))
   end do

   if (n_extra_data == 1) extra_data = 0.0

   energy = 0.0
   forces = 0.0
   do i=1, N
   do j=i, N

      dr_l0 = pos_l(:,i) - pos_l(:,j)
      dr_l0 = dr_l0 - floor(dr_l0+0.5)
      do dj1=-n_images,n_images
      dr_l(1) = dr_l0(1) + real(dj1, 8)
      do dj2=-n_images,n_images
      dr_l(2) = dr_l0(2) + real(dj2, 8)
      do dj3=-n_images,n_images
      dr_l(3) = dr_l0(3) + real(dj3, 8)
      if (i == j .and. dj1 == 0 .and. dj2 == 0 .and. dj3 == 0) cycle

         dr(1) = sum(cell(1,:)*dr_l)
         dr(2) = sum(cell(2,:)*dr_l)
         dr(3) = sum(cell(3,:)*dr_l)
	 dr_mag_sq = sum(dr*dr)
	 if (dr_mag_sq < cutoff_sq) then
	    dr_mag = sqrt(dr_mag_sq)
	    E_term = epsilon*(((sigma/dr_mag)**12 - (sigma/dr_mag)**6) - E_offset)
	    if (i == j) E_term = E_term * 0.5
	    energy = energy + E_term
	    if (n_extra_data == 1) then
	       extra_data(1,i) = extra_data(1,i) + 0.5*E_term
	       extra_data(1,j) = extra_data(1,j) + 0.5*E_term
	    endif
	    if (i /= j) then
	       forces(:,i) = forces(:,i) - epsilon*(-12.0*sigma**12/dr_mag**13 + 6.0*sigma**6/dr_mag**7)*(dr/dr_mag)
	       forces(:,j) = forces(:,j) + epsilon*(-12.0*sigma**12/dr_mag**13 + 6.0*sigma**6/dr_mag**7)*(dr/dr_mag)
	    endif
	 endif

      end do
      end do
      end do
   end do
   end do

end function ll_eval_forces

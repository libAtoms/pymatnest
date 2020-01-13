subroutine ll_init_model(N_params, params)
use example_bead_spring_polymer_params_mod
implicit none
   integer :: N_params
   double precision :: params(N_params)

   ! PARAMS: chain_length K_bond r0_bond K_angle epsilon sigma cutoff
   if (N_params /= 7) then
      print *, "example_bead_spring_polymer_model.F90 got 7 /= N_params ", N_params
      print *, "Need: chain_length K_bond r0_bond K_angle epsilon sigma cutoff"
      call exit(1)
   endif 

   chain_length = int(params(1))
   K_bond = params(2)
   r0_bond = params(3)
   K_angle = params(4)
   epsilon = params(5)
   sigma = params(6)
   cutoff = params(7)

   cutoff_sq = cutoff*cutoff

   E_offset = (sigma/cutoff)**12 - (sigma/cutoff)**6

end subroutine ll_init_model

subroutine ll_init_config(N, Z, pos, cell, Emax)
implicit none
   integer :: N
   integer :: Z(N)
   double precision :: pos(3,N), cell(3,3)
   double precision :: Emax
   return
end subroutine ll_init_config

double precision function ll_eval_energy(N, Z, pos, n_extra_data, extra_data, cell)
use example_mat_mod
use example_bead_spring_polymer_params_mod
implicit none
   integer :: N
   integer :: Z(N)
   double precision :: pos(3,N), cell(3,3)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data, N)

   integer :: i, j
   double precision :: dr(3), dr_mag, dr_mag_sq, dr_l(3), dr_l0(3), pos_l(3,N), drr(3)

   double precision :: cell_inv(3,3), E_term
   integer :: dj1, dj2, dj3

   integer :: n_images
   double precision cell_height(3), v_norm_hat(3)

   integer :: i_chain, j_chain, n_chains, i_bond

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
   i_chain = int((i-1)/chain_length)
   do j=i, N
      j_chain = int((j-1)/chain_length)
      dr_l0 = pos_l(:,i)-pos_l(:,j)
      dr_l0 = dr_l0 - floor(dr_l0+0.5)
      do dj1=-n_images,n_images
      dr_l(1) = dr_l0(1) + real(dj1, 8)
      do dj2=-n_images,n_images
      dr_l(2) = dr_l0(2) + real(dj2, 8)
      do dj3=-n_images,n_images
      dr_l(3) = dr_l0(3) + real(dj3, 8)
	 if (i == j .and. dj1 == 0 .and. dj2 == 0 .and. dj3 == 0) cycle !  no self-interaction
         if (i_chain == j_chain .and. abs(i-j) == 1 .and. dj1 == 0 .and. dj2 == 0 .and. dj3 == 0) cycle ! bond excluded

	 dr(1) = cell(1,1) * dr_l(1) + cell(1,2) * dr_l(2) + cell(1,3) * dr_l(3) ! sum(cell(1,:)*dr_l)
	 dr(2) = cell(2,1) * dr_l(1) + cell(2,2) * dr_l(2) + cell(2,3) * dr_l(3) ! sum(cell(2,:)*dr_l)
	 dr(3) = cell(3,1) * dr_l(1) + cell(3,2) * dr_l(2) + cell(3,3) * dr_l(3) ! sum(cell(3,:)*dr_l)
	 dr_mag_sq = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3) ! sum(dr*dr)
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

   n_chains = N/chain_length
   do i_chain=1, n_chains
      ! bonds
      do i_bond=(i_chain-1)*chain_length+1, i_chain*chain_length-1
         dr = pos(:,i_bond+1) - pos(:, i_bond)
         dr_mag = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
         ll_eval_energy = ll_eval_energy + K_bond * (dr_mag-r0_bond)*(dr_mag-r0_bond)
      end do
      ! angles
      do i_bond=(i_chain-1)*chain_length+2, i_chain*chain_length-1
         dr = pos(:,i_bond+1) - pos(:, i_bond)
         dr = dr / sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))

         drr = pos(:,i_bond-1) - pos(:, i_bond)
         drr = drr / sqrt(drr(1)*drr(1) + drr(2)*drr(2) + drr(3)*drr(3))

         ll_eval_energy = ll_eval_energy + K_angle * (1.0 + dr(1)*drr(1) + dr(2)*drr(2) + dr(3)*drr(3))
      end do
   end do

end function ll_eval_energy

integer function ll_move_atom_1(N, Z, pos, n_extra_data, extra_data, cell, d_i, d_pos, dEmax, dE)
use example_mat_mod
use example_bead_spring_polymer_params_mod
implicit none
   integer :: N
   integer :: Z(N)
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

   print *, "example_bead_spring_polymer_model does not implement move_atom_1 for MC"
   call exit(1)

!    do i=1, 3
!       v_norm_hat = cell(:,mod(i,3)+1) .cross. cell(:,mod(i+1,3)+1)
!       v_norm_hat = v_norm_hat / sqrt(sum(v_norm_hat**2))
!       cell_height(i) = abs(sum(v_norm_hat*cell(:,i)))
!    end do
!    n_images = ceiling(cutoff/minval(cell_height))
! 
!    call matrix3x3_inverse(cell, cell_inv)
!    ! into lattice coodinates 
!    do i=1, N
!        pos_l(1,i) = sum(cell_inv(1,:)*pos(:,i))
!        pos_l(2,i) = sum(cell_inv(2,:)*pos(:,i))
!        pos_l(3,i) = sum(cell_inv(3,:)*pos(:,i))
!    end do
!    d_pos_l(1) = sum(cell_inv(1,:)*d_pos(:))
!    d_pos_l(2) = sum(cell_inv(2,:)*d_pos(:))
!    d_pos_l(3) = sum(cell_inv(3,:)*d_pos(:))
! 
! 
!    if (n_extra_data == 1 .and. allocated(new_extra_data)) then
!       if (any(shape(new_extra_data) /= shape(extra_data))) then
! 	 deallocate(new_extra_data)
!       endif
!    endif
!    if (n_extra_data == 1 .and. .not. allocated(new_extra_data)) then
!       allocate(new_extra_data(n_extra_data, N))
!    endif
! 
!    if (n_extra_data == 1) new_extra_data = extra_data
! 
!    dE = 0.0
!    i=d_i
!    do j=1,N
!       if (j == i) cycle
! 
!       dr_l0 = pos_l(:,i) - pos_l(:,j)
!       dr_l0 = dr_l0 - floor(dr_l0+0.5)
! 
!       drp_l0 = pos_l(:,i)+d_pos_l(:) - pos_l(:,j)
!       drp_l0 = drp_l0 - floor(drp_l0+0.5)
! 
!       do dj1=-n_images,n_images
!       dr_l(1) = dr_l0(1) + real(dj1, 8)
!       drp_l(1) = drp_l0(1) + real(dj1, 8)
!       do dj2=-n_images,n_images
!       dr_l(2) = dr_l0(2) + real(dj2, 8)
!       drp_l(2) = drp_l0(2) + real(dj2, 8)
!       do dj3=-n_images,n_images
!       dr_l(3) = dr_l0(3) + real(dj3, 8)
!       drp_l(3) = drp_l0(3) + real(dj3, 8)
! 
! 	 dr(1) = cell(1,1) * dr_l(1) + cell(1,2) * dr_l(2) + cell(1,3) + dr_l(3) ! sum(cell(1,:)*dr_l)
! 	 dr(2) = cell(2,1) * dr_l(1) + cell(2,2) * dr_l(2) + cell(2,3) + dr_l(3) ! sum(cell(2,:)*dr_l)
! 	 dr(3) = cell(3,1) * dr_l(1) + cell(3,2) * dr_l(2) + cell(3,3) + dr_l(3) ! sum(cell(3,:)*dr_l)
! 
! 	 drp(1) = cell(1,1) * drp_l(1) + cell(1,2) * drp_l(2) + cell(1,3) + drp_l(3) ! sum(cell(1,:)*drp_l)
! 	 drp(2) = cell(2,1) * drp_l(1) + cell(2,2) * drp_l(2) + cell(2,3) + drp_l(3) ! sum(cell(2,:)*drp_l)
! 	 drp(3) = cell(3,1) * drp_l(1) + cell(3,2) * drp_l(2) + cell(3,3) + drp_l(3) ! sum(cell(3,:)*drp_l)
! 
! 	 dr_mag_sq = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3) ! sum(dr*dr)
! 	 drp_mag_sq = drp(1)*drp(1) + drp(2)*drp(2) + drp(3)*drp(3) ! sum(drp*drp)
! 
! 	 if (dr_mag_sq < cutoff_sq) then
! 	    dr_mag = sqrt(dr_mag_sq)
! 	    dE = dE -  epsilon*(((sigma/dr_mag)**12 - (sigma/dr_mag)**6) - E_offset)
! 	    if (n_extra_data == 1) then
! 	       new_extra_data(1,i) = new_extra_data(1,i) - 0.5*((sigma/dr_mag)**12 - &
!                                      (sigma/dr_mag)**6 - E_offset)
! 	       new_extra_data(1,j) = new_extra_data(1,j) - 0.5*((sigma/dr_mag)**12 - &
!                                      (sigma/dr_mag)**6 - E_offset)
! 	    endif
! 	 endif
! 	 if (drp_mag_sq < cutoff_sq) then
! 	    drp_mag = sqrt(drp_mag_sq)
! 	    dE = dE + epsilon*(((sigma/drp_mag)**12 - (sigma/drp_mag)**6) - E_offset)
! 	    if (n_extra_data == 1) then
! 	       new_extra_data(1,i) = new_extra_data(1,i) + 0.5*epsilon*(((sigma/drp_mag)**12 - &
!                                      (sigma/drp_mag)**6) - E_offset)
! 	       new_extra_data(1,j) = new_extra_data(1,j) + 0.5*epsilon*(((sigma/drp_mag)**12 - &
!                                      (sigma/drp_mag)**6) - E_offset)
! 	    endif
! 	 endif
! 
!       end do
!       end do
!       end do
!    end do
! 
!    if (dE < dEmax) then ! accept
!       pos(:,i) = pos(:,i) + d_pos(:)
!       if (n_extra_data == 1) extra_data = new_extra_data
!       ll_move_atom_1 = 1
!    else ! reject
!       dE = 0.0
!       ll_move_atom_1 = 0
!    endif

end function ll_move_atom_1

function ll_eval_forces(N, Z, pos, n_extra_data, extra_data, cell, forces) result(energy)
use example_mat_mod
use example_bead_spring_polymer_params_mod
implicit none
   integer :: N
   integer :: Z(N)
   double precision :: pos(3,N), cell(3,3), forces(3,N)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data, N)
   double precision :: energy ! result

   integer :: i, j
   double precision :: dr(3), dr_mag, dr_mag_sq, dr_l(3), dr_l0(3), pos_l(3,N), drr(3), drr_mag
   double precision :: cell_inv(3,3), E_term
   integer :: dj1, dj2, dj3

   integer n_images
   double precision cell_height(3), v_norm_hat(3), f_term(3)

   integer :: i_chain, j_chain, n_chains, i_bond

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
   i_chain = int((i-1)/chain_length)
   do j=i, N
      j_chain = int((j-1)/chain_length)

      dr_l0 = pos_l(:,i) - pos_l(:,j)
      dr_l0 = dr_l0 - floor(dr_l0+0.5)
      do dj1=-n_images,n_images
      dr_l(1) = dr_l0(1) + real(dj1, 8)
      do dj2=-n_images,n_images
      dr_l(2) = dr_l0(2) + real(dj2, 8)
      do dj3=-n_images,n_images
      dr_l(3) = dr_l0(3) + real(dj3, 8)
         if (i == j .and. dj1 == 0 .and. dj2 == 0 .and. dj3 == 0) cycle ! no self-interaction
         if (i_chain == j_chain .and. abs(i-j) == 1 .and. dj1 == 0 .and. dj2 == 0 .and. dj3 == 0) cycle ! bond excluded

	 dr(1) = cell(1,1) * dr_l(1) + cell(1,2) * dr_l(2) + cell(1,3) * dr_l(3) ! sum(cell(1,:)*dr_l)
	 dr(2) = cell(2,1) * dr_l(1) + cell(2,2) * dr_l(2) + cell(2,3) * dr_l(3) ! sum(cell(2,:)*dr_l)
	 dr(3) = cell(3,1) * dr_l(1) + cell(3,2) * dr_l(2) + cell(3,3) * dr_l(3) ! sum(cell(3,:)*dr_l)
	 dr_mag_sq = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3) ! sum(dr*dr)
	 if (dr_mag_sq < cutoff_sq) then
	    dr_mag = sqrt(dr_mag_sq)
	    E_term = epsilon*((sigma/dr_mag)**12 - (sigma/dr_mag)**6 - E_offset)
	    if (i == j) E_term = E_term * 0.5
	    energy = energy + E_term
	    if (n_extra_data == 1) then
	       extra_data(1,i) = extra_data(1,i) + 0.5*E_term
	       extra_data(1,j) = extra_data(1,j) + 0.5*E_term
	    endif
	    if (i /= j) then
               f_term = epsilon*(-12.0*sigma**12/dr_mag**13 + 6.0*sigma**6/dr_mag**7)*(dr/dr_mag)
	       forces(:,i) = forces(:,i) - f_term
	       forces(:,j) = forces(:,j) + f_term
	    endif
	 endif

      end do
      end do
      end do
   end do
   end do

   n_chains = N/chain_length
   do i_chain=1, n_chains
      ! bonds
      do i_bond=(i_chain-1)*chain_length+1, i_chain*chain_length-1
         dr = pos(:,i_bond+1) - pos(:, i_bond)
         dr_mag_sq = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
         dr_mag = sqrt(dr_mag_sq)
         energy = energy + K_bond * (dr_mag-r0_bond)*(dr_mag-r0_bond)
         f_term = K_bond*2.0*(dr_mag-r0_bond)*dr/dr_mag 
         forces(:,i_bond)   = forces(:,i_bond)   + f_term
         forces(:,i_bond+1) = forces(:,i_bond+1) - f_term
      end do
      ! angles
      do i_bond=(i_chain-1)*chain_length+2, i_chain*chain_length-1
         dr = pos(:,i_bond+1) - pos(:, i_bond)
         dr_mag = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
         dr = dr / dr_mag

         drr = pos(:,i_bond-1) - pos(:, i_bond)
         drr_mag = sqrt(drr(1)*drr(1) + drr(2)*drr(2) + drr(3)*drr(3))
         drr = drr / drr_mag

         energy = energy + K_angle * (1.0 + dr(1)*drr(1) + dr(2)*drr(2) + dr(3)*drr(3))
         ! dr
         f_term = K_angle * (drr/dr_mag - (dr(1)*drr(1) + dr(2)*drr(2) + dr(3)*drr(3)) * dr/dr_mag)
         forces(:,i_bond)   = forces(:, i_bond)   + f_term
         forces(:,i_bond+1) = forces(:, i_bond+1) - f_term
         ! drr
         f_term = K_angle * (dr/drr_mag - (dr(1)*drr(1) + dr(2)*drr(2) + dr(3)*drr(3)) * drr/drr_mag)
         forces(:,i_bond)   = forces(:, i_bond)   + f_term
         forces(:,i_bond-1) = forces(:, i_bond-1) - f_term
      end do
   end do

end function ll_eval_forces

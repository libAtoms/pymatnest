function fortran_seed_size()
   integer :: fortran_seed_size

   call random_seed(size=fortran_seed_size)
end function fortran_seed_size

subroutine fortran_set_seed(n_seed, seed)
   integer :: seed(n_seed)

   call random_seed(put=seed)
end subroutine fortran_set_seed


subroutine fortran_MC_atom_velo(N, vel, mass, n_steps, step_size, KEmax, final_KE, n_accept)
   implicit none
   integer :: N
   double precision :: vel(3,N), mass(N)
   integer :: n_steps
   double precision :: step_size, KEmax, final_KE
   integer :: n_accept

   integer :: d_i
   double precision :: d_r, KE, dKE, d_vel(3)


   integer :: i_step, i_at, t_i
   integer :: order(N)

   n_accept = 0
   KE = 0.5*sum(spread(mass,1,3)*vel**2)

   do i_step=1, n_steps

      do i_at=1, N
	 order(i_at) = i_at
      end do
      do i_at=1, N-1
	 call random_number(d_r); d_i = floor(d_r*(N-i_at+1))+i_at
	 if (d_i /= i_at) then
	    t_i = order(i_at)
	    order(i_at) = order(d_i)
	    order(d_i) = t_i
	 endif
      end do

      do i_at=1, N
	 d_i = order(i_at)
	 call random_number(d_vel)
	 d_vel = 2.0*step_size*(d_vel-0.5)

	 dKE = 0.5*mass(d_i)*(sum((vel(:,d_i)+d_vel(:))**2) - sum(vel(:,d_i)**2))
	 if (KE + dKE < KEmax) then
	    vel(1:3,d_i) = vel(1:3,d_i) + d_vel(1:3)
	    KE = KE + dKE
	    n_accept = n_accept + 1
	 endif
      end do
   end do

   final_KE = KE

end subroutine fortran_MC_atom_velo

subroutine fortran_MC_atom(N, pos, vel, mass, cell, n_steps, step_size_pos, step_size_vel, Emax, final_E, &
			   n_accept_pos, n_accept_vel)
   implicit none
   integer :: N
   double precision :: pos(3,N), vel(3,N), mass(N), cell(3,3)
   integer :: n_steps
   double precision :: step_size_pos, step_size_vel, Emax, final_E
   integer :: n_accept_pos, n_accept_vel

   logical :: do_vel
   integer :: d_i
   double precision :: d_r, E, dE, d_pos(3), d_vel(3)

   double precision, external :: ll_eval_energy, ll_eval_denergy_1

   integer :: i_step, i_at, t_i
   integer :: order(N)
   double precision :: vel_pos_rv

   do_vel = (step_size_vel /= 0.0)

   n_accept_pos = 0
   n_accept_vel = 0
   E = ll_eval_energy(N, pos, cell)
   if (do_vel) then
      E = E + 0.5*sum(spread(mass,1,3)*vel**2)
   endif

   do i_step=1, n_steps

      do i_at=1, N
	 order(i_at) = i_at
      end do
      do i_at=1, N-1
	 call random_number(d_r); d_i = floor(d_r*(N-i_at+1))+i_at
	 if (d_i /= i_at) then
	    t_i = order(i_at)
	    order(i_at) = order(d_i)
	    order(d_i) = t_i
	 endif
      end do

      do i_at=1, N
	 d_i = order(i_at)
	 if (do_vel) then
	    call random_number(d_vel)
	    d_vel = 2.0*step_size_vel*(d_vel-0.5)
	    call random_number(vel_pos_rv)
	 endif

	 if (do_vel .and.  vel_pos_rv < 0.5) then
	    dE = 0.5*mass(d_i)*(sum((vel(:,d_i)+d_vel(:))**2) - sum(vel(:,d_i)**2))
	    if (E + dE < Emax) then
	       vel(1:3,d_i) = vel(1:3,d_i) + d_vel(1:3)
	       E = E + dE
	       n_accept_vel = n_accept_vel + 1
	    endif
	 endif

	 call random_number(d_pos)
	 d_pos = 2.0*step_size_pos*(d_pos-0.5)
	 dE = ll_eval_denergy_1(N, pos, cell, d_i, d_pos)
	 if (E + dE < Emax) then
	    pos(1:3,d_i) = pos(1:3,d_i) + d_pos(1:3)
	    E = E + dE
	    n_accept_pos = n_accept_pos + 1
	 endif

	 if (do_vel .and.  vel_pos_rv >= 0.5) then
	    dE = 0.5*mass(d_i)*(sum((vel(:,d_i)+d_vel(:))**2) - sum(vel(:,d_i)**2))
	    if (E + dE < Emax) then
	       vel(1:3,d_i) = vel(1:3,d_i) + d_vel(1:3)
	       E = E + dE
	       n_accept_vel = n_accept_vel + 1
	    endif
	 endif

      end do
   end do

   final_E = E

end subroutine fortran_MC_atom

subroutine fortran_MD_atom_NVE(N, pos, vel, mass, cell, n_steps, timestep, final_E, debug)
   implicit none
   integer :: N
   double precision :: pos(3,N), vel(3,N), mass(N), cell(3,3)
   integer :: n_steps
   double precision :: timestep, final_E
   integer :: debug

   integer i_step, i
   double precision :: forces(3,N), acc(3,N), PE

   double precision, external :: ll_eval_forces, ll_eval_energy

   ! initialize accelerations
   PE = ll_eval_forces(N, pos, cell, forces)
   do i=1, 3
      acc(i,:) = forces(i,:) / mass(:)
   end do

   if (debug > 0) print *, "initial PE KE E ", PE, 0.5*sum(spread(mass,1,3)*vel**2), &
      PE+0.5*sum(spread(mass,1,3)*vel**2)

   do i_step=1, n_steps
      ! Verlet part 1
      vel = vel + 0.5*timestep*acc
      pos = pos + timestep*vel

      ! new accelerations at t+dt
      PE = ll_eval_forces(N, pos, cell, forces)
      do i=1, 3
	 acc(i,:) = forces(i,:) / mass(:)
      end do

      ! Verlet part 2
      vel = vel + 0.5*timestep*acc

      if (debug > 0) print *, "step PE KE E ", i_step, PE, 0.5*sum(spread(mass,1,3)*vel**2), &
	 PE+0.5*sum(spread(mass,1,3)*vel**2)

   end do

   final_E = PE + 0.5*sum(spread(mass,1,3)*vel**2)

end subroutine fortran_MD_atom_NVE

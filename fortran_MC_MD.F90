function fortran_seed_size()
   integer :: fortran_seed_size

   call random_seed(size=fortran_seed_size)
end function fortran_seed_size

subroutine fortran_set_seed(n_seed, seed)
   integer :: seed(n_seed)

   call random_seed(put=seed)
end subroutine fortran_set_seed


subroutine fortran_MC_atom_velo(N, vel, mass, n_steps, step_size, KEmax, final_KE, n_try, n_accept)
   implicit none
   integer :: N
   double precision :: vel(3,N), mass(N)
   integer :: n_steps
   double precision :: step_size, KEmax, final_KE
   integer :: n_try, n_accept

   integer :: d_i
   double precision :: d_r, KE, dKE, d_vel(3)


   integer :: i_step, i_at, t_i
   integer :: order(N)

   n_try = N*n_steps
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

subroutine fortran_MC_atom(N, Z, pos, vel, mass, n_extra_data, extra_data, cell, n_steps, &
                           step_size_pos, step_size_vel, Emax, KEmax, final_E, n_try, n_accept_pos, n_accept_vel)
   implicit none
   integer :: N
   integer :: Z(N)
   double precision :: pos(3,N), vel(3,N), mass(N), cell(3,3)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data,N)
   integer :: n_steps
   double precision :: step_size_pos, step_size_vel, Emax, KEmax, final_E
   integer :: n_try, n_accept_pos, n_accept_vel

   logical :: do_vel
   integer :: d_i
   double precision :: d_r, E, dE, KE, dKE, d_pos(3), d_vel(3)

   double precision, external :: ll_eval_energy
   integer, external :: ll_move_atom_1

   integer :: i_step, i_at, t_i
   integer :: order(N)
   double precision :: vel_pos_rv

   do_vel = (step_size_vel /= 0.0)

   n_try = N*n_steps
   n_accept_pos = 0
   n_accept_vel = 0
   E = ll_eval_energy(N, Z, pos, n_extra_data, extra_data, cell)
   if (do_vel) then
      KE = 0.5*sum(spread(mass,1,3)*vel**2)
      E = E + KE
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
            dKE = 0.5*mass(d_i)*(sum((vel(:,d_i)+d_vel(:))**2) - sum(vel(:,d_i)**2))
            dE = dKE
            if (E + dE < Emax .and. (KEmax <= 0.0 .or. KE + dKE < KEmax)) then
               vel(1:3,d_i) = vel(1:3,d_i) + d_vel(1:3)
               E = E + dE
               KE = KE + dKE
               n_accept_vel = n_accept_vel + 1
            endif
         endif

         call random_number(d_pos)
         d_pos = 2.0*step_size_pos*(d_pos-0.5)
         n_accept_pos = n_accept_pos + ll_move_atom_1(N, Z, pos, n_extra_data, extra_data, cell, d_i, d_pos, Emax-E, dE)
         E = E + dE

         if (do_vel .and.  vel_pos_rv >= 0.5) then
            dKE = 0.5*mass(d_i)*(sum((vel(:,d_i)+d_vel(:))**2) - sum(vel(:,d_i)**2))
            dE = dKE
            if (E + dE < Emax .and. (KEmax <= 0.0 .or. KE + dKE < KEmax)) then
               vel(1:3,d_i) = vel(1:3,d_i) + d_vel(1:3)
               E = E + dE
               KE = KE + dKE
               n_accept_vel = n_accept_vel + 1
            endif
         endif

      end do
   end do

   final_E = E

end subroutine fortran_MC_atom

subroutine fortran_GMC_atom(N, Z, pos, mass, n_extra_data, extra_data, cell, n_steps, &
                          step_size_pos, Emax, final_E, n_try, n_accept, d_pos_hat, debug)
   implicit none
   integer :: N
   integer :: Z(N)
   double precision :: pos(3,N), mass(N), cell(3,3)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data,N)
   integer :: n_steps
   double precision :: step_size_pos, Emax, final_E
   integer :: n_try, n_accept
   double precision :: d_pos_hat(3,N)
   integer :: debug

   double precision :: E, d_pos(3,N), cur_pos(3,N), Fhat(3,N), d_pos_reflect(3,N), &
                       last_good_pos(3,N), last_good_d_pos(3,N)

   double precision, external :: ll_eval_energy, ll_eval_forces

   integer :: n_reflect, n_reverse
   integer :: i_step, i
   logical :: reflected_last, reverse, reflect, step_passed

   d_pos = d_pos_hat * step_size_pos

   n_reflect = 0
   n_reverse = 0

   do i_step=1, n_steps

      last_good_pos = pos
      last_good_d_pos = d_pos
      reflected_last = .false. 
      step_passed = .true.

      pos = pos + d_pos

      reflect_reverse: do i = 1, 2 ! reflect and then reverse if needed
        E = ll_eval_energy(N, Z, pos, n_extra_data, extra_data, cell)
        if (isnan(E) .or. E >= Emax) then ! reflect or reverse
            step_passed = .false.
          
            if (.not. reflected_last) then
              reflect = .true.
              reverse = .false.
            else
              reverse = .true.
              reflect = .false.
            end if
     
            if (reflect) then
              E = ll_eval_forces(N, Z, pos, n_extra_data, extra_data, cell, Fhat)
              Fhat = Fhat / sqrt(sum(Fhat*Fhat))
     
              d_pos_reflect = d_pos - 2.0*Fhat*sum(Fhat*d_pos)
              pos = last_good_pos + d_pos + d_pos_reflect
              d_pos = d_pos_reflect
     
              reflected_last = .true.
              n_reflect = n_reflect + 1
            else if (reverse) then
              pos = last_good_pos
              d_pos = - last_good_d_pos
              ! undo reflection that just happened, replace it with a reversal
              n_reflect = n_reflect - 1 
              n_reverse = n_reverse + 1
            endif
   
        else ! coordinates (possibly after reflection) are OK. Save them.
          exit reflect_reverse
        end if
      end do reflect_reverse

   end do

   if (n_reverse > 0) E = ll_eval_energy(N, Z, pos, n_extra_data, extra_data, cell)

   ! set n_accept_pos so that acceptance rate is what we need (see above)
   n_try = n_reflect + n_reverse
   n_accept = n_reflect
   final_E = E

end subroutine fortran_GMC_atom

subroutine fortran_GMC_noreverse_atom(N, Z, pos, mass, n_extra_data, extra_data, cell, n_steps, &
                           step_size_pos, Emax, final_E, n_try, n_accept, d_pos_hat, debug)
   implicit none
   integer :: N
   integer :: Z(N)
   double precision :: pos(3,N), mass(N), cell(3,3)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data,N)
   integer :: n_steps
   double precision :: step_size_pos, Emax, final_E
   integer :: n_try, n_accept
   double precision :: d_pos_hat(3,N)
   integer :: debug

   double precision :: E, d_pos(3,N), cur_pos(3,N), Fhat(3,N)

   double precision, external :: ll_eval_energy, ll_eval_forces

   integer :: i_step

   E = ll_eval_energy(N, Z, pos, n_extra_data, extra_data, cell)
   final_E = E ! in case we reject
   d_pos = d_pos_hat * step_size_pos

   n_try = n_steps*N

   ! move a copy of pos, in case we reject
   cur_pos = pos
   do i_step=1, n_steps
      cur_pos = cur_pos + d_pos
      E = ll_eval_energy(N, Z, cur_pos, n_extra_data, extra_data, cell)
      if (isnan(E)) then ! can't reflect from NaN, forces aren't well defined
          exit
      endif
      if (E >= Emax) then ! reflect
          E = ll_eval_forces(N, Z, cur_pos, n_extra_data, extra_data, cell, Fhat)
          Fhat = Fhat / sqrt(sum(Fhat*Fhat))
          d_pos = d_pos - 2.0*Fhat*sum(Fhat*d_pos)
      end if
   end do

   if (E < Emax) then ! NaN should be false also
     final_E = E
     pos = cur_pos
     n_accept = n_steps*N
   else
     n_accept = 0
   endif

end subroutine fortran_GMC_noreverse_atom

subroutine fortran_MD_atom_NVE(N, Z, pos, vel, mass, n_extra_data, extra_data, cell, n_steps, timestep, final_E, debug)
   implicit none
   integer :: N
   integer :: Z(N)
   double precision :: pos(3,N), vel(3,N), mass(N), cell(3,3)
   integer :: n_extra_data
   double precision :: extra_data(n_extra_data,N)
   integer :: n_steps
   double precision :: timestep, final_E
   integer :: debug

   integer i_step, i
   double precision :: forces(3,N), acc(3,N), PE

   double precision, external :: ll_eval_forces, ll_eval_energy

   ! initialize accelerations
   PE = ll_eval_forces(N, Z, pos, n_extra_data, extra_data, cell, forces)
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
      PE = ll_eval_forces(N, Z, pos, n_extra_data, extra_data, cell, forces)
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

function fortran_MC_atom(N, pos, cell, n_steps, step_size, Emax, final_E) result(n_accept)
   implicit none
   integer :: N
   double precision :: pos(3,N), cell(3,3)
   integer :: n_steps
   double precision :: step_size, Emax, final_E
   ! result
   integer :: n_accept

   integer :: d_i
   double precision :: d_r, E, dE, d_pos(3)

   double precision, external :: ll_eval_energy, ll_eval_denergy_1

   integer :: i_step, i_at, t_i
   integer :: order(N)

   do i_at=1, N
      order(i_at) = i_at
   end do
   do i_at=1, N-1
      call random_number(d_r); d_i = floor(d_r*N)+1
      if (d_i /= i_at) then
	 t_i = order(i_at)
	 order(i_at) = order(d_i)
	 order(d_i) = t_i
      endif
   end do

   n_accept = 0
   E = ll_eval_energy(N, pos, cell)
   do i_step=1, n_steps
   do i_at=1, N
      d_i = order(i_at)
      call random_number(d_pos)
      d_pos = 2.0*step_size*(d_pos-0.5)
      dE = ll_eval_denergy_1(N, pos, cell, d_i, d_pos)
      if (E + dE < Emax) then
	 pos(1:3,d_i) = pos(1:3,d_i) + d_pos(1:3)
	 E = E + dE
	 n_accept = n_accept + 1
      endif
      if (E > Emax) print *, "fortran_MC_atom detected E > Emax ", E, Emax !DEBUG
   end do
   end do

   final_E = E

end function fortran_MC_atom

subroutine fortran_MD_atom_NVE(N, pos, vel, mass, cell, n_steps, timestep, final_E)
   implicit none
   integer :: N
   double precision :: pos(3,N), vel(3,N), mass(N), cell(3,3)
   integer :: n_steps
   double precision :: timestep, final_E

   integer i_step, i
   double precision :: forces(3,N), acc(3,N)

   double precision, external :: ll_eval_forces, ll_eval_energy

   ! initialize accelerations
   final_E = ll_eval_forces(N, pos, cell, forces)
   do i=1, 3
      acc(i,:) = forces(i,:) / mass(:)
   end do

   do i_step=1, n_steps
      ! Verlet part 1
      vel = vel + 0.5*timestep*acc
      pos = pos + timestep*vel

      ! new forces
      final_E = ll_eval_forces(N, pos, cell, forces)
      do i=1, 3
	 acc(i,:) = forces(i,:) / mass(:)
      end do

      ! Verlet part 2
      vel = vel + 0.5*timestep*acc
   end do

   final_E = ll_eval_energy(N, pos, cell)
end subroutine fortran_MD_atom_NVE

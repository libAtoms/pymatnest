function fortran_MC(N, pos, cell, n_steps, step_size, Emax, final_E) result(n_accept)
   implicit none
   integer :: N
   double precision :: pos(3,N), cell(3,3)
   integer :: n_steps
   double precision :: step_size, Emax, final_E
   ! result
   integer :: n_accept

   integer :: d_i
   double precision :: d_r, E, dE, d_pos(3)

   double precision, external :: ll_eval_energy, ll_eval_energy_1

   integer i_step, i_at

   n_accept = 0
   E = ll_eval_energy(N, pos, cell)
   do i_step=1, n_steps
   do i_at=1, N
      call random_number(d_r)
      d_i = floor(d_r*N)+1
      call random_number(d_pos)
      d_pos = 2.0*step_size*(d_pos-0.5)
      dE = ll_eval_energy_1(N, pos, cell, d_i, d_pos)
      if (E + dE < Emax) then
	 pos(1:3,d_i) = pos(1:3,d_i) + d_pos(1:3)
	 E = E + dE
	 n_accept = n_accept + 1
      endif
      if (E > Emax) print *, "fortran_MC detected E > Emax ", E, Emax !DEBUG
   end do
   end do

   final_E = E

end function fortran_MC

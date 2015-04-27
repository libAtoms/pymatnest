double precision function ll_eval_energy(N, pos, cell)
implicit none
   integer :: N
   double precision :: pos(3,N), cell(3,3)

   integer :: i, j
   double precision :: dr
   double precision :: E_offset  = 1.0/3.0**12 - 1.0/3.0**6

   ll_eval_energy = 0.0
   do i=1, N
   do j=i+1, N
      dr = sqrt(sum((pos(:,i)-pos(:,j))**2))
      if (dr < 3.0) then
	 ll_eval_energy = ll_eval_energy + ((1.0/dr**12 - 1.0/dr**6) - E_offset)
      endif
   end do
   end do
end function ll_eval_energy

double precision function ll_eval_energy_1(N, pos, cell, d_i, d_pos)
implicit none
   integer :: N
   double precision :: pos(3,N), cell(3,3)
   integer :: d_i
   double precision :: d_pos(3)

   double precision, external :: ll_eval_energy

   double precision E0, E1, pos_save(3)

   E0 = ll_eval_energy(N, pos, cell)
   pos_save(:) = pos(:,d_i)
   pos(:,d_i) = pos(:,d_i) + d_pos(:)
   E1 = ll_eval_energy(N, pos, cell)
   pos(:,d_i) = pos_save(:)

   ll_eval_energy_1 = E1-E0

end function ll_eval_energy_1

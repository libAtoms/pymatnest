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
!
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

   ll_eval_energy = 0.0

end function ll_eval_energy

double precision function ll_eval_denergy_1(N, pos, cell, d_i, d_pos)
use mat_mod
implicit none
   integer :: N
   double precision :: pos(3,N), cell(3,3)
   integer :: d_i
   double precision :: d_pos(3)

   ll_eval_denergy_1 = 0.0

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

   energy = 0.0
   forces = 0.0

end function ll_eval_forces

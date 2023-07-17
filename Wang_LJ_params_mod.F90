module Wang_LJ_params_mod
   integer, parameter :: N_Z = 2
   double precision :: epsilon(N_Z,N_Z)
   double precision :: sigma(N_Z, N_Z)
   double precision :: cutoff(N_Z,N_Z), cutoff_sq(N_Z,N_Z)
   double precision :: alpha(N_Z,N_Z)
   integer :: mu, nu

end module Wang_LJ_params_mod

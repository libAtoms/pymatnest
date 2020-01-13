module example_LJ_params_mod
   integer, parameter :: N_Z = 2
   double precision :: epsilon(N_Z,N_Z)
   double precision :: sigma(N_Z, N_Z)
   double precision :: cutoff(N_Z,N_Z), cutoff_sq(N_Z,N_Z)
   double precision :: E_offset(N_Z, N_Z)
end module example_LJ_params_mod

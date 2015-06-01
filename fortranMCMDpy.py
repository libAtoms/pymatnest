import ctypes
from numpy.ctypeslib import ndpointer
import numpy as np
import os, ase

class fortran_MC_MD:
   def __init__(self, model_so):
      if model_so.find("/") == 0:
	 self.model_lib = ctypes.CDLL(model_so, mode=ctypes.RTLD_GLOBAL)
      else:
	 if 'PYMATNEST_PATH' in os.environ:
	    self.model_lib = ctypes.CDLL(os.environ['PYMATNEST_PATH']+"/"+model_so, mode=ctypes.RTLD_GLOBAL)
	 else:
	    self.model_lib = ctypes.CDLL(os.path.dirname(__file__)+"/"+model_so, mode=ctypes.RTLD_GLOBAL)


      # ll_init_model

      # ll_init_config
      self.model_lib.ll_init_config_.argtypes = [ctypes.c_void_p, # N
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # pos
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # cell
	 ctypes.c_void_p ] # Emax

      # ll_eval_energy
      self.model_lib.ll_eval_energy_.restype = ctypes.c_double
      self.model_lib.ll_eval_energy_.argtypes = [ctypes.c_void_p, # N
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # pos
	 ctypes.c_void_p, # n_extra_data
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # extra_data
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")] # cell

      # ll_move_atom_1
      self.model_lib.ll_move_atom_1_.restype = ctypes.c_int
      self.model_lib.ll_move_atom_1_.argtypes = [ctypes.c_void_p, # N
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # pos
	 ctypes.c_void_p, # n_extra_data
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # extra_data
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # cell
	 ctypes.c_void_p, # d_i
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # d_pos
	 ctypes.c_void_p, # dE_max
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS") ] # dE

      # ll_eval_forces
      self.model_lib.ll_eval_forces_.restype = ctypes.c_double
      self.model_lib.ll_eval_forces_.argtypes = [ctypes.c_void_p, # N
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # pos
	 ctypes.c_void_p, # n_extra_data
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # extra_data
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # cell
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")] # forces

      if 'PYMATNEST_PATH' in os.environ:
	 self.lib = ctypes.CDLL(os.environ['PYMATNEST_PATH']+"/fortran_MC_MD.so")
      else:
	 self.lib = ctypes.CDLL(os.path.dirname(__file__)+"/fortran_MC_MD.so")

      # fortran_seed_size
      self.lib.fortran_seed_size_.restype = ctypes.c_int

      # fortran_set_seed
      self.lib.fortran_set_seed_.argtypes = [ctypes.c_void_p, # n_seed
	 ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")] # seed

      # fortran_MC_atom_velo
      self.lib.fortran_mc_atom_velo_.argtypes = [ctypes.c_void_p, # N
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # velo
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # masses
	 ctypes.c_void_p, # n_steps
	 ctypes.c_void_p, # step_size
	 ctypes.c_void_p, # KEmax
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # final_KE
	 ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")] # n_accept

      # fortran_MC_atom
      self.lib.fortran_mc_atom_.argtypes = [ctypes.c_void_p, # N
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # pos
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # velo
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # masses
	 ctypes.c_void_p, # n_extra_data
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # extra_data
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # cell
	 ctypes.c_void_p, # n_steps
	 ctypes.c_void_p, # step_size_pos
	 ctypes.c_void_p, # step_size_velo
	 ctypes.c_void_p, # Emax
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # final_E
	 ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"), # n_accept
	 ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")] # n_accept_velo

      # fortran_MD_atom_NVE
      self.lib.fortran_md_atom_nve_.argtypes = [ctypes.c_void_p, # N
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # pos
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # vel
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # mass
	 ctypes.c_void_p, # n_extra_data
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # extra_data
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # cell
	 ctypes.c_void_p, # n_steps
	 ctypes.c_void_p, # timestep
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # final_E
	 ctypes.c_void_p] # debug

# MODEL ###############################################################################

   def init_model(self):
      self.model_lib.ll_init_model_()

   def init_config(self, at, Emax):
      n = ctypes.c_int(len(at))
      Emax = ctypes.c_double(Emax)
      self.model_lib.ll_init_config_(ctypes.byref(n), at.get_positions(), at.get_cell(), ctypes.byref(Emax))

   def eval_energy(self, at):
      n = ctypes.c_int(len(at))
      if 'ns_extra_data' in at.arrays:
	 n_extra_data_c = ctypes.c_int(at.arrays['ns_extra_data'].shape[1])
	 extra_data = at.arrays['ns_extra_data']
      else:
	 n_extra_data_c = ctypes.c_int(0)
	 extra_data=np.zeros( (1) )
      return self.model_lib.ll_eval_energy_(ctypes.byref(n), at.get_positions(), ctypes.byref(n_extra_data_c), extra_data, at.get_cell())

   def move_atom_1(self, at, d_i, d_pos, dEmax):
      n = ctypes.c_int(len(at))
      if 'ns_extra_data' in at.arrays:
	 n_extra_data_c = ctypes.c_int(at.arrays['ns_extra_data'].shape[1])
	 extra_data = at.arrays['ns_extra_data']
      else:
	 n_extra_data_c = ctypes.c_int(0)
	 extra_data=np.zeros( (1) )
      d_i = ctypes.c_int(d_i)
      dEmax = ctypes.c_double(dEmax)
      dE = np.zeros( (1) )
      pos = at.get_positions()
      accept = self.model_lib.ll_move_atom_1_(ctypes.byref(n), pos, ctypes.byref(n_extra_data_c), extra_data,
					     at.get_cell(), ctypes.byref(d_i), d_pos, ctypes.byref(dEmax), dE)
      if accept > 0:
	 at.set_positions(pos)
	 if n_extra_data_c.value > 0:
	    at.arrays['ns_extra_data'][:,:] = extra_data
      return (accept, dE[0])

   def eval_forces(self, at, forces):
      n = ctypes.c_int(len(at))
      if 'ns_extra_data' in at.arrays:
	 n_extra_data_c = ctypes.c_int(at.arrays['ns_extra_data'].shape[1])
	 extra_data = at.arrays['ns_extra_data']
      else:
	 n_extra_data_c = ctypes.c_int(0)
	 extra_data=np.zeros( (1) )
      return self.model_lib.ll_eval_forces_(ctypes.byref(n), at.get_positions(), ctypes.byref(n_extra_data_c), extra_data, at.get_cell(), forces)

# MC/MD WALK ###############################################################################

   def seed_size(self):
      return self.lib.fortran_seed_size_()

   def set_seed(self, seed):
      n_seed = ctypes.c_int(len(seed))
      self.lib.fortran_set_seed_(ctypes.byref(n_seed), seed)

   def MC_atom_walk_velo(self, at, n_steps, step_size, KEmax):
      n = ctypes.c_int(len(at))
      n_steps = ctypes.c_int(n_steps)
      step_size = ctypes.c_double(step_size)
      KEmax = ctypes.c_double(KEmax)
      velo = at.get_velocities()
      n_accept = np.zeros( (1), dtype=np.int32)
      final_KE = np.zeros( (1), dtype=np.float64) 
      self.lib.fortran_mc_atom_velo_(ctypes.byref(n), velo, at.get_masses(),
	 ctypes.byref(n_steps), ctypes.byref(step_size),
	 ctypes.byref(KEmax), final_KE, n_accept)
      at.set_velocities(velo)
      return (n_accept[0], final_KE[0])

   def MC_atom_walk(self, at, n_steps, step_size_pos, Emax, step_size_velo=None):
      n = ctypes.c_int(len(at))
      n_steps = ctypes.c_int(n_steps)
      step_size_pos = ctypes.c_double(step_size_pos)
      if step_size_velo is None:
	 step_size_velo_c = ctypes.c_double(0.0)
      else:
	 step_size_velo_c = ctypes.c_double(step_size_velo)
      Emax = ctypes.c_double(Emax)
      pos = at.get_positions()
      if step_size_velo is None:
	 velo = np.zeros( (1), dtype=np.float64 )
      else:
	 velo = at.get_velocities()
      n_accept_pos = np.zeros( (1), dtype=np.int32)
      n_accept_velo = np.zeros( (1), dtype=np.int32)
      final_E = np.zeros( (1), dtype=np.float64) 
      if 'ns_extra_data' in at.arrays:
	 n_extra_data_c = ctypes.c_int(at.arrays['ns_extra_data'].shape[1])
	 extra_data = at.arrays['ns_extra_data']
      else:
	 n_extra_data_c = ctypes.c_int(0)
	 extra_data=np.zeros( (1) )
      self.lib.fortran_mc_atom_(ctypes.byref(n), pos, velo, at.get_masses(), ctypes.byref(n_extra_data_c), 
	 extra_data, at.get_cell(),
	 ctypes.byref(n_steps), ctypes.byref(step_size_pos), ctypes.byref(step_size_velo_c),
	 ctypes.byref(Emax), final_E, n_accept_pos, n_accept_velo)
      at.set_positions(pos)
      if n_extra_data_c.value > 0:
	 at.arrays['ns_extra_data'][:,:] = extra_data
      if step_size_velo is None:
	 return (n_accept_pos[0], final_E[0])
      else:
	 at.set_velocities(velo)
	 return (n_accept_pos[0], n_accept_velo[0], final_E[0])

   def MD_atom_NVE_walk(self, at, n_steps, timestep, debug):
      n = ctypes.c_int(len(at))
      n_steps = ctypes.c_int(n_steps)
      timestep = ctypes.c_double(timestep)
      pos = at.get_positions()
      vel = at.get_velocities()
      debug = ctypes.c_int(debug)
      final_E = np.zeros( (1), dtype=np.float64) 
      if 'ns_extra_data' in at.arrays:
	 n_extra_data_c = ctypes.c_int(at.arrays['ns_extra_data'].shape[1])
	 extra_data = at.arrays['ns_extra_data']
      else:
	 n_extra_data_c = ctypes.c_int(0)
	 extra_data= np.zeros( (1) )
      n_accept = self.lib.fortran_md_atom_nve_(ctypes.byref(n), 
	 pos, vel, at.get_masses(), ctypes.byref(n_extra_data_c), extra_data, at.get_cell(),
	 ctypes.byref(n_steps), ctypes.byref(timestep), final_E, ctypes.byref(debug))
      at.set_positions(pos)
      at.set_velocities(vel)
      if n_extra_data_c.value > 0:
	 at.arrays['ns_extra_data'][:,:] = extra_data
      return final_E[0]

import ctypes
from numpy.ctypeslib import ndpointer
import os, ase

class fortran_MC_MD:
   def __init__(self, model_so):
      self.model_lib = ctypes.CDLL(os.path.dirname(__file__)+"/"+model_so, mode=ctypes.RTLD_GLOBAL)

      self.model_lib.ll_init_config_.argtypes = [ctypes.c_void_p, # N
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # pos
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # cell
	 ctypes.c_void_p ] # Emax

      self.model_lib.ll_eval_energy_.restype = ctypes.c_double
      self.model_lib.ll_eval_energy_.argtypes = [ctypes.c_void_p, # N
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # pos
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")] # cell

      self.model_lib.ll_eval_denergy_1_.restype = ctypes.c_double
      self.model_lib.ll_eval_denergy_1_.argtypes = [ctypes.c_void_p, # N
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # pos
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # cell
	 ctypes.c_void_p, # d_i
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")] # d_pos

      self.model_lib.ll_eval_forces_.restype = ctypes.c_double
      self.model_lib.ll_eval_forces_.argtypes = [ctypes.c_void_p, # N
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # pos
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # cell
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")] # forces

      self.lib = ctypes.CDLL(os.path.dirname(__file__)+"/fortran_MC_MD.so")

      self.lib.fortran_mc_atom_.restype = ctypes.c_int
      self.lib.fortran_mc_atom_.argtypes = [ctypes.c_void_p, # N
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # pos
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # cell
	 ctypes.c_void_p, # n_steps
	 ctypes.c_void_p, # step_size
	 ctypes.c_void_p, # Emax
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")] # final_E

      self.lib.fortran_md_atom_nve_.argtypes = [ctypes.c_void_p, # N
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # pos
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # vel
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # mass
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # cell
	 ctypes.c_void_p, # n_steps
	 ctypes.c_void_p, # timestep
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")] # final_E

   def init_model(self):
      self.model_lib.ll_init_model_()

   def init_config(self, at, Emax):
      n = ctypes.c_int(len(at))
      Emax = ctypes.c_double(Emax)
      self.model_lib.ll_init_config_(ctypes.byref(n), at.get_positions(), at.get_cell(), ctypes.byref(Emax))

   def eval_energy(self, at):
      n = ctypes.c_int(len(at))
      return self.model_lib.ll_eval_energy_(ctypes.byref(n), at.get_positions(), at.get_cell())

   def eval_denergy_1(self, at, d_i, d_pos):
      n = ctypes.c_int(len(at))
      d_i = ctypes.c_int(d_i)
      return self.model_lib.ll_eval_denergy_1_(ctypes.byref(n), at.get_positions(), at.get_cell(), ctypes.byref(d_i), d_pos)

   def eval_forces(self, at, forces):
      n = ctypes.c_int(len(at))
      return self.model_lib.ll_eval_forces_(ctypes.byref(n), at.get_positions(), at.get_cell(), forces)

   def MC_atom_walk(self, at, n_steps, step_size, Emax, final_E):
      n = ctypes.c_int(len(at))
      n_steps = ctypes.c_int(n_steps)
      step_size = ctypes.c_double(step_size)
      Emax = ctypes.c_double(Emax)
      pos = at.get_positions()
      n_accept = self.lib.fortran_mc_atom_(ctypes.byref(n), 
	 pos, at.get_cell(),
	 ctypes.byref(n_steps), ctypes.byref(step_size), ctypes.byref(Emax), final_E)
      at.set_positions(pos)
      return n_accept

   def MD_atom_NVE_walk(self, at, n_steps, timestep, final_E):
      n = ctypes.c_int(len(at))
      n_steps = ctypes.c_int(n_steps)
      timestep = ctypes.c_double(timestep)
      pos = at.get_positions()
      vel = at.get_velocities()
      n_accept = self.lib.fortran_md_atom_nve_(ctypes.byref(n), 
	 pos, vel, at.get_masses(), at.get_cell(),
	 ctypes.byref(n_steps), ctypes.byref(timestep), final_E)
      at.set_positions(pos)
      at.set_velocities(vel)

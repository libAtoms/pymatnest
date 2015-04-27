import ctypes
from numpy.ctypeslib import ndpointer
import os, ase

class fortran_MC_MD:
   def __init__(self, model_so):
      self.model_lib = ctypes.CDLL(os.path.dirname(__file__)+"/"+model_so, mode=ctypes.RTLD_GLOBAL)
      self.model_lib.ll_eval_energy_.restype = ctypes.c_double
      self.model_lib.ll_eval_energy_.argtypes = [ctypes.c_void_p, # n
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # pos
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")] # cell
      self.model_lib.ll_eval_energy_1_.restype = ctypes.c_double
      self.model_lib.ll_eval_energy_1_.argtypes = [ctypes.c_void_p, # N
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # pos
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # cell
	 ctypes.c_void_p, # d_i
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")] # d_pos

      self.lib = ctypes.CDLL(os.path.dirname(__file__)+"/fortran_MC_MD.so")
      self.lib.fortran_mc_.restype = ctypes.c_int
      self.lib.fortran_mc_.argtypes = [ctypes.c_void_p, # N
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # pos
	 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), # cell
	 ctypes.c_void_p, # n_steps
	 ctypes.c_void_p, # step_size
	 ctypes.c_void_p] # Emax

   def eval_energy(self, at):
      print "eval_energy calling ll_eval_energy_()"
      n = ctypes.c_int(len(at))
      return self.model_lib.ll_eval_energy_(ctypes.byref(n), at.get_positions(), at.get_cell())

   def MC_walk(self, at, n_steps, step_size, Emax):
      n = ctypes.c_int(len(at))
      n_steps = ctypes.c_int(n_steps)
      step_size = ctypes.c_double(step_size)
      Emax = ctypes.c_double(Emax)
      return self.lib.fortran_mc_(ctypes.byref(n), 
	 at.get_positions(), at.get_cell(),
	 ctypes.byref(n_steps), ctypes.byref(step_size), ctypes.byref(Emax))


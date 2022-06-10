#!/usr/bin/env python

import os, sys, ase, ase.io
#import quippy
#from quippy import *
import ctypes
from numpy.ctypeslib import ndpointer
import numpy as np
import fortranMCMDpy

# chnage the model type here
model_name = 'example_Jaglasmooth_model.so'

n_arg = len(sys.argv)

at_this = ase.io.read(sys.argv[1],format="extxyz")

f_MC_MD = fortranMCMDpy.fortran_MC_MD(model_name)
#params = np.array([0.0])
params = np.array([float(_s) for _s in sys.argv[2:5]])
print(params)

f_MC_MD.init_model(params)

energy = f_MC_MD.eval_energy(at_this)
volume = at_this.get_volume()
cell=at_this.get_cell()

if (n_arg > 6):
   pressure = float(sys.argv[6])
else:
   pressure = 0.0
#   print("# energy volume enthalpy")

print("EVHA", energy/len(at_this), volume/len(at_this), energy + pressure * volume, volume/len(at_this)/cell[2,2])

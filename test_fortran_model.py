#!/usr/bin/env python

import ase.io
import sys
import numpy as np

if len(sys.argv) != 3:
    sys.stderr.write("Usage: %s model.so model_params < input.xyz\n" % sys.argv[0])
    sys.exit(1)

at = ase.io.read(sys.stdin, format="extxyz")

import fortranMCMDpy
FORTRAN_model = sys.argv[1]
FORTRAN_model_params = sys.argv[2]
f_MC_MD = fortranMCMDpy.fortran_MC_MD(FORTRAN_model)
params = np.array([ float(x) for x in FORTRAN_model_params.split() ])
f_MC_MD.init_model(params)

f0 = np.zeros((len(at),3))
e0 = f_MC_MD.eval_forces(at, f0)
print "e0 ", e0
print "f0 ", f0
print ""
pos_0 = at.get_positions()
for i_at in range(len(at)):
    for i_cart in range(3):

        pos_pert = pos_0.copy()
        for i_dx in range(8):
            dx = 10.0**(-i_dx)

            pos_pert[i_at, i_cart] = pos_0[i_at, i_cart] + dx
            at.set_positions(pos_pert)
            ep = f_MC_MD.eval_energy(at)

            pos_pert[i_at, i_cart] = pos_0[i_at, i_cart] - dx
            at.set_positions(pos_pert)
            em = f_MC_MD.eval_energy(at)
            print i_at, i_cart, dx, (ep-em)/(2.0*dx), f0[i_at,i_cart], (ep-em)/(2.0*dx)+ f0[i_at,i_cart]
        print ""

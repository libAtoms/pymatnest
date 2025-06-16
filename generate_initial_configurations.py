#!/usr/bin/env python

import sys

import numpy as np
import ase, ase.io
from ase.units import kB

RNG = np.random.default_rng(813756) # 1 = 12345, 2 = 5, 3 = 37
print(RNG)
KEMAX_MAX_T = 10000


def main(n_walkers=int(sys.argv[4]) , at="He", num_adatom = int(sys.argv[3])):
    slabs = []
    for _ in range(n_walkers):
        slab=ase.io.read(filename=sys.argv[1] ,format="extxyz")
        cell=slab.get_cell()

        for n in range(num_adatom):
            new_adatom = np.zeros(3)
            new_adatom[0] = np.random.rand(1)*cell[0,0]
            new_adatom[1] = np.random.rand(1)*cell[1,1]
            new_adatom[2] = np.random.rand(1)*(cell[2,2]-19.92125615) + 9.92125615

            slab.append(ase.Atom(at,position=new_adatom))
        
        masses = np.repeat(1, slab.positions.shape[0])
        masses[-num_adatom:] = np.repeat(1, num_adatom)
        slab.set_masses(masses)
        slab.set_momenta(np.zeros(slab.positions.shape))
        slab.info["KEmax"] = 1.5 * len(slab) * kB * KEMAX_MAX_T
        slab.info["iter"] = -1 
        slab.info["volume"] = slab.get_volume()
        slabs.append(slab)
    ase.io.write(images=slabs, filename=sys.argv[2], format="extxyz")

if __name__ == "__main__":
    main()

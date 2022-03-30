""" Generate (111) surfaces of fcc LJ solids """
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from ase.build import fcc111, add_adsorbate
from ase.io import write
from ase.units import kB
from ase.visualize.plot import plot_atoms

RNG = np.random.default_rng(5)  # 1 = 12345, 2 = 5, 3 = 37
PLT = False
KEMAX_MAX_T = 10000
set_wall = True


# Direction (x, y, z) normal to which a reflecting wall is set at fractional
# positions of 0 and 1 along that direction. Should only be used for surface
# slabs, centered in the simulation cell


def walkers(n_walkers=240, at="H", size=(2, 2, 3), sig=2.5, v=10, nat=4, buf=2):
    slabs = []
    for _ in range(n_walkers):
        nnd = 2 ** (1 / 6) * sig
        a = nnd * np.sqrt(2)
        slab = fcc111(at, size, a=a, vacuum=v, orthogonal=True, periodic=False)
        a = slab.cell[0][0]
        b = slab.cell[1][1]
        for i in range(nat):
            h = RNG.uniform(low=buf, high=v - buf)  # top surface
            p = (RNG.uniform(high=a), RNG.uniform(high=b))
            add_adsorbate(slab=slab, adsorbate=at, height=h, position=p)
        slab.set_masses(np.repeat(1, slab.positions.shape[0]))
        slab.set_momenta(np.zeros(slab.positions.shape))
        slab.info["KEmax"] = 1.5 * len(slab) * kB * KEMAX_MAX_T
        slab.info["iter"] = -1
        slab.info["volume"] = slab.get_volume()
        slab.info["set_wall"] = set_wall
        slabs.append(slab)
    return slabs


def main():
    # nats = np.arange(12) + 1
    nats = [8]
    n_walkers_per_particle = 28  # * 2
    n_iter = 50000
    for i, nat in enumerate(nats):
        # dir_name = "{:02d}-LJ_111-{}p-{}wpp-{}i".format(
        #    i + 1, nat, n_walkers_per_particle, n_iter)
        # dir_name = "15-double-wpp/08-LJ_111-8p-28wpp-50000i"
        dir_name = "16-change-seed/08-LJ_111-8p-28wpp-50000i"
        # os.system("mkdir -p " + dir_name)
        size = (2, 2, 3)
        keep_atoms_fixed = np.prod(size)
        start_species = keep_atoms_fixed + nat
        n_walkers = n_walkers_per_particle * nat
        n_iter_times_fraction_killed = int(np.ceil(n_iter / n_walkers))
        with open(dir_name + "/input", "w") as fnew:
            with open("input") as fold:
                for line in fold:
                    if "start_species" in line:
                        fnew.write("start_species=1 {}\n".format(start_species))
                    elif "keep_atoms_fixed" in line:
                        fnew.write("keep_atoms_fixed={}\n".format(
                            keep_atoms_fixed))
                    elif "n_walkers" in line:
                        fnew.write("n_walkers={}\n".format(n_walkers))
                    elif "n_iter_times_fraction_killed" in line:
                        fnew.write("n_iter_times_fraction_killed={}\n".format(
                            n_iter_times_fraction_killed))
                    else:
                        fnew.write(line)
        slabs = walkers(n_walkers=n_walkers, at="H", size=size, sig=2.5, v=20,
                        nat=nat, buf=2)
        write(filename=dir_name + "/slab.extxyz", images=slabs, format="extxyz")


if __name__ == "__main__":
    main()

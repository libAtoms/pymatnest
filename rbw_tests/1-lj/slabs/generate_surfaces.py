""" Generate (111) surfaces of fcc LJ solids """
import sys

import matplotlib.pyplot as plt
import numpy as np
from ase.build import fcc111, add_adsorbate
from ase.io import write
from ase.units import kB
from ase.visualize.plot import plot_atoms

RNG = np.random.default_rng(12345)  # 1 = 12345, 2 = 5, 3 = 37
PLT = False
KEMAX_MAX_T = 10000


def main(n_walkers=240, at="H", size=(2, 2, 3), sig=2.5, v=10, nat=4, buf=2):
    slabs = []
    for _ in range(n_walkers):
        nnd = 2 ** (1 / 6) * sig
        a = nnd * np.sqrt(2)
        slab = fcc111(at, size, a=a, vacuum=v, orthogonal=True, periodic=False)
        a = slab.cell[0][0]
        b = slab.cell[1][1]
        zlo = slab.positions[:, 2].min()
        zhi = slab.positions[:, 2].max()
        # for _ in range(nat):
        #    height = RNG.uniform(low=buf, high=vac * 2 - buf)
        #    pos = (RNG.uniform(high=a), RNG.uniform(high=b))
        #    add_adsorbate(slab=slab, adsorbate=at, height=height, position=pos)
        # slab.wrap()
        for i in range(nat):
            h_top = RNG.uniform(low=buf, high=v - buf)  # top surface
            h_bot = RNG.uniform(low=-zhi + buf, high=zlo - zhi - buf)  # bottom
            p_top = (RNG.uniform(high=a), RNG.uniform(high=b))
            p_bot = (RNG.uniform(high=a), RNG.uniform(high=b))
            add_adsorbate(slab=slab, adsorbate=at, height=h_top, position=p_top)
            add_adsorbate(slab=slab, adsorbate=at, height=h_bot, position=p_bot)
        slab.set_masses(np.repeat(1, slab.positions.shape[0]))
        slab.set_momenta(np.zeros(slab.positions.shape))
        slab.info["KEmax"] = 1.5 * len(slab) * kB * KEMAX_MAX_T
        slab.info["iter"] = -1
        slab.info["volume"] = slab.get_volume()
        slabs.append(slab)
    write(filename="slab-1.extxyz", images=slabs, format="extxyz")
    if PLT:
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(7, 7))
        plot_atoms(slabs[10], axs[0, 0], radii=0.7, rotation='90x,90y,0z')
        plot_atoms(slabs[20], axs[0, 1], radii=0.7, rotation='90x,90y,0z')
        plot_atoms(slabs[30], axs[1, 0], radii=0.7, rotation='90x,90y,0z')
        plot_atoms(slabs[40], axs[1, 1], radii=0.7, rotation='90x,90y,0z')
        fig.savefig("slab.png")


if __name__ == "__main__":
    main()

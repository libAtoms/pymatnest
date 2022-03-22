""" Generate clean (111) surface of fcc LJ solid """

import numpy as np
from ase.build import fcc111
from ase.calculators.lj import LennardJones
from ase.io import write
from ase.units import kB
from ase.visualize.plot import plot_atoms
from matplotlib import pyplot as plt

PLT = False
KEMAX_MAX_T = 10000


def main(at="H", size=(2, 2, 5), sig=2.5, vac=10):
    nnd = 2 ** (1 / 6) * sig
    a = nnd * np.sqrt(2)
    vacuum = vac - nnd * np.sqrt(2) / np.sqrt(3)
    slab = fcc111(at, size, a=a, vacuum=vacuum, orthogonal=True, periodic=True)
    slab.wrap()
    slab.set_masses(np.repeat(1, slab.positions.shape[0]))
    slab.set_momenta(np.zeros(slab.positions.shape))
    slab.info["KEmax"] = 1.5 * len(slab) * kB * KEMAX_MAX_T
    slab.info["iter"] = -1
    slab.info["volume"] = slab.get_volume()
    write(filename="clean.extxyz", images=slab, format="extxyz")
    if PLT:
        fig, ax = plt.subplots(figsize=(3.33, 3.33))
        plot_atoms(slab, ax, radii=0.7, rotation='90x,90y,0z')
        fig.savefig("clean.png")

    # calculate potential energy
    calc = LennardJones(sigma=2.5, epsilon=0.1)
    slab.set_calculator(calc)
    print(slab.get_potential_energy())


if __name__ == "__main__":
    main()

""" Calculate order parameters from output trajectory """

from ase.io import read
import numpy as np


def main():
    traj = read("output.traj.0.extxyz", ":")
    U_751 = -39.4586
    U_1251 = -34.6206
    for i, x in enumerate(traj):
        d_751 = np.abs(U_751 - x.info["ns_energy"])
        d_1251 = np.abs(U_1251 - x.info["ns_energy"])
        if d_751 < 1E-2:
            print(i)
        if d_1251 < 1.1E-2:
            print(i)
    i_x_751 = 709
    i_x_1251 = 458
    x_751 = traj[i_x_751]
    x_1251 = traj[i_x_1251]


if __name__ == "__main__":
    main()

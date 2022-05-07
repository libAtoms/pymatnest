""" Analyze nested sampling simulations """
import os
import sys
from glob import glob

import matplotlib.pyplot as plt
import numpy as np

M = 1
n = 200
D = 10
PREFIX = "../../../ns_analyse "
SUFFIX = " -M {} -n {} -D {}".format(M, n, D)


def main():
    paths = sorted(glob("../rotate-bug/2-fixed/**/output.energies", recursive=True))
    paths_w = sorted(glob("../rotate-bug/3-fixed-wall/**/output.energies", recursive=True))
    paths_2 = sorted(glob("../rotate-bug/4-doubl-walk/**/output.energies", recursive=True))

    # coverage analysis
    fig, axs = plt.subplots(ncols=3, figsize=(7, 7 / 3), sharey=True)
    for path, path_w, path_2 in zip(paths[:9], paths_w, paths_2):
        iax = int(path.split("-")[1][-1]) - 1

        # without walls
        outp = os.popen(PREFIX + path + SUFFIX).read().split("\n")[2:-1]
        outa = np.array([x.split() for x in outp])[:, [0, 4]].astype(float)
        temp = outa[:, 0]
        cv = outa[:, 1]
        axs[iax].plot(temp, cv)

        # with walls
        outp_w = os.popen(PREFIX + path_w + SUFFIX).read().split("\n")[2:-1]
        outa_w = np.array([x.split() for x in outp_w])[:, [0, 4]].astype(float)
        temp_w = outa_w[:, 0]
        cv_w = outa_w[:, 1]
        axs[iax].plot(temp_w, cv_w, ls="--")

        # double walkers
        outp_2 = os.popen(PREFIX + path_2 + SUFFIX).read().split("\n")[2:-1]
        outa_2 = np.array([x.split() for x in outp_2])[:, [0, 4]].astype(float)
        temp_2 = outa_2[:, 0]
        cv_2 = outa_2[:, 1]
        axs[iax].plot(temp_2, cv_2, ls=":")
    axs[0].set_ylabel("Heat capacity (eV/K)")
    axs[0].set_xlabel("Temperature (K)")
    axs[1].set_xlabel("Temperature (K)")
    axs[2].set_xlabel("Temperature (K)")
    axs[0].set_title("(a) Coverage = 8/8")
    axs[1].set_title("(b) Coverage = 7/8")
    axs[2].set_title("(c) Coverage = 9/8")
    plt.tight_layout()
    plt.savefig("coverage.png", dpi=300)

    # vacuum analysis
    fig, axs = plt.subplots(ncols=3, figsize=(7, 7 / 3), sharey=True)
    for i, (path, path_w, path_2) in enumerate(zip(
            paths[9:18], paths_w[9:18], paths_2[9:18])):
        iax = int(i / 3)

        # without walls
        outp = os.popen(PREFIX + path + SUFFIX).read().split("\n")[2:-1]
        outa = np.array([x.split() for x in outp])[:, [0, 4]].astype(float)
        temp = outa[:, 0]
        cv = outa[:, 1]
        axs[iax].plot(temp, cv)

        # with walls
        outp_w = os.popen(PREFIX + path_w + SUFFIX).read().split("\n")[2:-1]
        outa_w = np.array([x.split() for x in outp_w])[:, [0, 4]].astype(float)
        temp_w = outa_w[:, 0]
        cv_w = outa_w[:, 1]
        axs[iax].plot(temp_w, cv_w, ls="--")

        # double walkers
        outp_2 = os.popen(PREFIX + path_2 + SUFFIX).read().split("\n")[2:-1]
        outa_2 = np.array([x.split() for x in outp_2])[:, [0, 4]].astype(float)
        temp_2 = outa_2[:, 0]
        cv_2 = outa_2[:, 1]
        axs[iax].plot(temp_2, cv_2, ls=":")
    axs[0].set_ylabel("Heat capacity (eV/K)")
    axs[0].set_xlabel("Temperature (K)")
    axs[1].set_xlabel("Temperature (K)")
    axs[2].set_xlabel("Temperature (K)")
    axs[0].set_title("(a) Vacuum = 20 Å")
    axs[1].set_title("(b) Vacuum = 10 Å")
    axs[2].set_title("(c) Vacuum = 40 Å")
    plt.tight_layout()
    plt.savefig("vacuum.png", dpi=300)

    # area analysis
    fig, ax = plt.subplots(figsize=(3.33, 3.33))

    # without walls
    for i, path in enumerate(paths[18:]):
        outp = os.popen(PREFIX + path + SUFFIX).read().split("\n")[2:-1]
        outa = np.array([x.split() for x in outp])[:, [0, 4]].astype(float)
        temp = outa[:, 0]
        cv = outa[:, 1]
        ax.plot(temp, cv)

    # with walls
    for i, path in enumerate(paths_w[18:]):
        outp_w = os.popen(PREFIX + path + SUFFIX).read().split("\n")[2:-1]
        outa_w = np.array([x.split() for x in outp_w])[:, [0, 4]].astype(float)
        temp_w = outa_w[:, 0]
        cv_w = outa_w[:, 1]
        ax.plot(temp_w, cv_w, ls="--")
    ax.set_ylabel("Heat capacity (eV/K)")
    ax.set_xlabel("Temperature (K)")
    plt.tight_layout()
    plt.savefig("area.png", dpi=300)


if __name__ == "__main__":
    main()

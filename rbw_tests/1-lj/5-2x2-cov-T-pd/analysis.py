""" Analyze nested sampling simulations """
import os
import sys
from glob import glob

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

M = 1
n = 200
D = 10
PREFIX = "/Users/robertwexler/PycharmProjects/pymatnest/ns_analyse "
SUFFIX = " -M {} -n {} -D {}".format(M, n, D)


def main():
    paths = sorted(glob("*/output.energies"))
    paths_double = sorted(glob("13-double-n_iter/*/output.energies"))
    paths_triple = sorted(glob("14-triple-n_iter/*/output.energies"))
    paths_walker = sorted(glob("15-double-wpp/*/output.energies"))
    paths_rseeds = sorted(glob("16-change-seed/*/output.energies"))
    fig, axs = plt.subplots(figsize=(3.33, 3.33))
    cmap = cm.get_cmap("tab20c", 20)
    colors = [cmap(i) for i in range(len(paths))]
    j = 0
    min_cv = []
    for i, path in enumerate(paths):
        # 50K iterations
        nat = int(path.split("-")[2][:-1])
        if False:
            wpp = int(path.split("-")[3][:-3])
            n_iter = int(path.split("-")[4].split("/")[0][:-1])
            outp = os.popen(PREFIX + path + SUFFIX).read().split("\n")[2:-1]
            outa = np.array([x.split() for x in outp])[:, [0, 4]].astype(float)
            temp = outa[:, 0]
            cv = outa[:, 1] / nat * 1000
            axs.plot(temp,
                     cv + nat,
                     color="k",
                     alpha=0.25,
                     lw=1,
                     label=nat,
                     zorder=100 - nat)

        # 100K iterations
        if nat < 11:
            path_double = paths_double[i]
            outp_double = os.popen(
                PREFIX + path_double + SUFFIX).read().split("\n")[2:-1]
            outa_double = np.array(
                [x.split() for x in outp_double])[:, [0, 4]].astype(float)
            temp_double = outa_double[:, 0]
            cv_double = outa_double[:, 1] / nat * 1000
            axs.plot(temp_double,
                     cv_double + nat,
                     color=colors[i],
                     # color="k",
                     # ls="--",
                     # alpha=0.5,
                     # lw=1,
                     label=nat,
                     zorder=100 - nat)
            min_cv.append(np.min(cv_double) + nat)

        # 150K iterations
        if nat >= 11:
            path_triple = paths_triple[j]
            outp_triple = os.popen(
                PREFIX + path_triple + SUFFIX).read().split("\n")[2:-1]
            outa_triple = np.array(
                [x.split() for x in outp_triple])[:, [0, 4]].astype(float)
            temp_triple = outa_triple[:, 0]
            cv_triple = outa_triple[:, 1] / nat * 1000
            axs.plot(temp_triple,
                     cv_triple + nat,
                     color=colors[i],
                     # color="k",
                     # ls=":",
                     # alpha=1,
                     # lw=1,
                     label=nat,
                     zorder=100 - nat)
            j += 1
            min_cv.append(np.min(cv_triple) + nat)

        # double number of walkers per particle for two monolayers
        if nat == 8:
            path_walker = paths_walker[j]
            outp_walker = os.popen(
                PREFIX + path_walker + SUFFIX).read().split("\n")[2:-1]
            outa_walker = np.array(
                [x.split() for x in outp_walker])[:, [0, 4]].astype(float)
            temp_walker = outa_walker[:, 0]
            cv_walker = outa_walker[:, 1] / nat * 1000
            axs.plot(temp_walker,
                     cv_walker + nat,
                     # color=colors[i],
                     color="k",
                     ls="--",
                     alpha=0.5,
                     # lw=1,
                     label=nat,
                     zorder=100 - nat)

        # change random seed for two monolayers
        if nat == 8:
            path_rseeds = paths_rseeds[j]
            outp_rseeds = os.popen(
                PREFIX + path_rseeds + SUFFIX).read().split("\n")[2:-1]
            outa_rseeds = np.array(
                [x.split() for x in outp_rseeds])[:, [0, 4]].astype(float)
            temp_rseeds = outa_rseeds[:, 0]
            cv_rseeds = outa_rseeds[:, 1] / nat * 1000
            axs.plot(temp_rseeds,
                     cv_rseeds + nat,
                     # color=colors[i],
                     color="k",
                     ls=":",
                     alpha=0.5,
                     # lw=1,
                     label=nat,
                     zorder=100 - nat)

    axs.set_xlabel("Temperature (K)")
    # axs.set_ylabel("Heat capacity (meV K$^{-1}$ particle$^{-1}$)")
    axs.set_ylabel("Number of free particles")
    axs.set_yticks(min_cv)
    axs.set_yticklabels(np.arange(1, 13))
    # plt.legend(fontsize="small")
    plt.tight_layout()
    plt.savefig("pd-2.png", dpi=300)


if __name__ == "__main__":
    main()

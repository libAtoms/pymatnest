""" Analyze nested sampling simulations """
import os
import sys
from glob import glob

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm

M = 1
n = 2000
D = 1
PREFIX = "/Users/robertwexler/PycharmProjects/pymatnest/ns_analyse "
SUFFIX = " -M {} -n {} -D {}".format(M, n, D)


def main():
    paths = sorted(glob("*/output.energies"))
    fig, axs = plt.subplots(figsize=(3.33, 3.33))
    cmap = cm.get_cmap("tab20c", 20)
    colors = [cmap(i) for i in range(len(paths))]
    df = pd.DataFrame()
    for i, path in enumerate(paths):
        nat = int(path.split("-")[2][:-1])
        outp = os.popen(PREFIX + path + SUFFIX).read().split("\n")[2:-1]
        outa = np.array([x.split() for x in outp])[:, [0, 4]].astype(float)
        temp = outa[:, 0]
        cv = outa[:, 1] / nat * 1000
        axs.plot(temp,
                 cv + nat / 4,
                 color=colors[i],
                 label=nat,
                 zorder=100 - nat)
        if i == 0:
            df["temp"] = temp
        df[nat] = cv
    df.to_csv("cv_vs_temp.csv", index=False)
    axs.set_xlabel("Temperature (K)")
    axs.set_ylabel("Number of free particles")
    axs.set_yticks([])
    plt.tight_layout()
    plt.savefig("pd.png", dpi=300)


if __name__ == "__main__":
    main()

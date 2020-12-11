"""
Plotting functions.
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

# to access from inside and outside snakemake
rc = "matplotlibrc"
if not os.path.isfile(rc):
    rc = "../" + rc

plt.style.use(["bmh", rc])


__author__ = "Fabian Neumann (KIT)"
__copyright__ = f"Copyright 2020 {__author__}, GNU GPL 3"


def plot_histograms(truth, predictions, fn=None):

    if not isinstance(predictions, list):
        predictions = [predictions]

    fig, axes = plt.subplots(3, 3, figsize=(8, 8))
    for i, c in enumerate(truth.columns):
        ax = axes[int(i / 3)][i % 3]
        ax.set_title(c)
        ax.set_ylim([0, 0.1])
        bins = np.arange(-100, 101, 5)
        if c == "tsc":
            ax.set_xlim([-1.5, 1.5])
            ax.set_ylim([0, 3])
            ax.set_xlabel("bn EUR p.a.")
            bins = np.arange(-1.5, 1.6, 0.1)
        elif c == "transmission":
            ax.set_xlim([-25, 25])
            ax.set_xlabel("TWkm")
            bins = np.arange(-25, 26, 1)
        elif c in ["battery", "H2"]:
            ax.set_xlim([-25, 25])
            bins = np.arange(-25, 26, 1)
        elif c == "gini":
            ax.set_xlim([-0.1, 0.1])
            ax.set_xlabel("gini index [-]")
            ax.set_ylim([0, 100])
            bins = np.arange(-0.1, 0.11, 0.01)
        else:
            ax.set_xlim([-100, 100])
            ax.set_xlabel("GW")
        for j, p in enumerate(predictions):
            lookup = {0: "train", 1: "test"}
            if j in lookup.keys():
                j = lookup[j]
            (p - truth)[c].plot.hist(
                ax=ax, label=f"pred. {j}", alpha=0.4, bins=bins, density=True
            )
        ax.grid(False)
    plt.tight_layout()
    axes[2, 2].legend()
    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")


def plot_sobol(data, fn=None):
    fig, ax = plt.subplots(figsize=(5, 8))
    sns.heatmap(
        data,
        square=True,
        cmap="Blues",
        vmax=1,
        vmin=0,
        annot=True,
        fmt=".2f",
        cbar=False,
    )
    plt.ylabel("inputs")
    plt.xlabel("outputs")
    if fn is not None:
        plt.savefig(fn, bbox_inches="tight")

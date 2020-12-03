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
        ax.set_ylim([0, 0.025])
        bins = np.arange(-200, 201, 10)
        if c == "tsc":
            ax.set_xlim([-5.5, 5.5])
            ax.set_ylim([0, 0.7])
            ax.set_xlabel("bn EUR p.a.")
            bins = np.arange(-5.5, 5.6, 0.25)
        elif c == "transmission":
            ax.set_xlim([-75, 75])
            ax.set_xlabel("TWkm")
        elif c == "gini":
            ax.set_xlim([-0.2, 0.2])
            ax.set_xlabel("gini index [-]")
            ax.set_ylim([0, 20])
            bins = np.arange(-0.2, 0.25, 0.01)
        else:
            ax.set_xlim([-200, 200])
            ax.set_xlabel("GW")
        for j, p in enumerate(predictions):
            (p - truth)[c].plot.hist(
                ax=ax, label=f"pred. {j}", alpha=0.4, bins=bins, density=True
            )
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

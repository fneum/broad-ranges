"""
Helper functions and classes used within multiple snakemake rules.
"""

import pandas as pd
import numpy as np
import chaospy
import numpoly
import re

from sklearn.metrics import explained_variance_score, r2_score

__author__ = "Fabian Neumann (KIT)"
__copyright__ = f"Copyright 2020 {__author__}, GNU GPL 3"


class NamedJ:
    """Dictionary-like wrapper for joint random variable generator with names."""

    def __init__(self, distributions):
        self.J = self.J_from_dict(distributions.values())
        self.names = distributions.keys()
        self.mapping = {k: i for i, k in enumerate(self.names)}

    def __getitem__(self, attr):
        return self.J[self.mapping[attr]]

    def __repr__(self):
        return "\n".join([f"{k}: {self[k]}" for k in self.names])

    def J_from_dict(self, values):
        DD = []
        for v in values:
            D = getattr(chaospy, v["type"])
            DD.append(D(*v["args"]))
        return chaospy.J(*DD)

    def sample(self, size=100, rule="halton", fmt=3):
        samples = self.J.sample(size=size, rule=rule).round(fmt)
        index = [f"{n}-cost" for n in self.names]
        return pd.DataFrame(samples, index=index)


class NamedPoly:
    """Dictionary-like wrapper for vector numpoly polynomials with names."""

    def __init__(self, poly, names):
        self.poly = poly
        self.names = list(names)
        self.mapping = {k: i for i, k in enumerate(self.names)}

    def __getitem__(self, attr):
        return self.poly[self.mapping[attr]]

    def __repr__(self):
        r = numpoly.array_repr
        return "\n\n".join(
            [f"{k}: {r(self.poly[i])}" for i, k in enumerate(self.names)]
        )

    def __call__(self, args):
        return pd.DataFrame(self.poly(*args), index=self.names).squeeze()

    def __add__(self, other):
        assert self.names == other.names, "Names have to match!"
        poly = self.poly + other.poly
        return NamedPoly(poly, self.names)

    def __sub__(self, other):
        assert self.names == other.names, "Names have to match!"
        poly = self.poly - other.poly
        return NamedPoly(poly, self.names)

    def __mul__(self, other):
        assert self.names == other.names, "Names have to match!"
        poly = self.poly * other.poly
        return NamedPoly(poly, self.names)

    def __div__(self, other):
        assert self.names == other.names, "Names have to match!"
        poly = self.poly / other.poly
        return NamedPoly(poly, self.names)

    @classmethod
    def from_txt(cls, fn):
        text = open(fn).read()
        pattern = "\n# [A-Za-z\s\d]*\n"
        names = re.search(pattern, text)[0][3:-1].split(" ")
        poly = numpoly.loadtxt(fn)
        return cls(poly, names)

    def to_txt(self, fn, fmt="%.8f"):
        numpoly.savetxt(fn, self.poly, header=" ".join(self.names), fmt=fmt)


def load_dataset(fn, obj="cost", sense="min", eps=None, fixed="none", pos=None):

    assert not (obj == "cost" and fixed != "none"), "Incompatible choice!"

    raw_df = pd.read_csv(fn, index_col=0, header=list(range(10))).T

    levels = ["objective", "sense", "fixed"]
    df = raw_df.xs([obj, sense, fixed], level=levels, axis=0).copy()

    if fixed == "none":
        df.index = df.index.droplevel(level="position")
    else:
        if pos is not None:
            if not isinstance(pos, str):
                pos = str(pos)
            df = df.xs(pos, level="position")

    if obj == "cost":
        df.index = df.index.droplevel(level="epsilon")
    else:
        df = df.append(raw_df.xs(["cost", "min"], level=["objective", "sense"]))
        if eps is not None:
            if not isinstance(eps, str):
                eps = str(eps)
            df = df.xs(eps, level="epsilon")

    return df


def multiindex2df(multiindex):
    return pd.DataFrame(multiindex2array(multiindex), index=multiindex.names)


def multiindex2array(multiindex):
    return np.array([np.array(row).astype(float) for row in multiindex]).T


def calculate_errors(prediction, truth):
    kws = dict(multioutput="raw_values")
    diff = prediction - truth
    return pd.concat(
        {
            "mape": diff.abs().mean() / truth.mean() * 100,
            "mae": diff.abs().mean(),
            "r2": pd.Series(r2_score(truth, prediction, **kws), index=truth.columns),
            "variance_explained": pd.Series(
                explained_variance_score(truth, prediction, **kws), index=truth.columns
            ),
        },
        axis=1,
    )


def build_pce_prediction(model, samples):
    prediction = model(samples.values).clip(lower=0.0)
    prediction.columns = pd.MultiIndex.from_frame(samples.astype(str).T)
    return prediction.T

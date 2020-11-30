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

    def __call__(self, *args):
        return pd.Series(self.poly(*args), index=self.names)

    def __add__(self, other):
        assert self.names == other.names, "Names have to match!"
        poly = self.poly + other.poly
        return NamedPoly(poly, self.names)

    def __sub__(self, other):
        assert self.names == other.names, "Names have to match!"
        poly = self.poly - other.poly
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


def load_dataset(fn, obj="cost", sense="min"):

    raw_df = pd.read_csv(fn, index_col=0, header=list(range(8))).T

    df = raw_df.xs([obj, sense], level=["objective", "sense"], axis=0).copy()

    if obj == "cost":
        df.index = df.index.droplevel(level="epsilon")
    else:
        df = df.append(raw_df.xs(["cost", "min"], level=["objective", "sense"]))

    return df


# def load_dataset_old(fn):
#     df = pd.read_csv(fn, index_col=0, header=list(range(5))).T

#     df["offwind"] = df["offwind-ac"] + df["offwind-dc"]
#     df["wind"] = df["onwind"] + df["offwind"]
#     df["transmission"] = df["lines"] + df["links"] + 290
#     df.drop(
#         ["ror", "hydro", "PHS", "offwind-ac", "offwind-dc", "lines", "links"],
#         axis=1,
#         inplace=True,
#     )
#     return df


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


def build_ann_prediction(model, samples, mirror):
    prediction = model.predict(samples.T)
    return pd.DataFrame(prediction, index=mirror.index, columns=mirror.columns)


def build_pce_prediction(model, samples):
    prediction = samples.apply(lambda s: model(*s), result_type="expand")
    prediction.columns = pd.MultiIndex.from_frame(samples.astype(str).T)
    return prediction.T

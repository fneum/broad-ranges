"""
Build surrogate model.
"""

from sklearn.model_selection import train_test_split
from sklearn import linear_model as lm
import chaospy
import numpoly

import _helpers as h
import _plotters as p

__author__ = "Fabian Neumann (KIT)"
__copyright__ = f"Copyright 2020 {__author__}, GNU GPL 3"


def build_surrogate(order, distribution, train_set, sklearn_model=None):
    samples = h.multiindex2array(train_set.index)
    pce = chaospy.orth_ttr(order, distribution.J)
    surrogate = chaospy.fit_regression(pce, samples, train_set.values, sklearn_model)
    variables = train_set.columns
    return h.NamedPoly(surrogate, variables)


def build_sklearn_model(cf):
    with_sklearn = cf.get("with_sklearn", None)
    if not with_sklearn:
        return None
    method = with_sklearn.pop("method")
    return getattr(lm, method)(**with_sklearn)


def apply_multifidelity(model, filename, dimension, sense, order, distribution):
    if filename is None:
        return

    hf_dataset = h.load_dataset(filename, dimension, sense)
    hf_samples = h.multiindex2df(hf_dataset.index)
    lf_dataset = h.build_pce_prediction(model, hf_samples)

    scaling_factors = hf_dataset - lf_dataset
    scaling_function = build_surrogate(order, distribution, scaling_factors)

    model += scaling_function


if __name__ == "__main__":

    cf = snakemake.config
    uncertainties = cf["uncertainties"]
    epsilons = cf["scenarios"]["epsilon"]

    order = int(snakemake.wildcards.order)
    dimension = snakemake.wildcards.dimension
    sense = snakemake.wildcards.sense

    dataset = h.load_dataset(snakemake.input["low"], dimension, sense)

    if dimension != "cost":
        uncertainties["epsilon"] = dict(type="Uniform", args=[0, max(epsilons)])
    distribution = h.NamedJ(uncertainties)

    train_set, test_set = train_test_split(dataset, **cf["train_test_split"])

    # Model

    sklearn_model = build_sklearn_model(cf)

    model = build_surrogate(order, distribution, train_set, sklearn_model)

    model.to_txt(snakemake.output.polynomial, fmt="%.4f")

    # Correct with high fidelity model runs

    filename = snakemake.input.get("high", None)
    apply_multifidelity(model, filename, dimension, sense, order, distribution)

    # Evaluation

    train_samples = h.multiindex2df(train_set.index)
    train_predictions = h.build_pce_prediction(model, train_samples)

    test_samples = h.multiindex2df(test_set.index)
    test_predictions = h.build_pce_prediction(model, test_samples)

    p.plot_histograms(
        dataset, [train_predictions, test_predictions], fn=snakemake.output.plot
    )

    h.calculate_errors(train_predictions, train_set).to_csv(
        snakemake.output.train_errors, **cf["csvargs"]
    )
    h.calculate_errors(test_predictions, test_set).to_csv(
        snakemake.output.test_errors, **cf["csvargs"]
    )

"""
[description]
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


if __name__ == "__main__":

    cf = snakemake.config
    order = int(snakemake.wildcards.order)

    dataset = h.load_dataset(snakemake.input[0])
    distribution = h.NamedJ(cf["uncertainties"])

    train_set, test_set = train_test_split(dataset, **cf["train_test_split"])

    # Model

    sklearn_model = build_sklearn_model(cf)

    model = build_surrogate(order, distribution, train_set, sklearn_model)

    model.to_txt(snakemake.output.polynomial, fmt="%.4f")

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

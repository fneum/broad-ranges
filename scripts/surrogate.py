"""
Build surrogate model.
"""

from sklearn.model_selection import train_test_split
from sklearn import linear_model as lm
from numpoly import inner
import chaospy

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


def apply_multifidelity(model, kind, filename, dimension, sense, order, distribution):
    if filename is None:
        return

    hf_dataset = h.load_dataset(filename, dimension, sense)
    hf_samples = h.multiindex2df(hf_dataset.index)
    lf_dataset = h.build_pce_prediction(model, hf_samples)

    def additive_correction():
        scaling_factors = hf_dataset - lf_dataset
        scaling_function = build_surrogate(order, distribution, scaling_factors)
        return scaling_function

    def multiplicative_correction():
        scaling_factors = hf_dataset / lf_dataset
        scaling_function = build_surrogate(order, distribution, scaling_factors)
        return scaling_function

    def optimal_balance(func_a, func_b):
        nominator = inner(func_b, func_b)
        denominator = inner(func_a, func_a) + nominator
        omega = nominator / denominator
        return omega

    if kind == "additive":
        corrected_model = model + additive_correction()

    elif kind == "multiplicative":
        corrected_model = model * multiplicative_correction()

    elif kind == "hybrid":
        add_correction = additive_correction()
        mul_correction = multiplicative_correction()
        omega = optimal_balance(add_correction, mul_correction)
        corrected_model = omega * (model + add_correction) + (1 - omega) * (
            model * mul_correction
        )

    else:
        raise NotImplementedError(
            f"Multifidelity kind must be in ['additive', 'multiplicative', 'hybrid']. Is {kind}"
        )

    return corrected_model


if __name__ == "__main__":

    cf = snakemake.config
    uncertainties = cf["uncertainties"]
    epsilons = cf["nearoptimal"]["epsilon"]

    dimension = snakemake.wildcards.dimension
    sense = snakemake.wildcards.sense

    if dimension == "cost":
        cf_surrogate = cf["surrogate"]["cost"]
    else:
        cf_surrogate = cf["surrogate"]["epsilon"]

    dataset = h.load_dataset(snakemake.input["low"], dimension, sense)

    if dimension != "cost":
        uncertainties["epsilon"] = dict(type="Uniform", args=[0, max(epsilons)])
    distribution = h.NamedJ(uncertainties)

    train_set, test_set = train_test_split(dataset, **cf["train_test_split"])

    # Model

    order = cf_surrogate["order"]

    sklearn_model = build_sklearn_model(cf)
    model = build_surrogate(order, distribution, train_set, sklearn_model)
    model.to_txt(snakemake.output.low_polynomial, fmt="%.4f")

    # Evaluation (based on low-fidelity data)

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

    # Multifidelity

    filename = snakemake.input.get("high", None)
    cf_correction = cf_surrogate["correction"]
    approach = cf_correction["approach"]
    order = cf_correction["order"]

    model = apply_multifidelity(
        model, approach, filename, dimension, sense, order, distribution
    )

    model.to_txt(snakemake.output.high_polynomial, fmt="%.4f")

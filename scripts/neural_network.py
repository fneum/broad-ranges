"""[summary]"""

from sklearn import neural_network as ann
import pickle

import _helpers as h
import _plotters as p

__author__ = "Fabian Neumann (KIT)"
__copyright__ = f"Copyright 2020 {__author__}, GNU GPL 3"


def build_neural_network(train_set, params):
    samples = h.multiindex2array(train_set.index)
    model = ann.MLPRegressor(**params)
    model.fit(samples.T, train_set)
    return model


if __name__ == "__main__":

    cf = snakemake.config

    dataset = h.load_data(datafile)
    distribution = h.NamedJ(cf["uncertainties"])

    train_set, test_set = train_test_split(dataset, **cf["train_test_split"])

    # Model

    neural_network = build_neural_network()

    pickle.dump(neural_network, snakemake.output.ann)

    # Evaluation

    train_samples = h.multiindex2df(train_set.index)
    train_predictions = h.build_ann_prediction(neural_network, train_samples, train_set)

    test_samples = h.multiindex2df(test_set.index)
    test_predictions = h.build_ann_prediction(model, test_samples, test_set)

    p.plot_histograms(
        dataset, [train_predictions, test_predictions], fn=snakemake.output.plot
    )

    h.calculate_errors(train_predictions, train_set).to_csv(
        snakemake.output.train_errors, **cf["csvargs"]
    )
    h.calculate_errors(test_predictions, test_set).to_csv(
        snakemake.output.train_errors, **cf["csvargs"]
    )

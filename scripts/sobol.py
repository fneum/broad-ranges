"""[summary]"""

import chaospy
import numpoly
import pandas as pd

import _helpers as h
import _plotters as p

__author__ = "Fabian Neumann (KIT)"
__copyright__ = f"Copyright 2020 {__author__}, GNU GPL 3"


def calculate_sobol(surrogate, distribution, sobol="t", decimals=3):
    func = getattr(chaospy, f"Sens_{sobol}")
    sobol = func(surrogate.poly, distribution.J).round(decimals)
    return pd.DataFrame(sobol, index=distribution.names, columns=surrogate.names)


def calculate_sobol_m2(surrogate, distribution, decimals=3):
    sobol = chaospy.Sens_m2(surrogate, distribution.J).round(decimals)
    return pd.DataFrame(sobol, index=distribution.names, columns=distribution.names)


if __name__ == "__main__":

    cf = snakemake.config

    distribution = h.NamedJ(cf["uncertainties"])

    surrogate = h.NamedPoly.from_txt(snakemake.input[0])

    sobol = calculate_sobol(surrogate, distribution, snakemake.wildcards.sobol)

    sobol.to_csv(snakemake.output.data, **cf["csvargs"])

    p.plot_sobol(sobol, fn=snakemake.output.plot)

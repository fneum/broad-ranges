"""[summary]"""

import pypsa
import os
from pypsa.descriptors import nominal_attrs

__author__ = "Fabian Neumann (KIT)"
__copyright__ = f"Copyright 2020 {__author__}, GNU GPL 3"

# Add pypsa-eur(-mga) scripts to path for import
sys.path.insert(0, os.getcwd() + "/subworkflows/pypsa-eur/scripts")
sys.path.insert(0, os.getcwd() + "/subworkflows/pypsa-eur-mga/scripts")

from solve_network import (
    add_battery_constraints,
    add_EQ_constraints,
    add_CCL_constraints,
)

from generate_alternative import (
    process_objective_wildcard,
    define_mga_objective,
    define_mga_constraint,
    solve_network,
    to_pattern,
)

import logging

logger = logging.getLogger(__name__)


def add_fixed_capacity_constraint(n):

    component, pattern = snakemake.wildcards.fixedcarrier.split("+")
    position = float(snakemake.wildcards.position)
    attr = nominal_attrs[component]

    n_min = pypsa.Network(snakemake.input.min_network)
    n_max = pypsa.Network(snakemake.input.max_network)

    nom_min = n_min.df(component)[attr + "_opt"].filter(regex=to_regex(pattern)).sum()
    nom_max = n_max.df(component)[attr + "_opt"].filter(regex=to_regex(pattern)).sum()

    del n_min, n_max

    rhs = nom_min + position * (nom_max - nom_min)

    lhs = linexpr((1, get_var(n, c, attr))).sum()

    define_constraints(n, lhs, "==", rhs, "GlobalConstraint", "mu_fixed")

    logger.info(f"Fixing {pattern} {component} total capacity at {rhs/1e3} GW.")


def make_mga_model(n, sns):
    """Calls extra functionality modules."""

    obj_wc = snakemake.wildcards.objective.split("+")
    epsilon = float(snakemake.wildcards.epsilon)
    opts = snakemake.wildcards.opts.split("-")

    process_objective_wildcard(n, obj_wc)
    define_mga_objective(n)
    define_mga_constraint(n, sns, epsilon=epsilon, with_fix=True)

    add_battery_constraints(n)
    if "CCL" in opts and n.generators.p_nom_extendable.any():
        add_CCL_constraints(n, config)
    for o in opts:
        if "EQ" in o:
            add_EQ_constraints(n, o)

    if "fixedcarrier" in snakemake.wildcards:
        add_fixed_capacity_constraint(n)

    logger.info("finished extra functionality")


if __name__ == "__main__":

    n = pypsa.Network(snakemake.input[0])

    solve_opts = snakemake.config["solving"]["options"]
    opts = snakemake.wildcards.opts.split("-")

    # clip marginal cost
    threshold = 0.05
    n.generators.loc[n.generators.marginal_cost <= threshold, "marginal_cost"] = 0.0
    n.storage_units.loc[
        n.storage_units.marginal_cost <= threshold, "marginal_cost"
    ] = 0.0

    n = solve_network(
        n,
        config=snakemake.config,
        opts=opts,
        extra_functionality=make_mga_model,
        skip_objective=True,
    )

    n.export_to_netcdf(snakemake.output[0])

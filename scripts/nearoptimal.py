"""[summary]"""

import pypsa
import os

__author__ = "Fabian Neumann (KIT)"
__copyright__ = f"Copyright 2020 {__author__}, GNU GPL 3"

# Add pypsa-eur(-mga) scripts to path for import
sys.path.insert(0, os.getcwd() + "/subworkflows/pypsa-eur/scripts")
sys.path.insert(0, os.getcwd() + "/subworkflows/pypsa-eur-mga/scripts")

from solve_network import add_battery_constraints

from generate_alternative import (
    process_objective_wildcard,
    define_mga_objective,
    define_mga_constraint,
    solve_network,
)


def to_mga_model(n, sns):
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


if __name__ == "__main__":

    n = pypsa.Network(snakemake.input[0])

    solve_opts = snakemake.config["solving"]["options"]
    opts = snakemake.wildcards.opts.split("-")

    n = solve_network(
        n,
        config=snakemake.config,
        opts=opts,
        extra_functionality=to_mga_model,
        skip_objective=True,
    )

    n.export_to_netcdf(snakemake.output[0])

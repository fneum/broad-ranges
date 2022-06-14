# SPDX-FileCopyrightText: : 2021-2022 Fabian Neumann, Tom Brown
#
# SPDX-License-Identifier: MIT

"""[summary]"""

__author__ = "Fabian Neumann (KIT)"
__copyright__ = f"Copyright 2020 {__author__}, GNU GPL 3"

import pypsa
import sys
import os

# Add pypsa-eur scripts to path for import
sys.path.insert(0, os.getcwd() + "/subworkflows/pypsa-eur/scripts")

from solve_network import solve_network, prepare_network


def override_network(n):

    # need unique naming between links and lines
    # when objective includes lines and links
    n.lines.index = ["LN{}".format(i) for i in n.lines.index]
    n.links.index = ["LK{}".format(i) for i in n.links.index]

    ln_config = snakemake.config["lines"]
    n.lines.s_nom_max = n.lines.apply(
        lambda line: max(
            line.s_nom + ln_config["s_nom_add"],
            line.s_nom * ln_config["s_nom_factor"],
        ),
        axis=1,
    )

    lk_config = snakemake.config["links"]
    n.links.p_nom_max = float(lk_config["p_nom_max"])


if __name__ == "__main__":

    n = pypsa.Network(snakemake.input[0])

    override_network(n)

    solve_opts = snakemake.config["solving"]["options"]
    opts = snakemake.wildcards.opts.split("-")

    n = prepare_network(n, solve_opts=solve_opts)
    n = solve_network(n, config=snakemake.config, opts=opts)

    n.export_to_netcdf(snakemake.output[0])

import pypsa
import sys
import os

# Add pypsa-eur scripts to path for import
sys.path.insert(0, os.getcwd() + "/pypsa-eur/scripts")

from solve_network import solve_network, prepare_network


if __name__ == "__main__":

    n = pypsa.Network(snakemake.input[0])

    n = prepare_network(n)
    n = solve_network(n)

    n.export_to_netcdf(snakemake.output[0])

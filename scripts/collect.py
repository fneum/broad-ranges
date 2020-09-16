"""[summary]"""

import multiprocessing as mp
import pandas as pd
import pypsa

__author__ = "Fabian Neumann (KIT)"
__copyright__ = f"Copyright 2020 {__author__}, GNU GPL 3"


def assign_carriers(n):

    if "Load" in n.carriers.index:
        n.carriers = n.carriers.drop("Load")

    if "carrier" not in n.lines:
        n.lines["carrier"] = "AC"

    if n.links.empty:
        n.links["carrier"] = pd.Series(dtype=str)

    config = {
        "AC": {"color": "rosybrown", "nice_name": "HVAC Line"},
        "DC": {"color": "darkseagreen", "nice_name": "HVDC Link"},
    }
    for c in ["AC", "DC"]:
        if c in n.carriers.index:
            continue
        n.carriers = n.carriers.append(pd.Series(config[c], name=c))


def aggregate_costs(n, by_carrier=True):

    assign_carriers(n)

    components = dict(
        Link=("p_nom", "p0"),
        Generator=("p_nom", "p"),
        StorageUnit=("p_nom", "p"),
        Store=("e_nom", "p"),
        Line=("s_nom", None),
    )

    costs = {}
    for c in n.iterate_components(components.keys()):
        p_nom, p_attr = components[c.name]
        if c.df.empty:
            continue
        p_nom += "_opt"
        costs[(c.list_name, "capital")] = (
            (c.df[p_nom] * c.df.capital_cost).groupby(c.df.carrier).sum()
        )
        if p_attr is not None:
            p = c.pnl[p_attr].multiply(n.snapshot_weightings, axis=0).sum()
            if c.name == "StorageUnit":
                p = p.loc[p > 0]
            costs[(c.list_name, "marginal")] = (
                (p * c.df.marginal_cost).groupby(c.df.carrier).sum()
            )
    costs = pd.concat(costs) / 1e9  # bn EUR/a

    if by_carrier:
        costs = costs.groupby(level=2).sum()

    return costs


def retrieve_data(fn):

    n = pypsa.Network(fn)

    stats = pd.concat(
        [
            n.generators.groupby("carrier").p_nom_opt.sum() / 1e3,  # GW
            n.storage_units.groupby("carrier").p_nom_opt.sum() / 1e3,  # GW
        ]
    )

    lines = n.lines.eval("length * s_nom_opt / 1e6").sum()
    links = n.links.eval("length * p_nom_opt / 1e6").sum()
    stats["transmission"] = lines + links

    stats["offwind"] = stats["offwind-ac"] + stats["offwind-dc"]
    stats["wind"] = stats["onwind"] + stats["offwind"]

    stats["tsc"] = aggregate_costs(n).sum()

    to_drop = ["ror", "hydro", "PHS"]
    stats.drop(to_drop, inplace=True)

    stats.name = fn

    del n

    return stats


def parse2multiindex(df):
    def parse(fn):
        data = {}
        for o in fn[:-3].split("-"):
            s = o.split("+")
            if len(s) > 1:
                carrier = s[0]
                value = float(s[1])
                data[carrier] = value
        return pd.Series(data)

    df.columns = pd.MultiIndex.from_frame(pd.concat(map(parse, df.columns), axis=1).T)

    return df


if __name__ == "__main__":

    nprocesses = mp.cpu_count()

    with mp.Pool(processes=nprocesses) as pool:
        data = pool.map(retrieve_data, snakemake.input[0])

    df = parse2multiindex(pd.concat(data, axis=1))

    df.to_csv(snakemake.output.data)

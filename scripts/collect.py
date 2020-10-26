"""[summary]"""

import multiprocessing as mp
import pandas as pd
import pypsa
import re

__author__ = "Fabian Neumann (KIT)"
__copyright__ = f"Copyright 2020 {__author__}, GNU GPL 3"


def cumulative_share(n, by="bus"):
    # adapted from pypsa-eur-mga/scripts/extract_results.py

    n.loads["load"] = n.loads_t.p.multiply(n.snapshot_weightings, axis=0).sum()
    n.generators["energy"] = n.snapshot_weightings @ n.generators_t.p
    n.storage_units["energy"] = n.snapshot_weightings @ n.storage_units_t.inflow

    if by == "country":

        def bus2country(x):
            return x.bus[:2]

        n.loads["country"] = n.loads.apply(bus2country, axis=1)
        n.generators["country"] = n.generators.apply(bus2country, axis=1)
        n.storage_units["country"] = n.storage_units.apply(bus2country, axis=1)

    energy = (
        n.generators.groupby(by).energy.sum() + n.storage_units.groupby(by).energy.sum()
    )
    load = n.loads.groupby(by).load.sum()

    df = pd.concat([(energy / energy.sum()), (load / load.sum())], axis=1)
    df.sort_values(by="energy", inplace=True)

    return df.cumsum()


def get_gini(n):
    # adapted from pypsa-eur-mga/scripts/extract_results.py
    dfc = cumulative_share(n)
    return 1 - dfc.load.diff().multiply(dfc.energy.rolling(2).sum(), axis=0).sum()


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

    print(fn)

    n = pypsa.Network(fn)

    hvdc_b = n.links.carrier == "DC"

    stats = pd.concat(
        [
            n.generators.groupby("carrier").p_nom_opt.sum() / 1e3,  # GW
            n.storage_units.groupby("carrier").p_nom_opt.sum() / 1e3,  # GW
            n.links.loc[~hvdc_b].groupby("carrier").p_nom_opt.sum() / 1e3,  # GW
            n.stores.groupby("carrier").e_nom_opt.sum() / 1e3,  # GWh
        ]
    )

    hvac = n.lines.eval("length * s_nom_opt / 1e6").sum()
    hvdc = n.links.loc[hvdc_b].eval("length * p_nom_opt / 1e6").sum()
    stats["transmission"] = hvac + hvdc

    stats["offwind"] = stats["offwind-ac"] + stats["offwind-dc"]
    stats["wind"] = stats["onwind"] + stats["offwind"]

    stats["tsc"] = aggregate_costs(n).sum()

    stats["gini"] = get_gini(n)

    to_drop = ["ror", "hydro", "PHS", "offwind-ac", "offwind-dc"]
    stats.drop(to_drop, inplace=True)

    stats.name = fn

    del n

    return stats


def parse2multiindex(df):
    def parse(fn):
        data = {}
        fn_split = fn[:-3].split("_E")

        for o in fn_split[0].split("-"):
            s = o.split("+")
            if len(s) > 1:
                carrier = s[0] + "-cost"
                value = float(s[1])
                data[carrier] = value

        if len(fn_split) == 1:
            data["sense"] = "min"
            data["objective"] = "cost"
            data["epsilon"] = 0.0

        elif len(fn_split) == 2:
            eps_raw, obj_raw = fn_split[1].split("_")
            o = obj_raw.split("+")
            data["sense"] = o[2]
            data["objective"] = o[1] if o[1] else o[0][1:].lower()
            data["epsilon"] = float(eps_raw)

        else:
            raise NotImplementedError("Invalid filename.")

        match = re.findall(r"elec(\d+)", fn_split[0])
        if match:
            data["year"] = int(match[0])

        return pd.Series(data)

    df.columns = pd.MultiIndex.from_frame(pd.concat(map(parse, df.columns), axis=1).T)

    return df


if __name__ == "__main__":

    nprocesses = mp.cpu_count()

    with mp.Pool(processes=nprocesses) as pool:
        data = pool.map(retrieve_data, snakemake.input)

    df = parse2multiindex(pd.concat(data, axis=1))

    df.to_csv(snakemake.output.data)

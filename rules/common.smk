"""Helper functions that are called from within the workflow definitions."""

import sys, os

sys.path.insert(0, os.getcwd() + "/scripts")

import _helpers as h


def memory(w):
    """
    Calculate approximate memory requirements of solve_network rule.
    Adapted from https://github.com/PyPSA/pypsa-eur/blob/master/Snakefile.
    """
    factor = 3.0
    for o in w.opts.split("-"):
        m = re.match(r"^(\d+)h$", o, re.IGNORECASE)
        if m is not None:
            factor /= int(m.group(1))
            break
    request = int(0.35 * factor * (10000 + 195 * int(w.clusters)))
    minimum = 4000
    return max(minimum, request)


def experimental_design(w):

    scenarios = config["scenarios"][w.fidelity]
    scenarios.update(config["nearoptimal"])

    exp_design = config["experimental_design"][w.fidelity]
    fmt = exp_design.pop("round")

    uncertainty = h.NamedJ(config["uncertainties"])
    samples = uncertainty.J.sample(**exp_design).round(fmt)

    opts = scenarios["opts"]
    new_opts = []
    for s in samples.T:
        cost_set = dict(zip(uncertainty.names, s))
        cost_opts = "-".join([f"{c}+c{v}" for c, v in cost_set.items()])
        new_opts.append(f"{opts}-{cost_opts}")

    scenarios["opts"] = new_opts

    filename = "results/networks/elec_s_{clusters}_ec_lcopt_{opts}.nc"
    inputs = expand(filename, **scenarios)

    if scenarios["epsilon"]:
        filename = "results/networks/nearoptimal/elec_s_{clusters}_ec_lcopt_{opts}_E{epsilon}_O{objective}.nc"
        inputs += expand(filename, **scenarios)

    if "dependencies" in config.keys():
        dependencies = config["dependencies"]
        filename = "results/networks/dependencies/elec_s_{clusters}_ec_lcopt_{opts}_E{epsilon}_O{objective}_F{fixedcarrier}_P{position}.nc"
        partial = expand(
            filename,
            clusters=scenarios["clusters"],
            opts=scenarios["opts"],
            epsilon=dependencies["epsilon"],
            position=dependencies["position"],
            allow_missing=True,
        )
        for objective, fixedcarrier in dependencies["combinations"]:
            inputs += [
                fn.format(objective=objective, fixedcarrier=fixedcarrier)
                for fn in partial
            ]

    return inputs


def dependency_surrogates(w):
    if not "dependencies" in config.keys():
        return []

    dependencies = config["dependencies"]
    filename = (
        "results/pce/polynomial-high-{sense}-{dimension}-{epsilon}{fixed}{position}.txt"
    )
    partial = expand(
        filename,
        epsilon=dependencies["epsilon"],
        position=dependencies["position"],
        sense=["min", "max"],
        allow_missing=True
    )
    inputs = []
    for objective, fixedcarrier in dependencies["combinations"]:
        swept = objective.split("+")[1]
        fixed = fixedcarrier.split("+")[1]
        inputs += [fn.format(dimension=swept, fixed=fixed) for fn in partial]

    return inputs

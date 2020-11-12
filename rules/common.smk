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

    scenarios = config["scenarios"]

    uncertainty = h.NamedJ(config["uncertainties"])

    exp_design = config["experimental_design"]
    fmt = exp_design.pop("round")
    samples = uncertainty.J.sample(**exp_design).round(fmt)

    new_opts = []
    for opts in scenarios["opts"]:
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

    return inputs


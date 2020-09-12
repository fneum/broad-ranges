# adapted from https://github.com/PyPSA/pypsa-eur/blob/master/Snakefile
def memory(w):
    factor = 3.0
    for o in w.opts.split("-"):
        m = re.match(r"^(\d+)h$", o, re.IGNORECASE)
        if m is not None:
            factor /= int(m.group(1))
            break
    return int(factor * (10000 + 195 * int(w.clusters)))


def samples(w):
    return []

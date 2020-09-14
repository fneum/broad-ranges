configfile: "config.yaml"

subworkflow pypsaeur:
    workdir: "pypsa-eur"
    configfile: "config.pypsaeur.yaml"

include: "rules/common.smk"


rule solve_network:
    input: pypsaeur("networks/elec_s_{clusters}_ec_lcopt_{opts}.nc")
    output: "results/networks/elec_s_{clusters}_ec_lcopt_{opts}.nc"
    threads: 4
    resources: mem=memory
    script: "scripts/solve.py"


rule collect_samples:
    input: samples
    output: "results/training_set.csv"

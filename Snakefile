
# https://unix.stackexchange.com/questions/45583/argument-list-too-long-how-do-i-deal-with-it-without-changing-my-command
import os
os.system("ulimit -s 65536")

configfile: "config.yaml"

subworkflow pypsaeur:
    workdir: "subworkflows/pypsa-eur"
    configfile: "config.pypsaeur.yaml"

include: "rules/common.smk"


wildcard_constraints:
    fidelity="(high|low)",
    epsilon="[0-9\.]*",
    sobol="(t|m|m2)",
    dimension="(cost|wind|onwind|offwind|solar|storage|transmission|H2|battery)"


rule solve_network:
    input: pypsaeur("networks/elec_s_{clusters}_ec_lcopt_{opts}.nc")
    output: "results/networks/elec_s_{clusters}_ec_lcopt_{opts}.nc"
    threads: 4
    resources: mem=memory
    script: "scripts/solve.py"


rule solve_nearoptimal_network:
    input: rules.solve_network.output[0]
    output: "results/networks/nearoptimal/elec_s_{clusters}_ec_lcopt_{opts}_E{epsilon}_O{objective}.nc"
    threads: 4
    resources: mem=memory
    script: "scripts/nearoptimal.py"


if config["enable"]["collect_samples"]:
    rule collect_samples:
        input: experimental_design
        output: "results/dataset_{fidelity}.csv"
        threads: 24
        resources: mem=30000
        script: "scripts/collect.py"


# TODO placeholder for analysis
# rule analyse_full:
#     input: "results/dataset.csv"
#     output: []
#     threads: 1
#     resources: mem=8000
#     script: "scripts/analyse_full.py"


rule build_surrogate_model:
    input:
        **{fidelity: "results/dataset_" + fidelity +".csv" for fidelity in config["scenarios"].keys()}
    output:
        low_polynomial="results/pce/polynomial-low-{sense}-{dimension}-{epsilon}.txt",
        high_polynomial="results/pce/polynomial-high-{sense}-{dimension}-{epsilon}.txt",
        train_errors="results/pce/train-errors-{sense}-{dimension}-{epsilon}.csv",
        test_errors="results/pce/test-errors-{sense}-{dimension}-{epsilon}.csv",
        plot="results/pce/histogram-{sense}-{dimension}-{epsilon}.pdf"
    threads: 1
    resources: mem=8000
    script: "scripts/surrogate.py"


rule build_all_surrogates:
    input:
        epsilon=expand("results/pce/polynomial-high-{sense}-{dimension}-{epsilon}.txt",
               dimension=["onwind", "offwind", "wind", "solar", "transmission", "H2", "battery"],
               sense=["min", "max"],
               epsilon=config["nearoptimal"]["epsilon"]),
        cost="results/pce/polynomial-high-min-cost-0.0.txt"


# ruleorder: calculate_sensitivity_indices > calculate_nearoptimal_sensitivity_indices


rule calculate_sensitivity_indices:
    input: "results/pce/polynomial-high-min-cost.txt"
    output:
        data="results/graphics/sobol-min-cost.csv",
        plot="results/graphics/sobol-min-cost.pdf"
    threads: 1
    resources: mem=8000
    script: "scripts/sobol.py"


# rule calculate_nearoptimal_sensitivity_indices:
#     input: "results/pce/polynomial-high-{sense}-{dimension}.txt"
#     output:
#         data="results/graphics/sobol-{sense}-{dimension}.csv",
#         plot="results/graphics/sobol-{sense}-{dimension}.pdf"
#     threads: 1
#     resources: mem=8000
#     script: "scripts/nearoptimal_sobol.py"

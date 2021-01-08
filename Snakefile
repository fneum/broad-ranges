
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
    order="[0-9]*",
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
        polynomial="results/pce/polynomial-{order}-{sense}-{dimension}.txt",
        train_errors="results/pce/train-errors-{order}-{sense}-{dimension}.csv",
        test_errors="results/pce/test-errors-{order}-{sense}-{dimension}.csv",
        plot="results/pce/histogram-{order}-{sense}-{dimension}.pdf"
    threads: 1
    resources: mem=16000
    script: "scripts/surrogate.py"


ruleorder: calculate_sensitivity_indices > calculate_nearoptimal_sensitivity_indices


rule calculate_sensitivity_indices:
    input: "results/pce/polynomial-{order}-min-cost.txt"
    output:
        data="results/graphics/sobol-{order}-min-cost-{sobol}.csv",
        plot="results/graphics/sobol-{order}-min-cost-{sobol}.pdf"
    threads: 1
    resources: mem=8000
    script: "scripts/sobol.py"


rule calculate_nearoptimal_sensitivity_indices:
    input: "results/pce/polynomial-{order}-{sense}-{dimension}.txt"
    output:
        data="results/graphics/sobol-{order}-{sense}-{dimension}-{sobol}.csv",
        plot="results/graphics/sobol-{order}-{sense}-{dimension}-{sobol}.pdf"
    threads: 1
    resources: mem=8000
    script: "scripts/nearoptimal_sobol.py"

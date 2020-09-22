configfile: "config.yaml"

subworkflow pypsaeur:
    workdir: "subworkflows/pypsa-eur"
    configfile: "config.pypsaeur.yaml"

include: "rules/common.smk"


wildcard_constraints:
    order="[0-9]*",
    sobol="(t|m|m2)"


rule solve_network:
    input: pypsaeur("networks/elec_s_{clusters}_ec_lcopt_{opts}.nc")
    output: "results/networks/elec_s_{clusters}_ec_lcopt_{opts}.nc"
    threads: 4
    resources: mem=memory
    script: "scripts/solve.py"


if config["enable"]["collect_samples"]:
    rule collect_samples:
        input: experimental_design
        output:
            data="results/dataset.csv"
        threads: 1
        resources: mem=8000
        script: "scripts/collect.py"


rule analyse_full:
    input: "results/dataset.csv"
    output: [] # TODO
    threads: 1
    resources: mem=8000
    script: "scripts/analyse_full.py"


rule build_surrogate_model:
    input: "results/dataset.csv"
    output:
        polynomial="results/pce/polynomial-{order}.txt",
        train_errors="results/pce/train-errors-{order}.csv",
        test_errors="results/pce/test-errors-{order}.csv",
        plot="results/pce/histogram-{order}.pdf"
    threads: 1
    resources: mem=16000
    script: "scripts/surrogate.py"


rule calculate_sensitivity_indices:
    input: rules.build_surrogate_model.output.polynomial
    output:
        data="results/graphics/sobol-{order}-{sobol}.csv",
        plot="results/graphics/sobol-{order}-{sobol}.pdf"
    threads: 1
    resources: mem=8000
    script: "scripts/sobol.py"


rule build_neural_network:
    input: "results/dataset.csv"
    output:
        ann="results/ann/neural_network.pickle",
        train_errors="results/ann/train-errors.csv",
        test_errors="results/ann/test-errors.csv",
        plot="results/ann/histogram.pdf"
    threads: 1
    resources: mem=8000
    script: "scripts/neural_network.py"


cp config.technology-data.yaml subworkflows/technology-data/config.yaml
cd subworkflows/technology-data
rm -rf outputs
snakemake -j 1 -f compile_cost_assumptions
cp outputs/costs_2050.csv ../pypsa-eur/resources/costs.csv
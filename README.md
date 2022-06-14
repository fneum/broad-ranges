# Broad Ranges of Investment Configurations for Renewable Power Systems, Robust to Cost Uncertainty and Near-Optimality

[![REUSE status](https://api.reuse.software/badge/github.com/fneum/broad-ranges)](https://api.reuse.software/info/github.com/fneum/broad-ranges)
![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/fneum/broad-ranges?include_prereleases)
![Size](https://img.shields.io/github/repo-size/fneum/broad-ranges)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.6642651.svg)](https://doi.org/10.5281/zenodo.6642651)
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

This repository contains the entire scientific project, including code and manuscript.
## Abstract

To achieve ambitious CO$_2$ emission reduction targets quickly, the planning of
energy systems must accommodate societal preferences, e.g.~regarding
transmission reinforcements or onshore wind parks, and must also acknowledge
uncertainties of technology cost projections. To date, however, many models lean
towards only minimising system cost and only using a single set of cost
projections. Here, we address both criticisms in unison. While taking account of
cost uncertainties, we apply multi-objective optimisation techniques to explore
trade-offs in a fully renewable European electricity system between rising
system cost and the deployment of individual technologies for generating,
storing and transporting electricity. We identify boundary conditions that must
be met for cost-efficiency that are robust to how cost developments will unfold;
for instance, we find that some grid reinforcement and long-term storage
alongside large wind capacities appear essential. We reveal that near the
cost-optimum a broad spectrum of technologically diverse options exist, which
allows policymakers to make trade-offs regarding unpopular infrastructure.

## Repository Structure

- `data` contains additional input data sources
- `notebooks` contains results analysis and plotting scripts, as well as created graphics
- `scripts` contains Python scripts of the core workflow
- `results` contains workflow outputs
- `paper` contains the adjunct manuscript
- `subworkflows` contains links to dependent workflows including PyPSA-Eur and PyPSA-Eur-MGA
## Installation

```sh
mamba env create -f environment.yaml
```
## Usage

Set configuration in `config.yaml` and `config.pypsaeur.yaml`.

Run

```sh
snakemake -call build_all_surrogates
```

## License

Different licenses apply to different parts of the repository. See [specifications here](.reuse/dep5).
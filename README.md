# Surrogate Modelling for PyPSA-Eur

# Parametric Uncertainty

## Motivation

In optimisation theory the KKT equations can tell us about the sensitivity of the optimal point to small changes in the constraints. However, KKT does not tell us anything about the global behaviour of the objective function on the feasible space, which can be highly non-linear.

Flat directions near the optimum are very important for policymakers, because they allow them to make decisions based on other, non-economic criteria, without affecting the cost-effectiveness of the system.

First approach could be parameter sweeps in the space of costs and constraint constants (Schlachtberger, 2018)

Few parameters are really relevant, what determines the cost structure (Moret, 2017)

Focus on cost inputs and constraints affecting outputs relevant to public acceptance.

Explorative analysis.

Guiding robust and comprehensive policy.

What makes good policy advice?

Where in parameter space is it cheap?

- What's the danger zone?
- Or are 90% of the uncertainty space within 10% of the least cost solutions?

Whats the relation between GSA and stochastic and robust optimisation?

Look at measure of robustness: given a cost-optimal solution for one parameter set, how expensive would it be for other cost realisations?

non-differentiality, non-continuity, penny switching (in different regions), lifting the degeneracy, flipping points

Relation between complexity and uncertainty. Structural uncertainty vs parametric uncertainty.

complex models are often no better than simple models, when it comes to prediction

decision making under deep uncertainty: what are available strategies? what are vulnerabilities? what strategy reduces vulnerabilities.

explore more, communicate better

Models can (from Tom's lecture):

- easily give a false sense of exactness
- under- or overestimate rates of change
- underestimate social factors (concerns about nuclear, transmission, wind)
- extrapolate based on uncertain data (learning curves PV)
- focus on easy-to-solve rather than policy-relevant problems
- neglect uncertainty (short: weather forecasts, long: technology cost)
- neglect need for robustness (contingencies, attack)
- neglect complex interactions of markets and incentive structures (market power, lumpiness)
- neglect non-linearities and non-convexities (power flow, learning curves, behaviour)

## Uncertainty

Difference between parametric and structural uncertainty

- parametric: technology cost, useable potentials, discount rate, weather year
- structural: power flow approximation, spatial and temporal aggregation, technology simplification (e.g. UC)

Origin of cost uncertainty:

- technological learning and rate of deployment are unknown
- would be useful to identify regions in the parameter space which are very cheap
- can aim for them with policy; "controllable uncertainty"

Better to exclude discount rate in a capital intensive system (just a pre-factor to the objective).

There are many hidden assumptions!

Sources of uncertainty (a la Chris Dent):

- parametric (cost)
- condition (boundary, initial)
- functional (full feasible space not known)
- stochastic (randomness)
- solution (solver tolerance)
- structural (power flow)
- multi-model (across models)
- decision (what can be influenced)

## Cost Data

What are the distributions of technology cost projections?

- uniform:
- normal: (Mavromatidis, 2018)
- check sensitivity towards different distributions; e.g. normal vs uniform

Locate the status of learning curves for different technologies (fit learning parameter)

- upper limit: capacity to supply the whole world
- lower limit: today's capacity

Main challenge is quantifying the input uncertainties! (Moret, 2017)

take cost assumptions from https://github.com/PyPSA/pypsa-eur/pull/184

Tröndle paper

Danish Energy Agency

Use uniform distribution due to lack of better data (Moret, 2017)

Costs all the way to zero? E.g. what happens with very cheap solar?

Check how sensitive uncertainty propagation is wrt to the choice of distribution (e.g. normal vs uniform). What changes?

## Sensitivity Analysis

Difference between uncertainty and sensitivity:

- uncertainty analysis: to what extent does uncertainty exist in the outputs
- sensitivity analysis: what input parameters influence the outputs the most

scenario-based sensitivity analysis

- crude, local, subjective (Usher, 2015)
- more qualified and contextual description possible
- storyline with narrative elements (DeCarolis, 2017)

local sensitivity analysis

- parameter sweeps one dimension at a time (Schlachtberger, 2018)
- multi-dimensional sensitivity analysis analyse more dimensions at once
- main interdependencies may already be covered by single-parameter variation (Schlachtberger, 2018)
- "desirable to develop analytic results from simplified models that maintain the original model's sensitivities" (Schachtberger, 2018)
- use partial derivatives
- drawback: only small fraction of uncertainty space is sampled

global sensitivity analysis

- changing multiple paramters at a time; co-varying inputs; multi-parameter variations
- to identify interesting scenarios (Usher, 2015)
- "quantifying the effects of random input variables onto the variance of the response of a mathematical model" (Sudret, 2008)

a good overview: https://uncertainpy.readthedocs.io/en/latest/theory.html

## Local Sensitivity Analysis

- repeat cost-optimal scenarios for thesis
- argument: more nodes, 100% renewable, leads up to global sensitivity analysis
- redo with PyPSA-Eur-Sec

## Experimental Design

low-discrepancy MC sampling

- https://chaospy.readthedocs.io/en/master/sampling/sequences.html
- https://en.wikipedia.org/wiki/Low-discrepancy_sequence
- https://en.wikipedia.org/wiki/Quasi-Monte_Carlo_method

using halton sequence

## Surrogate Modelling with Polynomical Chaos Expansion

premise: outcome of original model cannot be obtained easily (e.g. computational constraints)

synonyms for surrogate: emulators, approximation models, response surface methods, metamodels

use cases: design space exploration, sensitivity analysis, what-if analyis

consider only simplified/aggregated outputs, junk complicated models

only input/output behaviour is important (link to machine learning)

scaling

- number of inputs: adds more dimensions to uncertainty space, curse of dimensionality
- number of outputs: should scale well as each output has its own polynomial, are independent

Use PCE to build surrogate models for calculating Sobol sensitivity indices analytically as a post-processing (Sudret, 2008)

computational cost of sensitivity indices reduces to that of estimating PCE coefficients, 2-3 orders of magnitude faster than traditional MC evaluation (Sudret, 2008)

The polynomials can be used for analytic calculations (like Sobol indices) and thereby have added value compared to pure machine learning approaches

polynomial regression = linear regression

- https://scikit-learn.org/stable/modules/linear_model.html#polynomial-regression-extending-linear-models-with-basis-functions

implement multifidelity approach

- Tröndle paper: high fidelity: 10 samples, 400 nodes, 4-hourly; low fidelity: 150 samples, 25 nodes, 4-hourly; no DC power flow

add transmission volume limit to input parameters

Are LHS/Halton/Sobol sampling better than gridded sampling?

Validation

- benchmark: model error <5% (Tröndle, 2020)

The mixture of distributions is not a distribution itself: https://en.wikipedia.org/wiki/Mixture_distribution

Talk to People

- Tim Tröndle, https://github.com/timtroendle/geographic-scale/
- Stefano Marelli, UQLab
- Till, PCE

## Surrogate Modelling with Neural Networks

as a simple benchmark case

Kaleb ML hints

- can build neural network separately for each output variable
- number of nodes = 2 * number of inputs
- normally 2-3 hidden layers
- 1000 samples could work, 100 could be too low

## Reconstruction

When using the surrogate model to determine the optimised system capacities
from a particular set of cost assumptions, we do not get the full original model outputs.

If we are interested in those, we can reconstruct these selectively
by adding the aggregate outcomes as `extra_functionality` constraints (e.g. onwind p_nom_opt = XXX GW) and run the full optimisation.

Caveat: due to inaccuracies of the surrogate model, one may need load shedding or other forms of slack.

## Conclusions

one of the best ways to make system cheaper is just to make offshore wind cheaper, which has high acceptance and big cost reduction potentials AND it has one of the biggest impacts on TSC.

by reducing capital cost of wind, not only expected costs decrease, but also the uncertainty band

# Near-Optimal

## Regional Analysis

more detailed regional resolution for technology aggregation (e.g. min/max solar in Southern Europe)

## More Search Directions

add additional search directions based on face normal of convex hull (Pedersen, 2020):

- https://github.com/TimToernes/MGA-PyPSA
- https://github.com/tulip-control/polytope


do C-Plots in 2D (or even 3D) with isolines corresponding to $\epsilon$ based on the convex hulls for all combinations of investment groups (do multiple dimensions at once)

## Sector Coupling

redo with sector-coupling version

## Twitter Comments

- "if we accept very slightly higher costs - insignificantly so in the big scheme of things - the flexibility to build an energy network that is politically acceptable increases substantially" @jmkorhon_en

- "Spanning a space of options for society and politics to work in instead of single least-cost-model-oitcomes. Should be a benchmark method for every study in my opinion." @alxrdk

- "bridging optimization modeling and actual, messy real-world decision making" @EmilDimanchev

- "running an optimizer to get the One True Solution can be really misleading" @ThatcherUlrich

- "Getting it done is more important than getting it right!" @BurkeKB

- "Sometimes people write a paper that is just begging to be written." @KenCaldeira

- "It excites and scares me because it’s a double edge sword. This is because something has to replace something. Usually the two are not mutually exclusive." @DrChrisClack

- "A whole new diverse world exists beyond the cost-optimal solution!" @neha__25

- "it dispells the myth of a unique path to cost-efficiency, legitimising attending other decision criteria *also* from cost perspective" @anunezjimenez

- "The idea of "optimality" is a seductive one, but in the real world (and especially in the real world of climate mitigation) almost NOTHING is optimal. It's very important to explore sensitivities to modest shifts in key drivers as this article does." @jgkoomey

- "the idea of an "optimal" energy system that "can be found" through modelling both appeals and repulses me. It's a major reason I remain sceptical about usefulness of optimization. This kind of papers, IMHO, boost these models policy relevance." @anunezjimenez

# Near-Optimal and Parametric Uncertainty

build surrogate including $\epsilon$ and sense as uncertain parameter

- higher sampling rate for low $\epsilon$ (0.5,1,2,4,8,12,16)

heatmap of tsc in c-plots? fix capacity of n-1 components, sweep the other?

interpretation of marrying near-optimal analysis with parametric uncertainty (what's the relation?)

- tilt of (1-$\epsilon$) constraint by varying parameters changes near-optimal feasible space
- could also represent different objectives rather than costs, e.g. area used etc. to represent social aspects


# Software

Chaospy

- https://github.com/jonathf/chaospy
- https://chaospy.readthedocs.io/en/master/

Uncertainpy

- https://uncertainpy.readthedocs.io/en/latest/
- https://github.com/simetenn/uncertainpy

Polychaos

- https://github.com/timueh/PolyChaos.jl
- https://timueh.github.io/PolyChaos.jl/dev/

SAlib

- https://github.com/SALib/SALib
- https://salib.readthedocs.io/en/latest/index.html

SMT

- https://github.com/SMTorg/smt
- https://smt.readthedocs.io/en/latest/
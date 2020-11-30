# Surrogate Modelling for PyPSA-Eur

# Parametric Uncertainty

## Uncertainty

Sources of uncertainty (a la Chris Dent):

- parametric (cost)
- condition (boundary, initial)
- functional (full feasible space not known)
- stochastic (randomness)
- solution (solver tolerance)
- structural (power flow)
- multi-model (across models)
- decision (what can be influenced)

Dimensions of uncertainty (Pye, 2018)

- technical (inexactness)
- methodological (unreliability)
- epistemological (ignorance)
- societal (social robustness)

Different time domains of uncertainty:

- short-term: operation
- long-term: planning

Difference between epistemic and aleatory uncertainty (Pfenninger, 2014)

- epistemic: more data/better models will reduce uncertainty
- aleatory:  more data/better models will **not** reduce uncertainty

Difference between parametric and structural uncertainty (extended from (DeCarolis, 2017))

- parametric: technology cost, useable potentials, discount rate, weather year, climate change, electricity demand, (behavioural/social)
- structural: power flow approximation, spatial and temporal aggregation, technology simplification (e.g. UC)

Origin of cost uncertainty:

- technological learning and rate of deployment are unknown
- would be useful to identify regions in the parameter space which are very cheap
- can aim for them with policy; "controllable uncertainty"

Better to exclude discount rate in a capital intensive system (just a pre-factor to the objective).

There are many hidden assumptions!

"Although the task of uncertainty characterization can itself
be seen as uncertain, exploring the nature of uncertainty is significantly more valuable than using deterministic, best-guess values." (Fraiture, 2020)

"the uncertainties inherent in the model structures and input parameters are at best underplayed and at worst ignored" (Yue, 2018)

## Cost Data

Main challenge is quantifying the input uncertainties! (Moret, 2017)

Two main sources of uncertainty wrt technological learning (Tröndle, 2020)

- future deployment rates unknown
- learning rate unknown

### Via Learning Curves

Locate the status of learning curves for different technologies (fit learning parameter)

- upper limit: capacity to supply the whole world
- lower limit: today's capacity

models "highly sensitive to uncertainty in the learning rates [...] due to the exponential relationship" (Mattssen, 2019)

To derive learning rates:

- IRENA renewable cost database
- www.energystorage.ninja

### Maximum-Entropy Approach

used in (Tröndle, 2020)

https://de.wikipedia.org/wiki/Maximum-Entropie-Methode

https://en.wikipedia.org/wiki/Principle_of_maximum_entropy

from Bayesian statistics

"assign an a-priori probability despite insufficient problem-specific information" (wiki)

take samples; among the set of all trial probability distributions take the one with maximal information entropy (wiki)

"the one that makes fewest assumptions about the true distribution of data" (wiki)

"entropy maximization with no testable information respects the universal constraint that the sum of the probabilities is one. Under this constraint, the maximum entropy probability distribution is the uniform distribution" (wiki)

### Distributions

What are the distributions of technology cost projections?

- uniform: (Moret, 2017) (Moret, 2016) (Shirizadeh, 2019) (Tröndle, 2020) (Pilpola, 2020) (Pizarro-Alonso, 2019) (Li, 2017) (Trutnevyte, 2013)
- normal: (Mavromatidis, 2018)
- triangle: (Li, 2020)

always independently sampled; easy but incorrect: e.g. offshore and onshore wind

The current spread of recent documented investment costs is close to uniform (Lopin, 2019)

"difficult and possibly misleading to associate a PDF to a parameter with unknown PDF" (Moret, 2016)

"The variance [..] was set according to the maturity of technologies" (Li, 2020)"

"To date, a general methodology for uncertainty characterization, assessing parameter uncertainty by type and degree [..] is missing." (Moret, 2016)

"Following a maximum entropy approach, we model [...] uncertainty with uniform distributions over ranges taken from the literature" (Tröndle, 2020)

### Ranges

"it is crucial to avoid an arbitrary a priori exclusion of parameters from the analysis" (Moret, 2017)

take pessimistic cost assumptions from `technology-data` with option `expectation: pessimistic`

- new with continuous updates since 2016
- no ranges for H2 pipeline, HVAC, HVDC, PHS, hydro, nuclear, ror, rooftop PV
- also checkout `./costcomparison.csv`

(Tröndle, 2020) supplementary material has almost always more conservative values than pessmisitic DEA database; data mostly from ETIP

| name           | min  | max  | unit      |
|----------------|------|------|-----------|
| discount       | 1.6  | 13.8 | % p.a.    |
| utility PV     | 280  | 580  | EUR/kW    |
| rooftop PV     | 760  | 1000 | EUR/kW    |
| onwind         | 800  | 1700 | EUR/kW    |
| offwind        | 1790 | 3270 | EUR/kW    |
| battery power  | 31   | 141  | EUR/kW    |
| battery energy | 36   | 166  | EUR/kWh   |
| H2 power       | 1123 | 2100 | EUR/kW    |
| H2 energy      | 6    | 12   | EUR/kWh   |
| transmission   | 700  | 1080 | EUR/MW/km |

JRC Energy Technology Reference Indicators (https://setis.ec.europa.eu/setis-output/energy-technology-reference-indicators) is relatively old and conservative from 2014

Danish Energy Agency technology cost database:

- lower and higher bounds "shall be interpreted as representing probabilities corresponding to a 90% confidence interval"

plus/minus 50% (Shirizadeh, 2019)

mostly less than plus/minus 25% (Pizarro-Alonso, 2019)

plus/minus 20% (Moret, 2017)

previous studies have considered "relatively narrow range[s] of uncertainties" (Li, 2017)

maybe: include transmission cost uncertainty by setting upper bound to 800 EUR/MW/km (which is twice our usual assumption; do this because Tröndle cost are much higher!)

maybe: add 25% to rooftop PV cost (due to missing uncertainty ranges in DEA, similar to utility PV)

Costs all the way to zero!

# Senstivity Analysis

## Motivation

In optimisation theory the KKT equations can tell us about the sensitivity of the optimal point to small changes in the constraints.
However, KKT does not tell us anything about the global behaviour of the objective function on the feasible space, which can be highly non-linear.

Flat directions near the optimum are very important for policymakers,
because they allow them to make decisions based on other,
non-economic criteria, without affecting the cost-effectiveness of the system.

First approach could be parameter sweeps in the space of costs and constraint constants (Schlachtberger, 2018)

Few parameters are really relevant, what determines the cost structure (Moret, 2017) (Usher, 2015)

Focus on cost inputs and constraints affecting outputs relevant to public acceptance.

Explorative analysis.

Guiding robust and comprehensive policy.

But what makes good policy advice?

Where in parameter space is it cheap?

- Or conversely: What's the danger zone?
- Or are 90% of the uncertainty space within 10% of the least cost solutions?
- These cheap areas should then be targeted with policy (e.g. drive technology learning by subsidies)
- smoothly leads up to future work with endogenous learning in multi-horizon investment planning!

Do I need a richer set of technologies? E.g.
solar thermal, nuclear, CAES, separate hydrogen storage to fuel cell, electrolysis, gas turbine, cavern storage, steel tank?

What is the relation between sensitivity analysis and stochastic and robust optimisation?

- sensitivity analysis: uncertainty wrapped around optimisation
- stochastic/robust optimisation: uncertainty embedded in optimisation
- Alternative: embed technological learning in optimisation (Heuberger, 2017) (Lopion, 2019), remaining uncertainty of learning rate, but "learning by doing" included

Look at measure of robustness: given a cost-optimal solution for one parameter set, how expensive would it be for other cost realisations? A form of hedging against uncertainty.

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

## Karush-Kuhn-Tucker Conditions

https://en.wikipedia.org/wiki/Karush%E2%80%93Kuhn%E2%80%93Tucker_conditions

strong duality: primal objective = dual objective

Conditions:

- stationarity
- primal feasibility
- dual feasibility
- complementary slackness

Example Applications:

- carbon cap <-> carbon price
- renewable generation target <-> required subsidy
- limited potentials <-> scarcity cost (onshore wind, transmission)

## Multi-Objective Optimisation

What we really do is calculating Pareto-fronts (e.g. when constraining line expansion) from multi-objective optimisation

https://en.wikipedia.org/wiki/Multi-objective_optimization

https://en.wikipedia.org/wiki/Pareto_efficiency

https://en.wikipedia.org/wiki/Maxima_of_a_point_set

## Categories of Sensitivity Analysis

no sensitivity analysis -> point-estimate method (Souroudi, 2013)

Difference between uncertainty and sensitivity analsis following (Usher, 2016):

- uncertainty characterization (analysis): which parameters are uncertain and to what extent
- uncertainty propagation: to what extent does uncertainty exist in the outputs
- sensitivity analysis: what input parameters influence the outputs the most

scenario-based sensitivity analysis

- crude, local, subjective (Usher, 2015)
- more qualified and contextual description possible
- storyline with narrative elements (DeCarolis, 2017)

local sensitivity analysis

- parameter sweeps one dimension at a time (Schlachtberger, 2018)
- main interdependencies may already be covered by single-parameter variation (Schlachtberger, 2018)
- "desirable to develop analytic results from simplified models that maintain the original model's sensitivities" (Schachtberger, 2018)
- use partial derivatives
- multi-dimensional sensitivity analysis analyse more dimensions at once
- drawback: only small fraction of uncertainty space is sampled
- LSA "may not identify influential parameters in planning (Pizarro-Alonso, 2019)

global sensitivity analysis

- changing multiple paramters at a time; co-varying inputs; multi-parameter variations
- to identify interesting scenarios (Usher, 2015)
- "quantifying the effects of random input variables onto the variance of the response of a mathematical model" (Sudret, 2008)
- elementary effects method (Morris screening, one-at-a-time) to counteract high computational requirements (Usher, 2015) (Pizarro-Alonso, 2019) (Moret, 2016)
- surrogate models using polynomial chaos expansion

a good overview: https://uncertainpy.readthedocs.io/en/latest/theory.html

## Local Sensitivity Analysis

a la (Schlachtberger, 2018), but with higher resolution and 100% renewable

use as build-up to global sensitivity analysis

- nodes
- temporal resolution (averaging, segmentation)
- cost
- potentials (offwind, onwind, solar)
- line expansion (lv1.0 ... lcopt)
- weather years
- equity requirements (done already!)
- emission reduction target

## Experimental Design

improvements to random sampling!

intendended to achieve "efficient coverage" (Usher, 2015)

low-discrepancy Monte-Carlo (MC) sampling

- https://chaospy.readthedocs.io/en/master/sampling/sequences.html
- https://en.wikipedia.org/wiki/Low-discrepancy_sequence
- https://en.wikipedia.org/wiki/Quasi-Monte_Carlo_method

using halton sequence

alternatives are Latin hypercube sampling (Tröndle, 2020) and Method of Morris (Usher, 2015)

Categories of sampling methods (following Palar, 2016):

- structured: orthogonal polynomial roots, Newton-Cotes, sparse grid
- unstructured: random, LHS, low-discrepancy

unstructured types typically used for uniform distributions

Discuss which are attractive and why!

## Global Sensitivity Analysis with PCE Surrogates

premise: outcome of original model cannot be obtained easily (e.g. computational constraints)

synonyms for surrogate: emulators, approximation models, response surface methods, metamodels

use cases: design space exploration, sensitivity analysis, what-if analyis

consider only simplified/aggregated outputs, junk complicated models

"reduce the number of deterministic evaluations while retaining accuracy" (Palar, 2016)

only input/output behaviour is important (link to machine learning)

scaling

- number of inputs: adds more dimensions to uncertainty space, curse of dimensionality [12 in (Tröndle, 2020), 36 in (Pilpola, 2020)]
- number of outputs: should scale well as each output has its own polynomial, are independent

### Polynomial Chaos Expansion (PCE)

Use PCE to build surrogate models for calculating Sobol sensitivity indices analytically as a post-processing (Sudret, 2008)

- Hilbert space technique
- idea: "expand random variables as a linear combination of orthogonal basis functions weighted by deterministic coefficients." (Mühlpfordt, 2019)
- "A polynomial chaos expansion for a random variable is what a classic Fourier series is for a periodic signal" (Mühlpfordt, 2019)
- "In practice the truncated PCE is used, meaning that only finitely many coefficients represent the random variable" (Mühlpford, 2019)
- steps: find polynomial basis, calculate polynomial coefficients
- techniques for finding coefficients:

  - spectral projection: project response onto basis functions using inner products and orthogonal polynomials (Palar, 2016)
  - point collocation: obtain coefficients by performing  regression (Palar, 2016)

want to do non-intrusive (i.e. wrapping around model) point collocation (i.e. MC sampling) based on regression

sparsity-of-effects principle:

- "system can be represented using a small number or statistically significant effects" (Berchier, 2016)
- in PCE: use only small number of polynomials from polynomial basis
- achieved by adding regularisation/penalty term to regression to favour sparse solutions (e.g. LARS regression)

over-sampling ratio (OSR):

- ratio between number of samples and polynomial basis (Palar, 2016)
- recommended is a value of 2 (Hoser, Walters, Balch via Palar, 2016)
- too high OSR: very coarse approximation; too low OSR: risk of over-fitting (Palar, 2016)

Relation of Polynomial Chaos and Principal Component Analysis:

- PCE is like PCA over an orthogonal vector space of polynomials

separate influential from non-influential parameters

computational cost of sensitivity indices reduces to that of estimating PCE coefficients, 2-3 orders of magnitude faster than traditional MC evaluation (Sudret, 2008)

The polynomials can be used for analytic calculations (like Sobol indices) and thereby have added value compared to pure machine learning approaches

polynomial regression = linear regression

- https://scikit-learn.org/stable/modules/linear_model.html#polynomial-regression-extending-linear-models-with-basis-functions

- Tröndle paper: high fidelity: 10 samples, 400 nodes, 4-hourly; low fidelity: 150 samples, 25 nodes, 4-hourly; no DC power flow

further uncertain input parameters:

- rountrip efficiency of hydrogen storage
- allowed transmission volume
- cost of transmission?

Validation

- benchmark: model error <5% (Tröndle, 2020)

The mixture of distributions is not a distribution itself: https://en.wikipedia.org/wiki/Mixture_distribution

Talk to People

- Tim Tröndle, https://github.com/timtroendle/geographic-scale/
- Stefano Marelli, UQLab
- Tillmann, PCE

### Multi-Fidelity Approach

following Jiant et al., 2019

only high-fidelity for surrogate model is time-consuming,
but only low-fidelity may result in distorted/inaccurate surrogate models

"If LF model is much cheaper to evaluate than the HF model, then significant
computational savings will be obtained" (Ng, 2012)

integrate both with multi-fidelity approach: correct LF model such that LF output matches HF values

build an own surrogate model for the correction function (i.e. "correction expansion") with identical polynomial basis

rule: "polynomial term in the correction expansion has to e a subset of the LF expansion" (Palar, 2016)

model fusion by joining both into a single expansion

here another advantage of low-discrepancy sampling: HF sample is subset of LF sample

"in cases with high correlation between LF an HF ($R^2 \geq 0.9$) the MF will very likely
increase the approximation accuracy relative to single HF" (Palar, 2016)

__Additive scaling approach__

Correct for the *difference* between HF samples and LF surrogate at HF sample points

$$
\tilde{f}_{MF}(\xi) = \tilde{f}_{LF} + \tilde{\alpha}(\xi)
$$

where $\tilde{\alpha}(\xi)$ is the additive scaling function. The sampled scaling factors are calculated with

$$
\alpha(\xi) = f_{HF}(\xi) - \tilde{f}_{LF}(\xi)
$$

Construct scaling function $\tilde{\alpha}(\xi)$ from scaling factors $\alpha(\xi)$ via some regression fit.

__Multiplicative scaling approach__

Correct for the *ratio* between HF samples and LF surrogate at HF sample points

$$
\tilde{f}_{MF}(\xi) = \tilde{f}_{LF} \cdot \tilde{\beta}(\xi)
$$

where $\tilde{\beta}(\xi)$ is the scaling function. The sampled scaling factors are calculated with

$$
\beta(\xi) = \frac{ f_{HF}(\xi) }{ \tilde{f}_{LF}(\xi) }
$$

Construct scaling function $\tilde{\beta}(\xi)$ from scaling factors $\beta(\xi)$ via some regression fit.

__Hybrid scaling approach__

Combine both additively via weighting factor $\omega$ as convex combination:

$$
\tilde{f}_{MF}(\xi) = \omega \left( \tilde{f}_{LF}(\xi) + \tilde{\alpha}(\xi) \right) + (1-\omega) \left( \tilde{f}_{LF}(\xi) \cdot \tilde{\beta}(\xi) \right)
$$

There is flexibility in the choice of $\omega$. Idea: compute $\omega$ "based on a regularisation of the combined correction function to prevent overfitting" (Ng, 2012); "minimise mean-square magnitude of add/mul correction

$$
\min_{\omega\in[0,1]} \lang\: \omega^2 \alpha^2(\xi),\: (1-\omega)^2 \beta^2(\xi) \:\rang
$$

leading to a balance between add/mul correction

$$
\omega = \frac{\lang \beta^2(\xi) \rang}{\lang \alpha^2(\xi) \rang + \lang \beta^2(\xi) \rang}
$$

## Surrogate Modelling with Neural Networks

as a simple benchmark case? rather not.

Kaleb ML hints

- can build neural network separately for each output variable
- number of nodes = 2 * number of inputs
- normally 2-3 hidden layers
- 1000 samples could work, 100 could be too low


## Conclusions

one of the best ways to make system cheaper is just to make offshore wind cheaper,
which has high acceptance and big cost reduction potentials AND it has one of the biggest impacts on TSC.

by reducing capital cost of wind, not only expected costs decrease, but also the uncertainty band (this is as expected)

# Near-Optimal

## Regional Analysis

more detailed regional resolution for technology aggregation (e.g. min/max solar in Southern Europe)

6 Regions variant:

- Western: GB|IE|FR|NL|BE|LU
- Central: DE|AT|CH
- Eastern: PL|CZ|HU|SK|LT|LV|EE
- Northern: NO|SE|DK|FI
- Southern: PT|ES|IT|GR|MT
- South-Eastern: RO|AL|MK|BG|RS|HR|SI|BA|ME

4 regions variant:

- Northern: NO|SE|DK|FI|LT|LV|EE
- Southern: PT|ES|IT|GR|AL|MK|BG|RS|HR|SI|BA|ME
- Western: CH|GB|IE|FR|NL|BE|LU
- Eastern: DE|AT|RO|PL|CZ|HU|SK

## More Weather Years

get a fuzzy boundary of c-plots based on repeated mga analysis for weather years 1979-2018

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

interpretation of marrying near-optimal analysis with parametric uncertainty (what's the relation?)

- tilt of (1-$\epsilon$) constraint by varying parameters changes near-optimal feasible space
- could also represent different objectives rather than costs, e.g. area used etc. to represent social aspects

"resilient SPORES" (Lombardi, 2020)

"SPORES are relatively insensitive to technology cost" (Lombardi, 2020)

Trutnevyte's EXPANSE combines MGA-type and MC-type uncertainty analysis in total 800 solutions (Li, 2017)

# Reconstruction and Disaggregation

When using the surrogate model to determine the optimised system capacities
from a particular set of cost assumptions, we do not get the full original model outputs.

If we are interested in those, we can reconstruct these selectively
by adding the aggregate outcomes as `extra_functionality` constraints (e.g. onwind p_nom_opt = XXX GW) and run the full optimisation.

Caveat: due to inaccuracies of the surrogate model, one may need load shedding or other forms of slack.

# Software

Chaospy

- https://github.com/jonathf/chaospy
- https://chaospy.readthedocs.io/en/master/

SAlib

- https://github.com/SALib/SALib
- https://salib.readthedocs.io/en/latest/index.html

Uncertainpy

- https://uncertainpy.readthedocs.io/en/latest/
- https://github.com/simetenn/uncertainpy

Polychaos.jl

- https://github.com/timueh/PolyChaos.jl
- https://timueh.github.io/PolyChaos.jl/dev/

SMT

- https://github.com/SMTorg/smt
- https://smt.readthedocs.io/en/latest/
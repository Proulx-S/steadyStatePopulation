# steadyStatePopulation

A minimal age-structured (Leslie-matrix / cohort-component) population model in MATLAB.
Age-dependent birth (fertility) and death (mortality) rates drive a population of yearly age
cohorts forward in time; fertility is auto-calibrated to a target net reproduction rate, which is
what makes the population settle into a genuine **steady state** rather than growing or shrinking
without bound — see [`populationModel.md`](populationModel.md) for why.

![age-structured population model](figures/populationModel.png)

- **top-left** — the age-dependent rates driving the model: mortality (red, rising with age, with
  an optional sharp senescence cliff late in life) and fertility (green, a triangular pulse over
  the fertile age window).
- **top-right** — total population over time, settling from its (deliberately non-stable) initial
  condition toward a steady long-run value.
- **bottom-left** — the initial (uniform) vs. final (stable) age distribution, directly comparing
  the starting and converged shapes.
- **bottom-right** — the full age-stratified population, normalized within each year, showing the
  age distribution's shape converge over time.

## Running it

Everything lives in a single script, [`doIt.m`](doIt.m). Open it in MATLAB and run it — no data,
no dependencies. All model parameters (age range, mortality/fertility shape, senescence cliff,
simulation horizon, target net reproduction rate) sit together in one `CONTROL PANEL` block near
the top; edit and re-run. The script ends by exporting the figure above to
`figures/populationModel.png`.

## Net reproduction rate (R0)

R0 is the expected number of offspring a person has over their entire lifetime, given the model's
mortality schedule — the standard demographic summary of whether a population is (more than)
replacing itself. It's computed from the *unscaled* fertility shape (the triangular pulse in the
top-left panel, before calibration) and the survivorship curve (the probability of surviving from
birth to each age, itself derived from the mortality curve): `R0 = sum(survivorship .*
fertilityShape)`. Fertility is then rescaled so this R0 hits `targetR0` (a `CONTROL PANEL`
parameter in `doIt.m`, default `1`): `fertility = fertilityShape * (targetR0 / R0)`. This
calibration is what the whole model hinges on — `targetR0 = 1` isn't just a plausible-looking
number, it's the exact condition for a long-run steady population (`>1` grows, `<1` declines); see
[`populationModel.md`](populationModel.md#8-why-r_01-gives-a-steady-state--the-eulerlotka-connection)
for the Euler–Lotka argument why.

## Formalism

[`populationModel.md`](populationModel.md) derives the model's equations in the order `doIt.m`
computes them — mortality, survival, fertility calibration, the cohort-update recursion — plus the
Euler–Lotka argument for why calibrating fertility to a net reproduction rate of 1 gives an exact
steady state, not an approximation.

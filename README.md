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
- **bottom-left** — the full age-stratified population, normalized within each year, showing the
  age distribution's shape converge over time.
- **bottom-right** — the initial (uniform) vs. final (stable) age distribution, directly comparing
  the starting and converged shapes.

## Running it

Everything lives in a single script, [`doIt.m`](doIt.m). Open it in MATLAB and run it — no data,
no dependencies. All model parameters (age range, mortality/fertility shape, senescence cliff,
simulation horizon, target net reproduction rate) sit together in one `CONTROL PANEL` block near
the top; edit and re-run. The script ends by exporting the figure above to
`figures/populationModel.png`.

## Formalism

[`populationModel.md`](populationModel.md) derives the model's equations in the order `doIt.m`
computes them — mortality, survival, fertility calibration, the cohort-update recursion — plus the
Euler–Lotka argument for why calibrating fertility to a net reproduction rate of 1 gives an exact
steady state, not an approximation.

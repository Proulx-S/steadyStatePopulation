# steadyStatePopulation

A minimal age-structured (Leslie-matrix / cohort-component) population model in MATLAB.
Age-dependent birth (fertility) and death (mortality) rates drive a population of yearly age
cohorts forward in time. By default (`rateSource='empirical'`) those rates are the world's actual
current age-specific fertility and mortality, so the model shows what they currently imply: the age
distribution's *shape* settles into the model's own stable form, but the *total* keeps growing,
since today's real rates aren't exactly at replacement. Switch to `rateSource='parameterized'` (see
[Rate source](#rate-source)) to instead calibrate fertility to a target net reproduction rate,
which is what makes the population settle into a genuine **steady state** — the property this
project is named for — rather than growing or shrinking without bound; see
[`populationModel.md`](populationModel.md) for why.

![age-structured population model](figures/populationModel.png)

- **top-left** — the age-dependent rates driving the model: mortality (red) and fertility (green)
  by age. With the default empirical rates these are the real, somewhat noisy age-band curves
  (a U-shaped mortality curve, low through adulthood then rising sharply past ~50; a fertility
  curve peaking in the late 20s); in parameterized mode these are smooth hand-tuned Gompertz and
  triangular shapes instead.
- **top-right** — total population over time. Its *shape* settles quickly; whether the *total*
  itself levels off depends on [rate source](#rate-source) — it doesn't with the default empirical
  rates (mild real-world growth, R0 ≈ 1.06), but does in parameterized mode with `targetR0=1`.
- **bottom-left** — the initial (current real-world) vs. final (model's own stable) age
  distribution, directly comparing the starting and converged shapes.
- **bottom-right** — the full age-stratified population, normalized within each year, showing the
  age distribution's shape converge over time.

## Running it

Everything lives in a single script, [`doIt.m`](doIt.m). Open it in MATLAB and run it — no data,
no dependencies. All model parameters (age range, mortality/fertility shape, senescence cliff,
simulation horizon, target net reproduction rate) sit together in one `CONTROL PANEL` block near
the top; edit and re-run. The script ends by exporting the figure above to
`figures/populationModel.png`.

## Rate source

`rateSource` (a `CONTROL PANEL` parameter, default `'empirical'`) picks between two
independent implementations of the mortality/fertility curves, selected by a `switch` in
`doIt.m` — everything downstream (survival, R0, the cohort simulation) uses whichever
`mortality(age)`/`fertility(age)` that switch produces, unchanged:

- **`'empirical'`** (default) — the current real-world age-specific rates for the World: fertility from the
  2023 age-specific fertility rate by 5-year age band, and mortality from the 2021 life table's
  probability of dying within each age band (both from
  [Our World in Data](https://ourworldindata.org/), sourced from UN World Population Prospects).
  Converted to this model's per-capita annual-rate convention (life-table band probability →
  annual hazard; per-1000-women fertility → per-capita, assuming a 50/50 sex split). Mortality
  beyond the data's oldest covered age (84) is extrapolated with a Gompertz (log-linear) fit to
  the last 5 empirical bands, rather than left undefined out to `ageMax`. Unlike the parameterized
  mode, fertility here is used *as-is* — not rescaled to `targetR0` — so whatever net reproduction
  rate the real data implies is what the model runs with (currently R0 ≈ 1.06, i.e. mild long-run
  growth, consistent with the world's fertility currently sitting slightly above replacement).
- **`'parameterized'`** — the hand-tuned Gompertz-plus-senescence-cliff mortality and triangular
  fertility described in [`populationModel.md`](populationModel.md), with fertility rescaled to
  hit `targetR0` exactly.

## Net reproduction rate (R0)

R0 is the expected number of offspring a person has over their entire lifetime, given the model's
mortality schedule — the standard demographic summary of whether a population is (more than)
replacing itself: `R0 = sum(survivorship .* fertilityShape)`, where `survivorship` is the
probability of surviving from birth to each age (itself derived from the mortality curve) and
`fertilityShape` is the *unscaled* fertility-by-age curve, before any calibration — NOT what's
plotted in the top-left panel, which shows the already-final `fertility`.

What happens next depends on [rate source](#rate-source). In `'parameterized'` mode, fertility is
rescaled so this R0 hits `targetR0` (a `CONTROL PANEL` parameter, default `1`):
`fertility = fertilityShape * (targetR0 / R0)`. This calibration is what that mode hinges on —
`targetR0 = 1` isn't just a plausible-looking number, it's the exact condition for a long-run
steady population (`>1` grows, `<1` declines); see
[`populationModel.md`](populationModel.md#8-why-r_01-gives-a-steady-state--the-eulerlotka-connection)
for the Euler–Lotka argument why. In `'empirical'` mode, `fertility = fertilityShape` is used
as-is — R0 is simply computed and reported (in the figure title), not targeted.

## Formalism

[`populationModel.md`](populationModel.md) derives the model's equations in the order `doIt.m`
computes them — mortality, survival, fertility calibration, the cohort-update recursion — plus the
Euler–Lotka argument for why calibrating fertility to a net reproduction rate of 1 gives an exact
steady state, not an approximation.

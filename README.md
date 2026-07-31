# steadyStatePopulation

A minimal age-structured (Leslie-matrix / cohort-component) population model in MATLAB.
Age-dependent birth (fertility) and death (mortality) rates drive a population of yearly age
cohorts forward in time. By default (`rateSource='parameterized'`) those rates are smooth Gompertz
(mortality) and triangular (fertility) curves — but not arbitrary ones: their shapes are fit by
least squares to the world's actual current age-specific rates (see [Rate
source](#rate-source)/[populationModel.md](populationModel.md) §10), so they closely track reality
while staying analytically simple. By default (`targetR0=[]`) they're also left unscaled, so the
model reports whatever net reproduction rate that fit implies (currently R0 ≈ 1.01) rather than
forcing one — set `targetR0=1` to instead calibrate fertility so the population settles into a
genuine **steady state**, the property this project is named for; see
[`populationModel.md`](populationModel.md) for why that particular calibration works.

![age-structured population model](figures/populationModel.png)

- **top-left** — the age-dependent rates driving the model: mortality (red) and fertility (green)
  by age, on a log axis (mortality alone spans ~4 orders of magnitude from age 0 to `ageMax`). Both
  rate sources show the same real U-shaped mortality pattern (high at birth, a childhood minimum,
  then rising) — the default parameterized curves are a smooth fit to it (populationModel.md §10);
  switch to `rateSource='empirical'` to see the real, unsmoothed age-band data instead.
- **top-right** — total population over time. Its *shape* settles quickly regardless; whether the
  *total* itself levels off depends on [`targetR0`](#net-reproduction-rate-r0) — it doesn't by
  default (`targetR0=[]`, mild growth at whatever R0 the rates actually imply), but does with
  `targetR0=1` set explicitly.
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

`rateSource` (a `CONTROL PANEL` parameter, default `'parameterized'`) picks between two
independent implementations of the mortality/fertility curves, selected by a `switch` in
`doIt.m` — everything downstream (survival, R0, the cohort simulation) uses whichever
`mortality(age)`/`fertility(age)` that switch produces, unchanged. Whether fertility gets rescaled
to `targetR0` is a *separate*, independent choice — see [Net reproduction rate
(R0)](#net-reproduction-rate-r0) — orthogonal to which rate source is picked here:

- **`'parameterized'`** (default) — smooth Gompertz-plus-senescence-cliff mortality (Eq. 1) and
  triangular fertility (Eq. 3) from [`populationModel.md`](populationModel.md). Its defaults aren't
  arbitrary: they're least-squares fit to the `'empirical'` curves below (§10), so even fully
  unscaled (`targetR0=[]`) they land close to the real net reproduction rate (R0 ≈ 1.01 vs. the
  real ≈ 1.06).
- **`'empirical'`** — the current real-world age-specific rates for the World: fertility from the
  2023 age-specific fertility rate by 5-year age band, and mortality from the 2021 life table's
  probability of dying within each age band (both from
  [Our World in Data](https://ourworldindata.org/), sourced from UN World Population Prospects).
  Converted to this model's per-capita annual-rate convention (life-table band probability →
  annual hazard; per-1000-women fertility → per-capita, assuming a 50/50 sex split). Mortality
  beyond the data's oldest covered age (84) is extrapolated with a Gompertz (log-linear) fit to
  the last 5 empirical bands, rather than left undefined out to `ageMax`.

## Net reproduction rate (R0)

R0 is the expected number of offspring a person has over their entire lifetime, given the model's
mortality schedule — the standard demographic summary of whether a population is (more than)
replacing itself: `R0 = sum(survivorship .* fertilityShape)`, where `survivorship` is the
probability of surviving from birth to each age (itself derived from the mortality curve) and
`fertilityShape` is the *unscaled* fertility-by-age curve, before any rescaling.

What happens next depends on `targetR0` (a `CONTROL PANEL` parameter, default empty `[]`) — this is
independent of [rate source](#rate-source), and works the same way for either:

- **`targetR0` empty** (default) — `fertility = fertilityShape` is used as-is; R0 is simply computed
  and reported (in the figure title), not targeted. Meaningful for both rate sources, since
  `fertilityShape` is already in real per-capita-rate units either way (parameterized mode's own
  peak-rate parameter, `fertilityPeakRate`, makes this so — see
  [`populationModel.md`](populationModel.md) Eq. 3).
- **`targetR0` a number** — fertility is rescaled so R0 hits it exactly:
  `fertility = fertilityShape * (targetR0 / R0)`. `targetR0 = 1` isn't just a plausible-looking
  number, it's the exact condition for a long-run steady population (`>1` grows, `<1` declines); see
  [`populationModel.md`](populationModel.md#8-why-r_01-gives-a-steady-state--the-eulerlotka-connection)
  for the Euler–Lotka argument why.

## Formalism

[`populationModel.md`](populationModel.md) derives the model's equations in the order `doIt.m`
computes them — mortality, survival, fertility calibration, the cohort-update recursion — plus the
Euler–Lotka argument for why calibrating fertility to a net reproduction rate of 1 gives an exact
steady state, not an approximation.

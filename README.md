# steadyStatePopulation

A minimal age-structured (Leslie-matrix / cohort-component; Leslie, 1945) population model in
MATLAB. Age-dependent birth (fertility) and death (mortality) rates drive a population of yearly
age cohorts forward in time. By default (`rateSource='conventional'`) those rates are literal
instances of two standard demographic laws — Siler's (1979) mortality model and Hadwiger's (1940)
fertility function — but not with arbitrary parameters: they're fit to the world's actual current
age-specific rates, shape AND net reproduction rate both (see [Rate
source](#rate-source)/[populationModel.md](populationModel.md) §10), so they closely track reality
while staying analytically simple. By default (`targetR0=[]`) they're also left unscaled, so the
model reports whatever net reproduction rate that fit implies (currently R0 ≈ 1.06, matching the
real data almost exactly) rather than forcing one — set `targetR0=1` to instead calibrate fertility
so the population settles into a genuine **steady state**, the property this project is named for;
see [`populationModel.md`](populationModel.md) for why that particular calibration works.

![age-structured population model](figures/populationModel.png)

- **top-left** — the age-dependent rates driving the model: mortality (red) and fertility (green)
  by age, on a log axis (mortality alone spans several orders of magnitude even just up to
  `ageMaxDisplay`). All three rate sources show the same real U-shaped mortality pattern (high at
  birth, a childhood minimum, then rising) — the default `'conventional'` curves are a smooth fit to
  it (populationModel.md §10); switch to `rateSource='empirical'` to see the real, unsmoothed
  age-band data instead, or `rateSource='rough'` for this project's own hand-built alternative fit.
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
no dependencies. All model parameters sit together in one `CONTROL PANEL` block near the top; edit
and re-run. The script ends by exporting the figure above to `figures/populationModel.png`.

Two things in the `CONTROL PANEL` worth knowing:

- **`ageMax` vs. `ageMaxDisplay`.** `ageMax` (large, 1000 by default) is how far the *simulation*
  tracks age classes — kept large so the plus-group (§6 of `populationModel.md`) holds negligible
  population, not to be looked at directly. `ageMaxDisplay` (100 by default) is a separate, purely
  cosmetic axis limit for the plots. They used to be conflated (`ageMax` alone controlled both), and
  the rough/conventional rate curves used to silently change shape with `ageMax` too (fixed — see
  `populationModel.md` §1a's note at the end).
- **Fitting happens live, every run.** Right after `ageVec` is built, `doIt.m` calls
  `empiricalRates` to get the current real-world mortality/fertility curves, then
  `fitMortalityRough`, `fitMortalityConventional`, `fitFertilityRough`, and
  `fitFertilityConventional` to fit all four forward models against them — regardless of which
  `rateSource` you've picked, so switching sources is instant and every fitted parameter struct
  (`roughMortParams`, `convMortParams`, `roughFertParams`, `convFertParams`) is always available in
  the workspace afterward to inspect or manually tweak, e.g. `roughMortParams.steepAge = 20;` (then
  re-run, or call `mortalityRough(ageVec, roughMortParams)` directly). If the empirical data in
  `empiricalRates` ever changes, every fitted default refits automatically — no separate manual
  re-fitting step (populationModel.md §10). Each of the four forward-model functions
  (`mortalityRough`, `mortalityConventional`, `fertilityRough`, `fertilityConventional`) also works
  standalone with just an age vector, falling back to hardcoded defaults from one such fit.

## Rate source

`rateSource` (a `CONTROL PANEL` parameter, default `'conventional'`) picks between three
independent implementations of the mortality/fertility curves, selected by a `switch` in
`doIt.m` — everything downstream (survival, R0, the cohort simulation) uses whichever
`mortality(age)`/`fertility(age)` that switch produces, unchanged. Whether fertility gets rescaled
to `targetR0` is a *separate*, independent choice — see [Net reproduction rate
(R0)](#net-reproduction-rate-r0) — orthogonal to which rate source is picked here:

- **`'conventional'`** (default) — literal instances of two standard demographic laws: Siler's
  (1979) competing-hazards mortality model (Eq. 1′) and Hadwiger's (1940) fertility function
  (Eq. 3′) from [`populationModel.md`](populationModel.md), implemented as `mortalityConventional`
  and `fertilityConventional`. Their defaults aren't arbitrary: they're a fit to the `'empirical'`
  curves below (§10) -- for fertility, shape AND net reproduction rate both, not shape alone -- so
  even fully unscaled (`targetR0=[]`) they land at essentially the real R0 (≈1.06).
- **`'rough'`** — this project's own hand-built alternative: a Gompertz-plus-senescence-cliff
  mortality curve (Eq. 1) and a triangular fertility pulse (Eq. 3), implemented as `mortalityRough`
  and `fertilityRough`. Fit to the same `'empirical'` curves the same way as `'conventional'`
  (§10); tracks the data about as well for mortality, slightly less well for fertility (§10) — kept
  as a simpler, closed-form comparison point.
- **`'empirical'`** — the current real-world age-specific rates for the World: fertility from the
  2023 age-specific fertility rate by 5-year age band (Our World in Data, 2026a, sourced from the
  UN Population Division, 2024a), and mortality from the 2021 life table's probability of dying
  within each age band (Our World in Data, 2026b, sourced from the **WHO Global Health Observatory**
  — not UN WPP; these two datasets, despite both coming through Our World in Data, have different
  original sources). Converted to this model's per-capita annual-rate convention (life-table band
  probability → annual hazard; per-1000-women fertility → per-capita, assuming a 50/50 sex split,
  itself a simplification — the real figure is closer to 48.8% female at birth; UN Population
  Division, 2024b). Mortality beyond the data's oldest covered age (84) is extrapolated with a
  Gompertz (log-linear) fit to the last 5 empirical bands, rather than left undefined out to
  `ageMax`. Full citations in [References](#references). This is the fitting *target* for both
  `'conventional'` and `'rough'` above (§10), via the `empiricalRates` function.

## Net reproduction rate (R0)

R0 is the expected number of offspring a person has over their entire lifetime, given the model's
mortality schedule — the standard demographic summary of whether a population is (more than)
replacing itself: `R0 = sum(survivorship .* fertilityShape)`, where `survivorship` is the
probability of surviving from birth to each age (itself derived from the mortality curve) and
`fertilityShape` is the *unscaled* fertility-by-age curve, before any rescaling.

What happens next depends on `targetR0` (a `CONTROL PANEL` parameter, default empty `[]`) — this is
independent of [rate source](#rate-source), and works the same way for either:

- **`targetR0` empty** (default) — `fertility = fertilityShape` is used as-is; R0 is simply computed
  and reported (in the figure title), not targeted. Meaningful for all three rate sources, since
  `fertilityShape` is already in real per-capita-rate units either way (`'rough'`'s own peak-rate
  parameter, `fertilityPeakRate`, and `'conventional'`'s Hadwiger scale parameter `a` both make this
  so — see [`populationModel.md`](populationModel.md) Eqs. 3, 3′).
- **`targetR0` a number** — fertility is rescaled so R0 hits it exactly:
  `fertility = fertilityShape * (targetR0 / R0)`. `targetR0 = 1` isn't just a plausible-looking
  number, it's the exact condition for a long-run steady population (`>1` grows, `<1` declines); see
  [`populationModel.md`](populationModel.md#8-why-r_01-gives-a-steady-state--the-eulerlotka-connection)
  for the Euler–Lotka argument why.

## Relation to the literature

The core method — Leslie-matrix cohort projection (Leslie, 1945) and the Euler–Lotka steady-state
condition (Sharpe & Lotka, 1911) — is standard, unmodified demographic theory. The curve *shapes*
feeding into it come in two flavors: `rateSource='rough'`'s mortality formula is
Siler (1979)/Heligman & Pollard (1980)-*inspired* (their classic infant + background + senescent
competing-hazards template) but not a literal instance of either, and its triangular fertility
shape is a simplification of the smooth, right-skewed curves standard in the field (Hadwiger, 1940;
Coale & Trussell, 1974; Schmertmann, 2003); `rateSource='conventional'` (the default) instead uses a
*literal* Siler mortality law and Hadwiger fertility function, matching the data about as well or
better (populationModel.md §10) while being direct instances of named, standard demographic laws
rather than simplified variants of them. [`populationModel.md`](populationModel.md) §11 works
through each comparison in detail — including sex structure, still unimplemented — and what's left
to bring the model even closer to current demographic practice.

## Formalism

[`populationModel.md`](populationModel.md) derives the model's equations in the order `doIt.m`
computes them — mortality, survival, fertility calibration, the cohort-update recursion — plus the
Euler–Lotka argument for why calibrating fertility to a net reproduction rate of 1 gives an exact
steady state, not an approximation.

## References

Coale, A. J., & Trussell, T. J. (1974). Model fertility schedules: Variations in the age structure
of childbearing in human populations. *Population Index*, 40(2), 185–258.

Hadwiger, H. (1940). Eine analytische Reproduktionsfunktion für biologische Gesamtheiten.
*Scandinavian Actuarial Journal*, 1940(3–4), 101–113.

Heligman, L., & Pollard, J. H. (1980). The age pattern of mortality. *Journal of the Institute of
Actuaries*, 107(1), 49–80.

Leslie, P. H. (1945). On the use of matrices in certain population mathematics. *Biometrika*,
33(3), 183–212.

Schmertmann, C. P. (2003). A system of model fertility schedules with graphically intuitive
parameters. *Demographic Research*, 9, 81–110.

Sharpe, F. R., & Lotka, A. J. (1911). A problem in age-distribution. *Philosophical Magazine*,
21(124), 435–438.

Siler, W. (1979). A competing-risk model for animal mortality. *Ecology*, 60(4), 750–757.

UN Population Division. (2024a). *World Population Prospects 2024, Online Edition*.
https://population.un.org/wpp/

UN Population Division. (2024b). Sex ratio at birth [Data set]. *World Population Prospects 2024*
– processed by Our World in Data. https://ourworldindata.org/grapher/sex-ratio-at-birth

Our World in Data. (2026a). Fertility rate by age group [Data set]. Sourced from UN Population
Division, World Population Prospects (2024a).
https://ourworldindata.org/grapher/fertility-rate-by-age-group

Our World in Data. (2026b). Probability of dying by age [Data set]. Sourced from World Health
Organization, Global Health Observatory.
https://ourworldindata.org/grapher/probability-of-dying-by-age

Central Intelligence Agency. (2021). *The World Factbook: Age structure (World)*, via IndexMundi.
https://www.indexmundi.com/world/age_structure.html

Full derivations and per-claim citations: [`populationModel.md`](populationModel.md).

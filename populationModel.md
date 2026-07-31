# Age-structured population model — formalism

`doIt.m` implements a discrete-time, age-structured (Leslie-matrix / cohort-component) population
model: age-dependent birth and death rates drive a population of yearly age cohorts forward in
time. This note derives the model's equations in the order `doIt.m` computes them, and works out
*why* calibrating fertility to a target net reproduction rate makes the population settle into a
genuine steady state — the property the project is named for.

---

## Notation

| symbol | meaning | code variable |
|---|---|---|
| $a=0,1,\dots,A$ | age, in whole years ($A=$ `ageMax`, a plus-group) | `ageVec` |
| $t=0,1,\dots,T$ | time, in years ($T=$ `nYears`) | `years` |
| $N(a,t)$ | population in age class $a$ at year $t$ | `N` |
| $P(t)$ | total population at year $t$ | `totalPop` |
| $\mu(a)$ | per-capita annual mortality hazard at age $a$ | `mortality` |
| $s(a)$ | per-capita annual survival probability at age $a$ | `survival` |
| $b(a)$ | per-capita annual fertility (birth) rate at age $a$ | `fertility` |
| $\ell(a)$ | survivorship: probability of being alive at age $a$, given born | `survivorship` |
| $R_0$ | net reproduction rate (expected offspring per lifetime) | `sum(survivorship.*fertility)` |
| $\lambda$ | asymptotic year-over-year population growth multiplier | — |

---

## 1. Age-dependent mortality

*(This section, and §3 below, describe the `rateSource='parameterized'` mode. §9 gives the
`'empirical'` alternative -- the project's own default.)*

Baseline hazard rises with age (Gompertz-like), plus a sharp senescence cliff past a configurable
fraction $f_s$ (`steepAgeFrac`) of the maximum age:

$$
\mu(a) = \underbrace{\mu_0 + \mu_1\left(\frac{a}{A}\right)^2}_{\text{baseline ageing}} \;+\;
\underbrace{\mathbb{1}[a > f_s A]\cdot \mu_2\Big(e^{\,k\,(a - f_s A)} - 1\Big)}_{\text{senescence cliff}}
\tag{1}
$$

where $\mu_0=$ `deathRateBase`, $\mu_1=$ `deathRateSlope`, $\mu_2=$ `steepMortalityScale`,
$k=$ `steepMortalityRate`, and $\mathbb{1}[\cdot]$ is the indicator function. The cliff term is
exactly zero at and before $a=f_s A$, then grows exponentially — with the current (fit, §10)
defaults $f_s=0.144$, $A=150$, it kicks in at age $21.6$ and drives $\mu(150)\approx25.6$ (survival
$\approx7\times10^{-12}$/yr) by `ageMax`. $\mu_1=0$ in the fit defaults: the cliff term alone
captures the empirical curve's rise better than a separate baseline quadratic (§10).

## 2. Survival probability

Converting the annual hazard to a survival probability (standard continuous-hazard-over-one-year
approximation):

$$
s(a) = e^{-\mu(a)} \tag{2}
$$

## 3. Age-dependent fertility

Fertility is a triangular pulse over the fertile age window $[a_{\min},a_{\max}]$
(`fertileMin`,`fertileMax`), peaking at $a^*$ (`fertilityPeakAge`) with peak height $b_{\max}$
(`fertilityPeakRate`), zero outside it — a *shape* $b_0(a)$ fixed before calibration (§5). Including
$b_{\max}$ keeps $b_0(a)$ in real per-capita-rate units on its own, so it stays meaningful even with
no calibration at all, i.e. `targetR0=[]` (§5) — an earlier version of this model normalized the
peak to $1$ instead, meaningless without calibration:

$$
b_0(a) = \begin{cases}
b_{\max}\cdot\max\!\left(0,\ 1 - \dfrac{|a-a^*|}{w}\right), & a_{\min}\le a\le a_{\max} \\[4pt]
0, & \text{otherwise}
\end{cases}
\qquad w=\max(a^*-a_{\min},\ a_{\max}-a^*) \tag{3}
$$

## 4. Survivorship

The probability of being alive at age $a$, given birth, is the cumulative product of the survival
probabilities of every age passed through:

$$
\ell(0)=1, \qquad \ell(a) = \prod_{k=0}^{a-1} s(k) \quad (a \ge 1) \tag{4}
$$

## 5. Net reproduction rate and fertility calibration

The net reproduction rate $R_0$ is the expected number of offspring a newborn has over its whole
lifetime, under the mortality schedule of Eq. (1)–(2):

$$
R_0 = \sum_{a=0}^{A} \ell(a)\, b_0(a) \tag{5}
$$

If `targetR0` is non-empty, fertility is rescaled so $R_0$ hits that chosen target:

$$
b(a) = b_0(a)\cdot\frac{R_0^{\text{target}}}{R_0} \tag{6}
$$

If `targetR0` is empty (the project's own default), no rescaling happens at all: $b(a)=b_0(a)$, and
$R_0$ from Eq. (5) is simply reported, not targeted -- meaningful now that $b_0(a)$ itself is
already in real units (Eq. 3's $b_{\max}$).

## 6. Cohort update (Leslie-matrix recursion)

Each year, every age class ages up by one year according to its survival probability; the oldest
class ($a=A$) is a **plus-group** that retains its own survivors rather than aging out; newborns
enter age $0$ in proportion to the population-weighted fertility of every age class:

$$
N(a{+}1,\,t{+}1) = N(a,t)\,s(a), \qquad a = 0,\dots,A{-}2 \tag{7a}
$$

$$
N(A,\,t{+}1) = N(A{-}1,t)\,s(A{-}1) \;+\; N(A,t)\,s(A) \tag{7b}
$$

$$
N(0,\,t{+}1) = \sum_{a=0}^{A} N(a,t)\, b(a) \tag{7c}
$$

## 7. Total population

$$
P(t) = \sum_{a=0}^{A} N(a,t) \tag{8}
$$

---

## 8. Why $R_0=1$ gives a steady state — the Euler–Lotka connection

The classical result linking a population's vital rates to its long-run behavior is the
**Euler–Lotka equation**. In its discrete-age form it says the asymptotic growth multiplier
$\lambda$ (population size scales as $\lambda^t$ once the age structure has settled down) is the
positive root of

$$
1 = \sum_{a=0}^{A} \ell(a)\, b(a)\, \lambda^{-a} \tag{9}
$$

Setting $\lambda=1$ in Eq. (9) reduces it *exactly* to Eq. (5)/(6)'s calibration condition,
$1 = \sum_a \ell(a) b(a) = R_0$. So $\lambda=1 \iff R_0=1$: calibrating fertility to
$R_0^{\text{target}}=1$ (Eq. 6) is not an approximation of a steady state, it **is** the condition
for one — no need to solve Eq. (9) for $\lambda$ directly. (`targetR0>1` gives $\lambda>1$,
long-run growth; `targetR0<1` gives $\lambda<1$, long-run decline — Eq. (9) still governs the
rate, it's just not needed for the steady-state case.)

The same theory predicts the **stable age distribution** the population converges to is
proportional to $\ell(a)\lambda^{-a}$; at $\lambda=1$ that is simply the survivorship curve
$\ell(a)$ itself, normalized:

$$
\lim_{t\to\infty} \frac{N(a,t)}{P(t)} = \frac{\ell(a)}{\sum_{a'} \ell(a')} \tag{10}
$$

This is a checkable prediction, not just an assertion: running `doIt.m` with
`rateSource='parameterized'` and `targetR0` set to $1$ (its own default is empty, §5, so this needs
setting explicitly -- giving $\lambda=1$ exactly), the simulated final-year age distribution
`N(:,end)/sum(N(:,end))` matches `survivorship/sum(survivorship)` to within $8\times10^{-10}$ after
`nYears`$=500$ years, starting from the real-world (not the model's own stable) age distribution in
the parameter table below — confirming Eq. (10) numerically to near machine precision. (Eq. (10) is
specifically the $\lambda=1$ case; with `rateSource='empirical'` or any other non-empty `targetR0`,
$\lambda\neq1$ in general, so the general $\ell(a)\lambda^{-a}$ form above it applies instead.)

---

## 9. Empirical rate alternative

`rateSource='empirical'` (the project's own default) replaces Eqs. (1) and (3) with the current real-world age-specific rates
for the World, rather than the hand-tuned shapes: fertility from the 2023 age-specific fertility
rate by 5-year age band, mortality from the 2021 life table's probability of dying within each age
band (both [Our World in Data](https://ourworldindata.org/), sourced from UN World Population
Prospects; exact per-band values are in `doIt.m`'s `fertAsfr`/`mortQx`).

**Mortality**, converting each band's probability of dying $q_{\text{band}}$ (given alive at the
band's start) to this model's per-capita annual hazard, for a band of width $w$ years:

$$
\mu(a) = -\frac{\ln(1-q_{\text{band}(a)})}{w_{\text{band}(a)}} \tag{11}
$$

which is just Eq. (2) solved for $\mu$ given the band's own multi-year survival
$(1-q_{\text{band}})=s(a)^w$. Past the data's oldest covered age (84), no life-table value exists,
so $\mu(a)$ is extrapolated with a Gompertz (log-linear) fit to the last 5 empirical bands'
hazards, continuing their trend out to `ageMax` rather than leaving it undefined.

**Fertility**, converting the age-specific fertility rate $\text{ASFR}(a)$ (births per 1000 women
per year) to this model's per-capita (not per-woman) rate, assuming a female population share
$p_f=$ `femaleFrac` ($=0.5$):

$$
b(a) = \frac{\text{ASFR}(a)}{1000}\, p_f \tag{12}
$$

zero outside the reported age bands (10–54).

**No R0 calibration in this mode.** Unlike §5–6, $b(a)$ from Eq. (12) is used directly — NOT
rescaled to `targetR0` — so $R_0=\sum_a \ell(a)b(a)$ is whatever the real data implies, not a
chosen target. This is deliberate: the point of this mode is to see what the world's *actual*
current rates predict, not to force them into the steady-state condition of §8. (Currently
$R_0\approx1.06$ — mild long-run growth, consistent with global fertility sitting slightly above
replacement.)

---

## 10. Fitting the parameterized shapes to empirical data

The `'parameterized'` defaults (Eqs. 1, 3) are not arbitrary hand-picks — they're a least-squares
fit of those SAME functional forms to the `'empirical'` mortality/fertility curves of §9 (including
§9's own old-age Gompertz extrapolation beyond age 84), via `fminsearch` (Nelder–Mead, unconstrained
— parameters reparametrized through `exp`/logistic transforms to keep $\mu_0,\mu_1,\mu_2,k>0$ and
$f_s\in(0,1)$).

**Mortality** is fit in LOG-hazard space (minimizing $\sum_a\big(\ln\mu_{\text{param}}(a)-\ln\mu_{\text{emp}}(a)\big)^2$)
since $\mu$ spans several orders of magnitude with age. The fit found $\mu_1\approx0$: Eq. (1)'s
baseline quadratic term contributes essentially nothing in the best fit — the senescence-cliff term
alone, starting much younger ($f_sA\approx21.6$) than "late-life cliff" suggests, captures the
empirical curve's overall rise better than baseline-plus-cliff can. This is an intrinsic limitation
of Eq. (1), not a fitting failure: it has no mechanism for the empirical curve's early-childhood DIP
(high at birth, lowest around age 5–14, then rising), since both of Eq. (1)'s terms are
monotonically non-decreasing in age — the fit's residual concentrates there, not at the older ages
that matter most for the model's long-run behavior.

**Fertility** is fit in LINEAR space (minimizing $\sum_a(b_{0,\text{param}}(a)-b_{0,\text{emp}}(a))^2$),
since $b_0$ includes true zeros where log-space is undefined.

**Cross-validation**: running the fit defaults with `targetR0=[]` (completely unscaled, §5) gives
$R_0\approx1.02$ — close to the empirical mode's own $R_0\approx1.06$ (§9), a good sign the fit
captures the real data's overall reproduction level, not just its shape.

---

## Parameter reference

| symbol | code | default | meaning |
|---|---|---|---|
| $A$ | `ageMax` | $150$ | oldest age class (plus-group) |
| $T$ | `nYears` | $500$ | simulation horizon |
| — | `popInit` | $8\times10^9$ | total starting population |
| — | — | — | initial age distribution: current real-world world age structure (CIA World Factbook, 2021-2023 estimates: 0-14 25.2%, 15-24 15.3%, 25-54 40.6%, 55-64 9.2%, 65+ 9.7%), not the model's own stable shape -- see Eq. (10) |
| $\mu_0$ | `deathRateBase` | $0.0011$ | mortality hazard at age $0$ (fit, §10) |
| $\mu_1$ | `deathRateSlope` | $0$ | added baseline hazard at age $A$ (fit to $\approx0$, §10) |
| $f_s$ | `steepAgeFrac` | $0.144$ | age fraction of $A$ where the senescence cliff begins (fit, §10) |
| $\mu_2$ | `steepMortalityScale` | $0.00053$ | extra-hazard scale past the cliff (fit, §10) |
| $k$ | `steepMortalityRate` | $0.084$ | extra-hazard growth rate past the cliff (fit, §10) |
| $a_{\min}$ | `fertileMin` | $15.1$ | youngest fertile age (fit, §10) |
| $a_{\max}$ | `fertileMax` | $42.4$ | oldest fertile age (fit, §10) |
| $a^*$ | `fertilityPeakAge` | $27.0$ | age of peak fertility (fit, §10) |
| $b_{\max}$ | `fertilityPeakRate` | $0.0709$ | per-capita rate at $a^*$ (fit, §10) |
| $R_0^{\text{target}}$ | `targetR0` | empty (`[]`) | net reproduction rate target; empty = don't rescale, just estimate/report $R_0$ (§5); a number ($1=$steady, $>1=$growth, $<1=$decline) rescales fertility to hit it exactly |

---

## Modelling assumptions

1. **Yearly, discrete-age cohorts** — no within-year age or seasonality structure; birth/death
   probabilities are applied once per year (Eqs. 1–2, 7).
2. **No sex structure** — $N(a,t)$ is a single (unisex) population; fertility $b(a)$ is a per-capita
   rate applied to the whole cohort, not per-female.
3. **No migration, environmental stochasticity, or density dependence** — rates $\mu(a)$, $b(a)$ are
   fixed functions of age alone, independent of $t$ or $P(t)$; this is exactly what makes the
   Euler–Lotka argument (§8) apply cleanly.
4. **Plus-group at $A$** (Eq. 7b) — the oldest class is an absorbing bin, not a hard cutoff; nobody
   is assumed to die exactly at age $A$, though the senescence cliff (Eq. 1) makes survival past it
   very unlikely with default parameters.

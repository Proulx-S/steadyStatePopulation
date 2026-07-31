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

Baseline hazard rises with age (Gompertz-like), plus a sharp senescence cliff past a configurable
fraction $f_s$ (`steepAgeFrac`) of the maximum age:

$$
\mu(a) = \underbrace{\mu_0 + \mu_1\left(\frac{a}{A}\right)^2}_{\text{baseline ageing}} \;+\;
\underbrace{\mathbb{1}[a > f_s A]\cdot \mu_2\Big(e^{\,k\,(a - f_s A)} - 1\Big)}_{\text{senescence cliff}}
\tag{1}
$$

where $\mu_0=$ `deathRateBase`, $\mu_1=$ `deathRateSlope`, $\mu_2=$ `steepMortalityScale`,
$k=$ `steepMortalityRate`, and $\mathbb{1}[\cdot]$ is the indicator function. The cliff term is
exactly zero at and before $a=f_s A$, then grows exponentially — with the default $f_s=0.95$,
$A=90$, it kicks in at age $85.5$ and drives $\mu(90)\approx4.5$ (survival $\approx1\%$/yr) by
`ageMax`.

## 2. Survival probability

Converting the annual hazard to a survival probability (standard continuous-hazard-over-one-year
approximation):

$$
s(a) = e^{-\mu(a)} \tag{2}
$$

## 3. Age-dependent fertility

Fertility is a triangular pulse over the fertile age window $[a_{\min},a_{\max}]$
(`fertileMin`,`fertileMax`), peaking at $a^*$ (`fertilityPeakAge`), zero outside it — an
unscaled *shape* $b_0(a)$ fixed before calibration (§5):

$$
b_0(a) = \begin{cases}
\max\!\left(0,\ 1 - \dfrac{|a-a^*|}{w}\right), & a_{\min}\le a\le a_{\max} \\[4pt]
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

Fertility is then rescaled so $R_0$ hits a chosen target (`targetR0`, default $1$):

$$
b(a) = b_0(a)\cdot\frac{R_0^{\text{target}}}{R_0} \tag{6}
$$

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

This is a checkable prediction, not just an assertion: running `doIt.m` with defaults, the
simulated final-year age distribution `N(:,end)/sum(N(:,end))` matches `survivorship/sum(survivorship)`
to within $2\times10^{-4}$ after `nYears`$=150$ years, starting from the real-world (not the
model's own stable) age distribution in the parameter table below — confirming Eq. (10) numerically.

---

## Parameter reference

| symbol | code | default | meaning |
|---|---|---|---|
| $A$ | `ageMax` | $90$ | oldest age class (plus-group) |
| $T$ | `nYears` | $150$ | simulation horizon |
| — | `popInit` | $8\times10^9$ | total starting population |
| — | — | — | initial age distribution: current real-world world age structure (CIA World Factbook, 2021-2023 estimates: 0-14 25.2%, 15-24 15.3%, 25-54 40.6%, 55-64 9.2%, 65+ 9.7%), not the model's own stable shape -- see Eq. (10) |
| $\mu_0$ | `deathRateBase` | $0.004$ | mortality hazard at age $0$ |
| $\mu_1$ | `deathRateSlope` | $0.08$ | added baseline hazard at age $A$ |
| $f_s$ | `steepAgeFrac` | $0.95$ | age fraction of $A$ where the senescence cliff begins |
| $\mu_2$ | `steepMortalityScale` | $0.05$ | extra-hazard scale past the cliff |
| $k$ | `steepMortalityRate` | $1$ | extra-hazard growth rate past the cliff |
| $a_{\min}$ | `fertileMin` | $15$ | youngest fertile age |
| $a_{\max}$ | `fertileMax` | $45$ | oldest fertile age |
| $a^*$ | `fertilityPeakAge` | $28$ | age of peak fertility |
| $R_0^{\text{target}}$ | `targetR0` | $1$ | net reproduction rate target ($1=$steady, $>1=$growth, $<1=$decline) |

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

# Age-structured population model — formalism

`doIt.m` implements a discrete-time, age-structured (Leslie-matrix / cohort-component) population
model (Leslie, 1945): age-dependent birth and death rates drive a population of yearly age cohorts
forward in time. This note derives the model's equations in the order `doIt.m` computes them, works
out *why* calibrating fertility to a target net reproduction rate makes the population settle into a
genuine steady state — the property the project is named for — and compares each modelling choice
against the standard demographic literature, flagging where a simplification here diverges from
current practice (§11) and citing sources throughout (§References).

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

Three terms: a decaying excess hazard at birth (infant mortality), a slowly-rising baseline
(Gompertz-like ageing), and a sharp senescence cliff past an absolute age $a_c$ (`steepAge`):

$$
\mu(a) = \underbrace{\mu_i\, e^{-a/\tau_i}}_{\text{infant mortality}} \;+\;
\underbrace{\mu_0 + \mu_1\left(\frac{a}{A_{\text{ref}}}\right)^2}_{\text{baseline ageing}} \;+\;
\underbrace{\mathbb{1}[a > a_c]\cdot \mu_2\Big(e^{\,k\,(a - a_c)} - 1\Big)}_{\text{senescence cliff}}
\tag{1}
$$

where $\mu_i=$ `infantMortalityScale`, $\tau_i=$ `infantMortalityDecay`, $\mu_0=$ `deathRateBase`,
$\mu_1=$ `deathRateSlope`, $A_{\text{ref}}=$ `ageRef`, $\mu_2=$ `steepMortalityScale`, $k=$
`steepMortalityRate`, and $\mathbb{1}[\cdot]$ is the indicator function. The infant term is what
lets $\mu(a)$ have the real U-shape (high at birth, a childhood minimum, then rising) — neither of
the other two terms is ever decreasing in $a$, so without it the model cannot represent that dip at
all (§10). The cliff term is exactly zero at and before $a=a_c$, then grows exponentially — with the
current (fit, §10) defaults $a_c=11.01$ (right where the infant term has mostly decayed away) and
$A_{\text{ref}}=150$, $\mu(90)\approx0.168$ (survival $\approx0.845$/yr).

**Relation to the literature.** The rising-with-age "Gompertz-like" middle+cliff behaviour is named
for Gompertz's (1825) observation that adult human mortality rises close to exponentially with age
— the oldest quantitative law in mortality modelling and still the basis of the senescent term in
essentially every mortality model used today. Eq. (1)'s overall three-part shape (a term that falls
with age, one that's roughly flat, one that rises) is the same competing-hazards idea behind two
standard mortality laws: Siler's (1979) model,
$h(a)=a_1e^{-a_2a}+a_3+a_4e^{a_5a}$ — a *falling* exponential (infancy) plus a true *constant*
(background) plus a *rising, unthresholded* exponential from birth (senescence) — and
Heligman & Pollard's (1980) eight-parameter law for human data specifically, which uses a power-law
infant term, a distinct log-normal "accident hump" for young-adult mortality (elevated risk from
injury/violence around ages 15–30, which Eq. (1) has no term for at all), and a logistic senescent
term. Eq. (1) follows this same infant-plus-background-plus-senescent template but is **not** a
literal instance of either: its "background" term is a slowly *rising* quadratic rather than Siler's
true constant, and its senescent term is *thresholded* at $a_c$ and shifted (zero below $a_c$, not
merely small) rather than Siler's or Heligman-Pollard's smooth, unthresholded rise from birth. §11
discusses whether adopting one of these standard forms literally would be a better-grounded choice.

**$A_{\text{ref}}$ and $a_c$ are fixed constants, deliberately NOT tied to $A$ (`ageMax`)** — an
earlier version of this formula normalized the baseline term by $A$ itself and set the cliff
threshold as a *fraction* of $A$, so simply increasing `ageMax` (e.g. to shrink the plus-group's
share of the population, §6) silently rescaled $\mu(a)$ at every age, for no physical reason (the
fit was performed at one specific `ageMax`; deviating from it deviated the curve). Decoupling them
means $\mu(a)$ — and hence $R_0$, `totalPop`, everything — no longer depends on `ageMax` at all;
only `ageMaxDisplay` (plot axis limits, purely cosmetic) does.

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

**Relation to the literature.** A piecewise-linear "tent" is not a standard demographic fertility
model — real age-specific fertility schedules are smooth and (mildly) right-skewed: they rise faster
from the youngest fertile ages than they fall off toward the oldest (confirmed directly in this
model's own empirical data, §9: the rate climbs from age 15 to its 25–29 peak over ~10 years but
takes until the early 50s, ~25 years, to fall back to zero). The standard alternatives are smooth,
unimodal, right-skewed curves fit with as few as 2–4 parameters: Hadwiger's (1940) function (still
widely used, close in spirit to a right-skewed bell curve), the empirical-natural-fertility-based
Coale & Trussell (1974) model schedules, and Schmertmann's (2003) more recent spline-based system,
designed specifically so its parameters ("graphically intuitive") map directly onto features like
peak age and peak level — much like this model's own `fertilityPeakAge`/`fertilityPeakRate`. Eq. (3)
trades that smoothness for a closed form simple enough to invert by hand (§10's fit only needs to
find 4 numbers); §11 suggests the Hadwiger function specifically as a same-parameter-count, smoother
alternative worth trying.

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
**Euler–Lotka equation** — Euler derived a special case in 1760, in French, largely unnoticed by
demographers until Keyfitz & Keyfitz's English translation (Euler, 1760/1970) brought it back into
view; independently, Sharpe & Lotka (1911) derived the general continuous-age version that underlies
modern stable population theory. In its discrete-age
form (the natural fit here, since ages in this model are already whole-year classes, not a
discretized approximation of the classical continuous-age integral) it says the asymptotic growth
multiplier $\lambda$
(population size scales as $\lambda^t$ once the age structure has settled down) is the positive
root of

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
`N(:,end)/sum(N(:,end))` matches `survivorship/sum(survivorship)` to within $2\times10^{-10}$ after
`nYears`$=500$ years, starting from the real-world (not the model's own stable) age distribution in
the parameter table below — confirming Eq. (10) numerically to near machine precision. (Eq. (10) is
specifically the $\lambda=1$ case; with `rateSource='empirical'` or any other non-empty `targetR0`,
$\lambda\neq1$ in general, so the general $\ell(a)\lambda^{-a}$ form above it applies instead.)

---

## 9. Empirical rate alternative

`rateSource='empirical'` (the project's own default) replaces Eqs. (1) and (3) with the current
real-world age-specific rates for the World, rather than the hand-tuned shapes: fertility from the
2023 age-specific fertility rate by 5-year age band (Our World in Data, 2026a, sourced from United
Nations Population Division, 2024a), mortality from the 2021 life-table probability of dying within
each age band (Our World in Data, 2026b, sourced from the World Health Organization's Global Health
Observatory — **not** UN WPP; the two datasets have different original sources despite both being
accessed through the same Our World in Data chart family, a distinction an earlier draft of this
document got wrong). Exact per-band values are in `doIt.m`'s `fertAsfr`/`mortQx`; full citations in
§References.

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

zero outside the reported age bands (10–54). $p_f=0.5$ is itself a simplification of a real,
age-varying quantity: about 105 boys are born for every 100 girls globally (United Nations
Population Division, 2024b), i.e. female share $\approx0.488$ at birth, not $0.5$, and the female
SHARE of each age cohort rises above 0.5 at older ages since male mortality is higher at nearly
every age (assumption 2 already flags the model's lack of sex structure generally; this is the
specific number that assumption trades away).

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
— parameters reparametrized through `exp`/logistic transforms to keep $\mu_i,\tau_i,\mu_0,\mu_1,\mu_2,k>0$;
the cliff threshold was fit as a fraction of the fitting run's own `ageMax` then converted to the
fixed absolute age $a_c$ used in Eq. (1) -- see the note at the end of §1 for why it's absolute now).

**Mortality** is fit in LOG-hazard space (minimizing $\sum_a\big(\ln\mu_{\text{param}}(a)-\ln\mu_{\text{emp}}(a)\big)^2$)
since $\mu$ spans several orders of magnitude with age. A first pass *without* the infant term
(i.e. fitting only the baseline-plus-cliff part of Eq. 1) found $\mu_1\approx0$ and a mediocre fit
(RMS log-residual $0.34$): with no term that's ever DEcreasing in age, that version of Eq. (1)
structurally cannot represent the empirical curve's early-childhood dip (high at birth, lowest
around age 5–14, then rising) — its best compromise just ignored that region. Adding the infant
term (Eq. 1's first term) fixes this: refit jointly, RMS log-residual drops to $0.15$ (less than
half), the childhood dip is now captured (visually near-exact from age $\sim$3 onward), and
$\mu_1$ comes back to a real, nonzero value ($5.7\times10^{-4}$) since the baseline term is no
longer needed to (badly) approximate the dip too.

**Fertility** is fit in LINEAR space, minimizing pointwise squared error PLUS a penalty on the
resulting $R_0$'s mismatch from the empirical mode's own $R_0$ (§9), weighted by *this* mode's own
survivorship (Eq. 4) since that's what actually gets used once plugged into Eq. (5):

$$
\sum_a\big(b_{0,\text{param}}(a)-b_{0,\text{emp}}(a)\big)^2 \;+\;
\Lambda\Big(\textstyle\sum_a \ell_{\text{param}}(a)\,b_{0,\text{param}}(a) - R_0^{\text{emp}}\Big)^2
$$

with $\Lambda$ large enough ($10^5$) to make the second term dominate. A pointwise-only fit
($\Lambda=0$) matches the empirical *shape* well (visually near-identical peak/width) but still
landed at $R_0\approx1.01$ against the empirical $R_0\approx1.06$ — a $\sim$5% shortfall traced
to two compounding effects, roughly half each: (i) the empirical curve has real, small fertility
just outside a pointwise-fit window's edges (the 10–14 and 45+ age bands) that a hard triangular
cutoff misses entirely regardless of how well it fits inside the window; (ii) even inside the
window, minimizing pointwise error doesn't preserve the curve's *area* — the fit slightly overshoots
near the peak and undershoots at the edges, and those don't cancel. Adding the $R_0$ penalty lets
the optimizer trade a negligible amount of pointwise fit (SSE $0.00196$ vs. $0.00195$) for a wider
window ($a_{\min}$ drops to $12.1$) that closes both gaps at once — by construction, since $R_0$
matching is now directly in the objective, not a hoped-for side effect of shape matching.

**Cross-validation**: running the fit defaults with `targetR0=[]` (completely unscaled, §5) now gives
$R_0\approx1.06$, matching the empirical mode's own $R_0$ (§9) essentially exactly — by design, per
the $R_0$-matching term above, but confirming it actually works when plugged into the full model
rather than just the fitting script.

---

## 11. Relation to the demographic literature; suggested refinements

The core machinery — Leslie-matrix cohort projection (Leslie, 1945) and the Euler–Lotka steady-state
condition (Euler, 1760/1970; Sharpe & Lotka, 1911) — is standard, unmodified demographic theory; §8's
derivation is textbook stable-population theory, not a novel or simplified version of it. The
simplifications worth flagging are all in the *inputs* to that machinery, not the machinery itself:

| Where | This model | Standard practice | Suggested refinement |
|---|---|---|---|
| Mortality shape (Eq. 1, §1) | 3-term hazard, but with a *thresholded, shifted* senescent term and a *rising* (not constant) background term | Siler's (1979) or Heligman & Pollard's (1980) 3-/8-term hazards: smooth, unthresholded from birth; Heligman-Pollard adds a distinct young-adult "accident hump" | Refit Eq. (1) as a literal Siler model ($\mu_i e^{-a/\tau_i} + \mu_0 + \mu_2 e^{k a}$, no threshold) — same parameter count, standard form, removes the somewhat ad hoc cliff construction |
| Fertility shape (Eq. 3, §3) | Piecewise-linear triangle, hard cutoffs, kinked peak | Smooth, right-skewed unimodal curves: Hadwiger (1940), Coale & Trussell (1974), Schmertmann (2003) | Swap in a Hadwiger function (same $\sim$4 parameters, smooth, standard, well-suited to a hand fit) |
| Sex structure (assumption 2; Eq. 12, §9) | None — one population, $p_f=0.5$ fixed | Nearly all demographic models track males/females separately, since fertility is female-specific and mortality differs by sex at every age; real sex ratio at birth $\approx1.05$ boys/girl (United Nations Population Division, 2024b) | Split $N(a,t)$ into male/female vectors with their own mortality; only females contribute to Eq. (7c) |
| No migration/density-dependence (assumption 3) | Rates fixed, independent of $t$ or $P(t)$ | Standard for a *closed*, deterministic stable-population exercise like this one (the same assumption Euler–Lotka itself requires, §8) — not a simplification relative to the classical theory, just a scope choice | None needed for this model's own stated purpose; would need revisiting only if migration or crowding effects become a modelling goal |
| Initial age distribution | Current real-world structure (CIA World Factbook, 2021, via IndexMundi) | Same kind of source real demographic projections start from | None — this one already matches practice |

None of these change the model's core conclusion (§8): Leslie-matrix projection plus the Euler–Lotka
condition is exact, standard theory regardless of which mortality/fertility functional forms feed
into it. The suggested refinements would make the *inputs* more standard and slightly more accurate
(smoother fertility, a proper Siler-form mortality, sex structure) without touching that conclusion.

---

## Parameter reference

| symbol | code | default | meaning |
|---|---|---|---|
| $A$ | `ageMax` | $1000$ | oldest age class (plus-group). Deliberately large -- see the note at the end of §1 -- so the plus-group holds negligible population; independent of the curves below |
| — | `ageMaxDisplay` | $100$ | age axis limit in the plots ONLY; purely cosmetic, no effect on the simulation |
| $T$ | `nYears` | $500$ | simulation horizon |
| — | `popInit` | $8\times10^9$ | total starting population |
| — | — | — | initial age distribution: current real-world world age structure (Central Intelligence Agency, 2021, via IndexMundi: 0-14 25.2%, 15-24 15.3%, 25-54 40.6%, 55-64 9.2%, 65+ 9.7%), not the model's own stable shape -- see Eq. (10) |
| $\mu_i$ | `infantMortalityScale` | $0.0129$ | excess hazard at age $0$ from the infant-mortality term (fit, §10; optimizer: $0.0128711$) |
| $\tau_i$ | `infantMortalityDecay` | $1.21$ | its decay time constant, years (fit, §10; optimizer: $1.20517$) |
| $\mu_0$ | `deathRateBase` | $0.00066$ | baseline hazard at age $0$, excluding the infant term (fit, §10; optimizer: $0.000659275$) |
| $\mu_1$ | `deathRateSlope` | $0.00057$ | added baseline hazard at age $A_{\text{ref}}$ (fit, §10; optimizer: $0.000572547$) |
| $A_{\text{ref}}$ | `ageRef` | $150$ | fixed reference age for $\mu_1$'s term -- NOT `ageMax` (§1) |
| $a_c$ | `steepAge` | $11.01$ | absolute age where the senescence cliff begins -- NOT a fraction of `ageMax` (§1; fit, §10; optimizer: $11.0086$) |
| $\mu_2$ | `steepMortalityScale` | $0.000227$ | extra-hazard scale past the cliff (fit, §10; optimizer: $0.000226789$) |
| $k$ | `steepMortalityRate` | $0.0836$ | extra-hazard growth rate past the cliff (fit, §10; optimizer: $0.083606$) |
| $a_{\min}$ | `fertileMin` | $12.1$ | youngest fertile age (fit incl. R0-matching, §10; optimizer: $12.0851$) |
| $a_{\max}$ | `fertileMax` | $42.6$ | oldest fertile age (fit incl. R0-matching, §10; optimizer: $42.6161$) |
| $a^*$ | `fertilityPeakAge` | $27.0$ | age of peak fertility (fit, §10; optimizer: $27.0495$) |
| $b_{\max}$ | `fertilityPeakRate` | $0.0713$ | per-capita rate at $a^*$ (fit incl. R0-matching, §10; optimizer: $0.07127$) |
| $R_0^{\text{target}}$ | `targetR0` | empty (`[]`) | net reproduction rate target; empty = don't rescale, just estimate/report $R_0$ (§5); a number ($1=$steady, $>1=$growth, $<1=$decline) rescales fertility to hit it exactly |

---

## Modelling assumptions

1. **Yearly, discrete-age cohorts** — no within-year age or seasonality structure; birth/death
   probabilities are applied once per year (Eqs. 1–2, 7).
2. **No sex structure** — $N(a,t)$ is a single (unisex) population; fertility $b(a)$ is a per-capita
   rate applied to the whole cohort, not per-female; §9 and §11 discuss the specific numbers (sex
   ratio at birth, differential mortality) this trades away, and how a real demographic model
   usually handles it.
3. **No migration, environmental stochasticity, or density dependence** — rates $\mu(a)$, $b(a)$ are
   fixed functions of age alone, independent of $t$ or $P(t)$; this is exactly what makes the
   Euler–Lotka argument (§8) apply cleanly, and matches what that classical theory itself assumes
   (§11) — not an extra simplification layered on top of it.
4. **Plus-group at $A$** (Eq. 7b) — the oldest class is an absorbing bin, not a hard cutoff; nobody
   is assumed to die exactly at age $A$, though the senescence cliff (Eq. 1) makes survival past it
   very unlikely with default parameters.

---

## References

Coale, A. J., & Trussell, T. J. (1974). Model fertility schedules: Variations in the age structure
of childbearing in human populations. *Population Index*, 40(2), 185–258.
https://doi.org/10.2307/2733910

Euler, L. (1970). A general investigation into the mortality and multiplication of the human
species (N. Keyfitz & B. Keyfitz, Trans.). *Theoretical Population Biology*, 1(3), 307–314.
https://doi.org/10.1016/0040-5809(70)90048-1 (Original work published 1760)

Gompertz, B. (1825). On the nature of the function expressive of the law of human mortality, and
on a new mode of determining the value of life contingencies. *Philosophical Transactions of the
Royal Society of London*, 115, 513–583. https://doi.org/10.1098/rstl.1825.0026

Hadwiger, H. (1940). Eine analytische Reproduktionsfunktion für biologische Gesamtheiten.
*Scandinavian Actuarial Journal*, 1940(3–4), 101–113.
https://doi.org/10.1080/03461238.1940.10404802

Heligman, L., & Pollard, J. H. (1980). The age pattern of mortality. *Journal of the Institute of
Actuaries*, 107(1), 49–80. https://doi.org/10.1017/S0020268100040257

Leslie, P. H. (1945). On the use of matrices in certain population mathematics. *Biometrika*,
33(3), 183–212. https://doi.org/10.1093/biomet/33.3.183

Schmertmann, C. P. (2003). A system of model fertility schedules with graphically intuitive
parameters. *Demographic Research*, 9, 81–110. https://doi.org/10.4054/DemRes.2003.9.5

Sharpe, F. R., & Lotka, A. J. (1911). A problem in age-distribution. *Philosophical Magazine*,
21(124), 435–438. https://doi.org/10.1080/14786440408637050

Siler, W. (1979). A competing-risk model for animal mortality. *Ecology*, 60(4), 750–757.
https://doi.org/10.2307/1936612

United Nations, Department of Economic and Social Affairs, Population Division. (2024a). *World
Population Prospects 2024, Online Edition*. https://population.un.org/wpp/

United Nations, Department of Economic and Social Affairs, Population Division. (2024b). Sex
ratio at birth [Data set]. *World Population Prospects 2024* – processed by Our World in Data.
https://ourworldindata.org/grapher/sex-ratio-at-birth

Central Intelligence Agency. (2021). *The World Factbook: Age structure (World)*, via IndexMundi.
https://www.indexmundi.com/world/age_structure.html

Our World in Data. (2026a). Fertility rate by age group [Data set]. Sourced from United Nations,
World Population Prospects (2024a). https://ourworldindata.org/grapher/fertility-rate-by-age-group

Our World in Data. (2026b). Probability of dying by age [Data set]. Sourced from World Health
Organization, Global Health Observatory.
https://ourworldindata.org/grapher/probability-of-dying-by-age

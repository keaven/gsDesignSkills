# Code Patterns for gsDesign

**Note**: These patterns target gsDesign >= 3.9.0 (development version at
github.com/keaven/gsDesign). Features marked **(dev)** may not be on CRAN yet.

## Table of Contents
1. [One-sided efficacy-only design](#one-sided-efficacy)
2. [Two-sided asymmetric design with non-binding futility](#two-sided-nonbinding)
3. [Harm bounds (test.type 7/8)](#harm-bounds) **(dev)**
4. [Time-to-event design with gsSurv](#time-to-event-gssurv)
5. [Calendar-based timing with gsSurvCalendar](#calendar-timing)
6. [Power computation with gsSurvPower](#gsssurvpower) **(dev)**
7. [Normal endpoint design](#normal-endpoint)
8. [Binomial endpoint design](#binomial-endpoint)
9. [Integer sample size rounding](#integer-rounding)
10. [Bound summary and reporting](#bound-summary)
11. [Plotting designs](#plotting)
12. [Spending functions](#spending-functions)
    - [Two-parameter spending function families](#two-parameter-sf)
    - [Choosing between spending function families](#choosing-sf)
    - [Worked example: boundary-driven spending for PFS](#boundary-driven-sf)
13. [Sequential p-values](#sequential-p-values)
14. [Conditional power and sample size re-estimation](#conditional-power)
15. [Updating bounds at analysis](#updating-bounds)
16. [Multiplicity graph with hGraph](#hgraph)
17. [sfExtremeValue2 for futility bound calibration](#sfextremevalue2)
18. [Selective futility bound testing with testLower](#testlower)
19. [Binomial test statistics with testBinomial and nBinomial](#testbinomial)

---

## One-sided efficacy-only design {#one-sided-efficacy}

The simplest design: one-sided testing with only an efficacy (upper) bound.

```r
library(gsDesign)

# Generic design (sample size as ratio of fixed design)
x <- gsDesign(
  k = 3,              # 3 analyses (2 interim + final)
  test.type = 1,      # 1-sided, efficacy only
  alpha = 0.025,      # 1-sided Type I error
  beta = 0.1,         # Type II error (90% power)
  sfu = sfLDOF,       # Lan-DeMets O'Brien-Fleming spending
  timing = c(0.5, 0.75)  # Information fractions at IA1, IA2
)
x
plot(x)
```

### test.type values

| test.type | Description |
|-----------|-------------|
| 1 | One-sided, efficacy only |
| 2 | Two-sided, symmetric |
| 3 | Two-sided, asymmetric, beta-spending, binding futility |
| 4 | Two-sided, asymmetric, beta-spending, non-binding futility |
| 5 | Two-sided, asymmetric, null-spending, binding futility |
| 6 | Two-sided, asymmetric, null-spending, non-binding futility |
| 7 | Two-sided, asymmetric, binding futility + binding harm **(dev)** |
| 8 | Two-sided, asymmetric, non-binding futility + non-binding harm **(dev)** |

## Two-sided asymmetric design with non-binding futility {#two-sided-nonbinding}

Non-binding futility (test.type = 4 or 6) means the trial can continue past
a futility crossing without inflating Type I error. This is the most common
choice for confirmatory trials.

```r
x <- gsDesign(
  k = 3,
  test.type = 4,         # Non-binding futility, beta-spending
  alpha = 0.025,
  beta = 0.1,
  timing = c(0.5, 0.75),
  sfu = sfLDOF,          # Efficacy: Lan-DeMets O'Brien-Fleming
  sfl = sfHSD,           # Futility: Hwang-Shih-DeCani
  sflpar = -2            # Moderate futility spending
)
x
gsBoundSummary(x)
```

### Choosing futility spending

- `sfHSD` with `sflpar = -2`: moderate futility spending (common default)
- `sfHSD` with `sflpar = 1`: aggressive early futility (Pocock-like)
- `sfHSD` with `sflpar = -4`: conservative futility (O'Brien-Fleming-like)
- `sfLDOF`: Lan-DeMets O'Brien-Fleming futility (very conservative early)
- `sfPower` with `sflpar = 0.5`: square-root spending

## Harm bounds (test.type 7/8) {#harm-bounds}

**(dev)** — requires gsDesign >= 3.9.0 from github.com/keaven/gsDesign.

Test types 7 and 8 add a **harm bound** (upper bound for experimental harm)
alongside the standard efficacy and futility bounds. This is useful when
monitoring for the possibility that the experimental treatment is worse
than control (e.g., OS harm in oncology trials).

### FDA guidance context

The FDA draft guidance *Approaches to Assessment of Overall Survival in Oncology Clinical Trials* (August 2025) recommends:

- All randomized oncology trials should be designed to assess OS, even when OS is not the primary endpoint
- Interim analyses of OS for futility or harm should be included, timed to limit exposure to potentially harmful therapeutics
- Timing of interim and final OS analyses should be **event-driven** rather than time-based
- Pre-specify a threshold (HR ≥ some value) to indicate potential harm, with justification
- Design for sufficient events to rule out harm with pre-specified precision (e.g., 95% CI for OS HR excludes the harm threshold)
- When OS is not the primary endpoint, also pre-specify: information fraction, boundaries (on p-value scale and translated to HR), and whether boundaries are binding
- For multiple endpoints with OS, pre-specify a multiplicity plan controlling study-wise Type I error

`test.type = 8` directly addresses these requirements by adding a non-binding harm bound to the group sequential design.

- `test.type = 7`: binding futility + binding harm
- `test.type = 8`: non-binding futility + non-binding harm (recommended)

### Example: OS harm assessment per FDA guidance

5-year study with 12-month enrollment and annual analyses. Efficacy testing
is skipped at Year 1 (too few events), but harm is monitored from the start.
Low control failure rate (median OS = 48 months) keeps sample size realistic.

```r
os_harm <- gsSurvCalendar(
  test.type = 8,              # Non-binding futility + non-binding harm
  alpha = 0.025,
  beta = 0.2,                 # 80% power
  sfu = sfLDOF,               # Efficacy spending
  sfl = sfHSD, sflpar = -2,   # Futility spending
  sfharm = sfHSD,             # Harm bound spending
  sfharmparam = -2,
  astar = 0.05,               # Harm bound Type I error under H0
  calendarTime = c(12, 24, 36, 48, 60),  # Annual analyses
  spending = "information",
  lambdaC = log(2) / 48,     # Control median OS = 48 months
  hr = 0.6,                   # Target HR
  eta = 0.001,
  gamma = 10,                 # Enrollment rate (patients/month)
  R = 12,                     # 12-month enrollment
  minfup = 48,                # 60 - 12 = 48 months minimum follow-up
  testUpper = c(FALSE, TRUE, TRUE, TRUE, TRUE),   # Skip efficacy at Year 1
  testHarm = c(TRUE, TRUE, TRUE, TRUE, TRUE)      # Harm at all analyses
) |> toInteger()

gsBoundSummary(os_harm, timename = "Month", tdigits = 0)
```

Key design choices:

- **Event-driven timing**: although analyses are scheduled annually,
  spending is based on information fraction (`spending = "information"`)
- **Efficacy skipped at Year 1**: too few events for meaningful efficacy
  assessment; harm monitoring starts immediately
- **`astar = 0.05`**: probability of crossing the harm bound under H0;
  controls the false-positive rate for declaring harm
- **Low failure rate**: median OS = 48 months ensures the trial collects
  enough events to rule out harm with adequate precision

### Selective bound testing

Control which analyses include each bound type with `testUpper`, `testLower`,
and `testHarm` (logical vectors of length `k`).

```r
# 3-analysis design: efficacy at all, futility at IA1-IA2, harm at IA1 only
x <- gsDesign(
  k = 3,
  test.type = 8,
  alpha = 0.025,
  beta = 0.1,
  timing = c(0.5, 0.75),
  sfu = sfLDOF,
  sfl = sfHSD, sflpar = -2,
  sfharm = sfHSD, sfharmparam = -2,
  testUpper = c(TRUE, TRUE, TRUE),
  testLower = c(TRUE, TRUE, FALSE),    # No futility at final
  testHarm = c(TRUE, FALSE, FALSE)     # Harm only at IA1
)
```

## Time-to-event design with gsSurv {#time-to-event-gssurv}

`gsSurv()` combines `nSurv()` (sample size) with `gsDesign()` (boundaries)
for a complete survival trial design. This is the recommended approach.

```r
x <- gsSurv(
  k = 3,                    # Number of analyses
  test.type = 4,            # Non-binding futility
  alpha = 0.025,
  beta = 0.1,
  timing = c(0.5, 0.75),   # Event fraction at interims
  sfu = sfLDOF,             # Efficacy spending
  sfl = sfHSD,              # Futility spending
  sflpar = -2,
  lambdaC = log(2) / 12,   # Control median = 12 months
  hr = 0.7,                 # Hazard ratio under H1
  eta = 0.001,              # Annual dropout rate
  gamma = c(2, 4, 6, 8),   # Enrollment rates by period
  R = c(3, 3, 3, 3),       # Duration of each enrollment period
  T = 36,                   # Total study duration (months)
  minfup = 12               # Minimum follow-up (months)
)
x
gsBoundSummary(x, timename = "Month", tdigits = 1)
```

### Key parameters for nSurv/gsSurv

- `lambdaC`: Control hazard rate. For exponential: `log(2) / median_survival`
- `hr`: Hazard ratio (experimental/control) under H1. HR < 1 means experimental is better
- `hr0`: Null hypothesis HR (default 1 for superiority; set < 1 for non-inferiority)
- `eta`: Dropout hazard rate (both arms unless `etaE` specified)
- `gamma`: Enrollment rates. Scalar for constant, vector for piecewise
- `R`: Duration of each enrollment period matching rows of `gamma`
- `T`: Total study duration. Set `T = NULL` to solve for T given `minfup`
- `minfup`: Minimum follow-up. Set `minfup = NULL` to solve for minfup given `T`

### Three ways to derive sample size and obtain desired power {#three-power-patterns}

There are three mutually exclusive approaches for `nSurv()`/`gsSurv()` to derive the sample size needed to achieve target power (Lachin & Foulkes, 1986). Exactly one of these determines the "free variable" that is adjusted:

**Pattern 1: Vary enrollment rate (recommended)**

Fix `T` and `minfup`. The enrollment rates in `gamma` are scaled proportionally to achieve target power. This is the Lachin-Foulkes method and is generally the most reliable approach.

```r
# gamma is scaled proportionally to achieve 90% power
x <- gsSurv(
  k = 3, test.type = 4, alpha = 0.025, beta = 0.1,
  lambdaC = log(2) / 12, hr = 0.7, eta = 0.001,
  gamma = c(2, 4, 6, 8),   # Relative rates; will be scaled
  R = c(3, 3, 3, 3),       # 12 months enrollment total
  T = 36,                   # Fixed study duration
  minfup = 24,              # Fixed minimum follow-up
  sfu = sfLDOF, sfl = sfHSD, sflpar = -2
)
# T = 36 = sum(R) + minfup = 12 + 24
# Output gamma values are scaled from input to achieve power
```

**Pattern 2: Vary enrollment duration**

Fix `gamma` rates and `minfup`, set `T = NULL`. Solves for how long enrollment must last (adjusting `R` and `T`).

```r
# Enrollment rates are fixed; solve for enrollment duration
x <- gsSurv(
  k = 3, test.type = 4, alpha = 0.025, beta = 0.1,
  lambdaC = log(2) / 12, hr = 0.7, eta = 0.001,
  gamma = 10,               # Fixed enrollment rate (patients/month)
  R = 18,                   # Starting value; will be adjusted
  T = NULL,                 # Solve for study duration
  minfup = 24,              # Fixed minimum follow-up
  sfu = sfLDOF, sfl = sfHSD, sflpar = -2
)
# Output: R and T are computed; gamma stays at 10
```

**Pattern 3: Vary follow-up duration (not recommended)**

Fix `gamma` and `R`, set `minfup = NULL`. Solves for minimum follow-up. **This method often fails** -- it produces an error if the fixed enrollment either over-powers the trial with no follow-up or under-powers it with infinite follow-up.

```r
# Fixed enrollment; solve for follow-up (use with caution)
x <- gsSurv(
  k = 3, test.type = 4, alpha = 0.025, beta = 0.1,
  lambdaC = log(2) / 12, hr = 0.7, eta = 0.001,
  gamma = 10,               # Fixed enrollment rate
  R = 18,                   # Fixed enrollment duration
  T = NULL,                 # Will be computed as sum(R) + minfup
  minfup = NULL,            # Solve for minimum follow-up
  sfu = sfLDOF, sfl = sfHSD, sflpar = -2
)
# May fail with: "Minimum follow-up greater than study duration"
```

### Piecewise exponential failure rates

```r
# Control median changes from 12 months (first 6 months) to 18 months (after)
x <- gsSurv(
  k = 3,
  test.type = 4,
  alpha = 0.025,
  beta = 0.1,
  lambdaC = log(2) / c(12, 18),  # Piecewise hazard rates
  S = 6,                          # Changepoint at 6 months
  hr = 0.7,
  eta = 0.001,
  gamma = 6,
  R = 12,
  T = 36,
  minfup = 12
)
```

### Stratified design

```r
# 2 strata with different control hazard rates
x <- gsSurv(
  k = 3,
  test.type = 4,
  alpha = 0.025,
  beta = 0.1,
  lambdaC = matrix(log(2) / c(10, 14), ncol = 2),  # columns = strata
  hr = 0.7,
  eta = 0.001,
  gamma = matrix(c(4, 6), ncol = 2),  # Different enrollment per stratum
  R = 12,
  T = 36,
  minfup = 12
)
```

## Calendar-based timing with gsSurvCalendar {#calendar-timing}

When analyses are triggered by calendar time rather than event counts.

```r
x <- gsSurvCalendar(
  test.type = 4,
  alpha = 0.025,
  beta = 0.1,
  sfu = sfLDOF,
  sfl = sfHSD,
  sflpar = -2,
  calendarTime = c(18, 26, 36),    # Calendar months from start
  spending = "information",         # or "calendar" for calendar-based spending
  lambdaC = log(2) / 12,
  hr = 0.7,
  eta = 0.001,
  gamma = c(2, 4, 6),
  R = c(3, 3, 6),
  minfup = 12
) %>% toInteger()

gsBoundSummary(x, timename = "Month", tdigits = 1)
```

**Key difference**: `spending = "information"` (default) uses information fraction
for alpha spending; `spending = "calendar"` uses calendar time fraction, which
spends less at early interims and saves more alpha for the final analysis.

## Power computation with gsSurvPower {#gsssurvpower}

**(dev)** — requires gsDesign >= 3.9.0 from github.com/keaven/gsDesign.

`gsSurvPower()` computes power for a group sequential survival design under
specified assumptions, without solving for sample size. Unlike `gsSurv()` which
solves for enrollment to achieve target power, `gsSurvPower()` takes fixed
assumptions and returns the resulting power.

### Quick start: reuse a design object

```r
design <- gsSurv(
  k = 3, test.type = 4, alpha = 0.025, sided = 1, beta = 0.1,
  sfu = sfHSD, sfupar = -4, sfl = sfHSD, sflpar = -2,
  lambdaC = log(2) / 12, hr = 0.7, eta = 0.01,
  gamma = 10, R = 16, minfup = 12, T = 28
)

# Power at the design HR (should reproduce ~90%)
pwr <- gsSurvPower(x = design, plannedCalendarTime = design$T)
pwr$power

# What-if: power under a worse HR
gsSurvPower(x = design, hr = 0.8, plannedCalendarTime = design$T)$power

# What-if: different alpha (e.g., from multiplicity reallocation)
gsSurvPower(x = design, alpha = 0.0125, plannedCalendarTime = design$T)$power
```

### Calendar time vs event-driven timing

```r
# Calendar-driven: fix analysis times, events recomputed under assumed HR
# A worse HR -> more events at same time -> different information fractions
pwr_cal <- gsSurvPower(x = design, hr = 0.8, plannedCalendarTime = design$T)

# Event-driven: fix event counts, calendar times adjust
# Information fractions stay constant -> matches gsDesign power plot
pwr_evt <- gsSurvPower(x = design, hr = 0.8, targetEvents = design$n.I)
```

### Complex analysis timing rules

Combine multiple criteria per analysis: calendar time floors, event targets,
maximum extensions, minimum time between analyses, and minimum follow-up.

```r
gsSurvPower(
  x = design,
  hr = 0.75,
  plannedCalendarTime = c(20, 28, NA),    # Floor time for IA1, IA2; final is event-driven
  targetEvents = c(NA, NA, 300),           # Target 300 events at final
  maxExtension = c(NA, NA, 6),             # Max 6-month extension for final
  minTimeFromPreviousAnalysis = c(NA, 6, 6),  # At least 6 months between analyses
  minFollowUp = c(6, NA, 12)              # Min 6-mo FU at IA1, 12-mo at final
)
```

### Hazard ratio roles: hr vs hr1

- `hr`: The assumed treatment effect for power computation (the "what-if")
- `hr1`: The design alternative that calibrates futility bounds (test.type 3, 4, 7, 8)

When `x` is provided, `hr1` defaults to `x$hr`, so futility bounds stay calibrated
to the original design even when evaluating power at a different `hr`.

```r
# Futility calibrated to HR = 0.7 (from design), power evaluated at HR = 0.8
gsSurvPower(x = design, hr = 0.8, plannedCalendarTime = design$T)$power

# Override futility calibration (unusual but possible)
gsSurvPower(x = design, hr = 0.8, hr1 = 0.8, plannedCalendarTime = design$T)$power
```

### Without a reference design

```r
gsSurvPower(
  k = 2, test.type = 4, alpha = 0.025, sided = 1,
  sfu = sfLDOF, sfl = sfHSD, sflpar = -2,
  lambdaC = log(2) / 6, hr = 0.65, eta = 0.01,
  gamma = 8, R = 18, ratio = 1,
  plannedCalendarTime = c(24, 36)
)$power
```

### targetN for enrollment scaling

Use `targetN` to rescale enrollment to hit a target sample size without
changing the relative enrollment ramp-up.

```r
# Original design enrolls ~160 patients
# What if we can only enroll 120?
gsSurvPower(x = design, targetN = 120, plannedCalendarTime = design$T)$power
```

### informationRates for spending control

Cap spending at planned information fractions to prevent over-spending when
events accrue faster than planned.

```r
gsSurvPower(
  x = design,
  hr = 0.8,
  plannedCalendarTime = design$T,
  informationRates = design$timing,   # Cap spending at planned fractions
  fullSpendingAtFinal = TRUE           # Ensure full alpha spent at final
)$power
```

### Harm bounds with gsSurvPower

```r
gsSurvPower(
  k = 3, test.type = 8,
  alpha = 0.025, sided = 1,
  sfu = sfLDOF, sfl = sfHSD, sflpar = -2,
  sfharm = sfHSD, sfharmparam = -2,
  lambdaC = log(2) / 12, hr = 0.7, eta = 0.01,
  gamma = 10, R = 16, minfup = 12,
  plannedCalendarTime = c(20, 28, 36),
  testHarm = c(TRUE, TRUE, FALSE)  # Harm bound at IA1 and IA2 only
)$power
```

### Variance method options

```r
# Default: Lachin-Foulkes (recommended)
gsSurvPower(x = design, method = "LachinFoulkes", plannedCalendarTime = design$T)$power

# Schoenfeld (simpler, slightly conservative)
gsSurvPower(x = design, method = "Schoenfeld", plannedCalendarTime = design$T)$power

# Freedman (accounts for non-proportional allocation)
gsSurvPower(x = design, method = "Freedman", plannedCalendarTime = design$T)$power
```

## Normal endpoint design {#normal-endpoint}

For continuous endpoints with known variance.

```r
# Fixed design sample size
n <- nNormal(
  delta1 = 0.5,     # Treatment difference under H1
  sd = 1.1,         # Standard deviation
  alpha = 0.025,
  beta = 0.1,
  ratio = 1          # 1:1 randomization
)

# Group sequential design
x <- gsDesign(
  k = 3,
  test.type = 4,
  n.fix = n,
  alpha = 0.025,
  beta = 0.1,
  timing = c(0.5, 0.75),
  sfu = sfLDOF,
  sfl = sfHSD,
  sflpar = -1,
  delta1 = 0.5,       # Natural parameter for bound display
  endpoint = "normal"
)
gsBoundSummary(x)
```

## Binomial endpoint design {#binomial-endpoint}

For comparing two proportions.

```r
library(gsDesign)

# Fixed design sample size (risk difference scale)
n.fix <- nBinomial(p1 = 0.15, p2 = 0.10, alpha = 0.025, beta = 0.1)

# Group sequential design
x <- gsDesign(
  k = 2,
  test.type = 4,
  n.fix = n.fix,
  alpha = 0.025,
  beta = 0.1,
  delta1 = 0.05,        # p1 - p2 under H1
  endpoint = "Binomial"
)
gsBoundSummary(x, deltaname = "p[C]-p[E]")

# Risk ratio scale
n.fix.rr <- nBinomial(p1 = 0.30, p2 = 0.15, scale = "RR")
x.rr <- gsDesign(
  k = 2,
  n.fix = n.fix.rr,
  delta1 = log(0.15 / 0.30),
  endpoint = "Binomial"
)
gsBoundSummary(x.rr, deltaname = "RR", logdelta = TRUE)
```

### Power table for binomial endpoint

```r
# Compute power across a grid of control rates and treatment effects
power_table <- binomialPowerTable(
  pC = c(0.8, 0.9),
  delta = seq(-0.05, 0.05, 0.025),
  n = 70,
  ratio = 1,
  alpha = 0.025
)
```

## Integer sample size rounding {#integer-rounding}

Always round to integers before reporting a design.

```r
# For survival designs (rounds events and sample sizes)
x <- gsSurv(...) %>% toInteger()

# With unequal randomization (rounds to multiples of ratio + 1)
x <- gsSurv(..., ratio = 2) %>% toInteger(ratio = 2)

# Round down at final analysis instead of up
x <- gsSurv(...) %>% toInteger(roundUpFinal = FALSE)
```

## Bound summary and reporting {#bound-summary}

```r
# Standard bound summary
gsBoundSummary(x)

# For survival designs, show time units
gsBoundSummary(x, timename = "Month", tdigits = 1)

# Show hazard ratio at bounds (log scale by default for survival)
gsBoundSummary(x, logdelta = TRUE, deltaname = "HR")

# Exclude certain statistics from summary
gsBoundSummary(x, exclude = c("B-value", "CP", "CP H1", "PP"))

# Show all available boundary descriptions
gsBoundSummary(x, exclude = NULL)

# One-line summary
summary(x)

# Output to gt table
gsBoundSummary(x) %>% as_gt()

# Output to RTF
gsBoundSummary(x) %>% as_rtf(file = "bounds.rtf")

# Output to LaTeX
xprint(xtable::xtable(gsBoundSummary(x), caption = "Boundary summary"))
```

## Plotting designs {#plotting}

`plot.gsDesign()` supports 7 plot types via `plottype`:

```r
library(ggplot2)

# 1: Z-value boundaries (default)
plot(x)
plot(x, plottype = "Z")

# 2: Power by effect size
plot(x, plottype = "power")

# 3: Approximate treatment effect at boundaries
plot(x, plottype = "thetahat")

# 4: Conditional power at boundaries
plot(x, plottype = "CP")

# 5: Spending functions
plot(x, plottype = "sf")

# 6: Expected sample size (ASN)
plot(x, plottype = "ASN")

# 7: B-values at boundaries
plot(x, plottype = "B-value")

# Hazard ratio plot (for survival designs)
plot(x, plottype = "HR")

# Risk ratio plot (for binomial designs)
plot(x, plottype = "RR")
```

## Spending functions {#spending-functions}

All spending functions take the same signature: `sf(alpha, t, param)` where
`alpha` is the total spending, `t` is the information fraction (0 to 1),
and `param` is the function-specific parameter.

### Common spending functions

```r
# Lan-DeMets O'Brien-Fleming (no parameter needed)
# Conservative early, aggressive late
sfLDOF(alpha = 0.025, t = c(0.5, 0.75, 1))

# Hwang-Shih-DeCani (gamma parameter)
# gamma < 0: more conservative (OBF-like); gamma > 0: more aggressive (Pocock-like)
sfHSD(alpha = 0.025, t = c(0.5, 0.75, 1), param = -4)  # Conservative
sfHSD(alpha = 0.025, t = c(0.5, 0.75, 1), param = 1)   # Aggressive

# Kim-DeMets power spending
# param = power exponent: 1 = uniform, 2 = quadratic, 3 = cubic
sfPower(alpha = 0.025, t = c(0.5, 0.75, 1), param = 3)

# Pointwise spending (custom knots)
# param = cumulative spending at each time point
sfPoints(alpha = 0.025, t = c(0.5, 0.75, 1), param = c(0.005, 0.015, 0.025))

# Generalized Lan-DeMets O'Brien-Fleming
# param > 0 controls degree (larger = more conservative)
sfLDOF(alpha = 0.025, t = c(0.5, 0.75, 1), param = 2)
```

### Comparing spending functions visually

```r
t <- seq(0, 1, 0.01)

# Compare multiple spending functions
par(mfrow = c(1, 1))
plot(t, sfLDOF(0.025, t)$spend, type = "l", ylab = "Cumulative spending",
     xlab = "Information fraction", main = "Spending function comparison")
lines(t, sfHSD(0.025, t, -4)$spend, col = 2)
lines(t, sfHSD(0.025, t, -2)$spend, col = 3)
lines(t, sfHSD(0.025, t, 1)$spend, col = 4)
legend("topleft",
  legend = c("sfLDOF", "sfHSD(-4)", "sfHSD(-2)", "sfHSD(1)"),
  col = 1:4, lty = 1
)
```

### Two-parameter spending function families {#two-parameter-sf}

**Reference**: Anderson KM, Clark JB. Fitting spending functions.
*Statist. Med.* 2010; 29:321–327.

Two-parameter families use the general form (equation 12 in the paper):

$$\alpha(t; a, b) = \alpha F(a + b F^{-1}(t))$$

where $F$ is a CDF defined on $(-\infty, \infty)$. Given two desired spending
points, the parameters $a$ and $b$ are found by solving a linear system:

$$F^{-1}(s_i) = a + b \times F^{-1}(t_i), \quad i = 0, 1$$

#### Parameterization modes

All two-parameter functions accept two calling conventions:

1. **Direct**: `param = c(a, b)` — specify internal parameters directly
2. **Point-fitting**: `param = c(t1, t2, u1, u2)` — specify two points where
   `sf(t1) = alpha * u1` and `sf(t2) = alpha * u2`. gsDesign solves for `a, b`
   internally. All four values must be in (0, 1) with `t1 < t2` and `u1 < u2`.

```r
# Point-fitting mode (most common): fit spending at two information fractions
# At IF = 0.1, spend 10% of alpha; at IF = 0.4, spend 40% of alpha
sfLogistic(alpha = 0.025, t = c(0.25, 0.5, 0.75, 1),
           param = c(0.1, 0.4, 0.10, 0.40))

# Direct mode: specify a and b directly
sfLogistic(alpha = 0.025, t = c(0.25, 0.5, 0.75, 1),
           param = c(0, 1))
```

#### Distribution choices

| Function | CDF $F$ | Properties |
|----------|---------|------------|
| `sfLogistic` | Logistic | General purpose; used in GUSTO V trial |
| `sfNormal` | Standard normal | Nearly identical to sfLogistic |
| `sfCauchy` | Cauchy | Flat between fitted points; **robust to timing changes** |
| `sfExtremeValue` | Extreme value: $F(x)=\exp(-e^{-x})$ | Conservative early spending |
| `sfExtremeValue2` | Flipped extreme value: $F(x)=1-\exp(-\exp(x))$ | Conservative early; **(dev)** |
| `sfTDist` | t-distribution (3 params: $a$, $b$, df) | df=1 → Cauchy; df=∞ → Normal |
| `sfBetaDist` | Incomplete beta ($a>0$, $b>0$) | Nonlinear fitting via `nlminb()` |

#### Code examples for each family

```r
t <- seq(0, 1, 0.01)
# All fitted to the same two points: (0.1, 0.10) and (0.4, 0.40)
pts <- c(0.1, 0.4, 0.10, 0.40)

# Logistic
sfLogistic(0.025, t, pts)$spend

# Normal
sfNormal(0.025, t, pts)$spend

# Cauchy (uses Cauchy CDF; flat between fitted points)
sfCauchy(0.025, t, pts)$spend

# Extreme value
sfExtremeValue(0.025, t, pts)$spend

# Flipped extreme value (dev)
sfExtremeValue2(0.025, t, pts)$spend

# t-distribution with 3 df (between Cauchy and Normal)
sfTDist(0.025, t, c(pts, 3))$spend   # 5th param = df

# Beta distribution (different parameterization: a, b > 0)
sfBetaDist(0.025, t, c(0.5, 2))$spend
```

#### Comparing two-parameter families visually

Adapted from Figure 1 of Anderson & Clark (2009):

```r
# All spending functions fit to same point at t=0.75: spend = alpha * 0.75^3
t <- seq(0, 1, 0.01)
alpha <- 0.025
pts <- c(0.4, 0.75, 0.004 / alpha, 0.75^3)

plot(t, sfCauchy(alpha, t, pts)$spend, type = "l", lty = 2,
     xlab = "Information fraction", ylab = "Cumulative spending",
     main = "Two-parameter spending families")
lines(t, sfLogistic(alpha, t, pts)$spend, lty = 1)
lines(t, sfNormal(alpha, t, pts)$spend, lty = 3)
lines(t, sfExtremeValue(alpha, t, pts)$spend, lty = 4)
legend("topleft",
  legend = c("Cauchy", "Logistic", "Normal", "Extreme value"),
  lty = c(2, 1, 3, 4))
```

### Choosing between spending function families {#choosing-sf}

#### One-parameter vs two-parameter

Use **one-parameter** when a single desired property suffices (e.g., "conservative
early spending like O'Brien-Fleming"). Use **two-parameter** when you need to
match desired boundaries at two specific interim analyses — for example, a
proof-of-concept futility bar at IA1 and an accelerated approval efficacy bar
at IA2.

#### Selection guide

| If you need... | Use | Why |
|---------------|-----|-----|
| Conservative early, default | `sfLDOF` | Standard O'Brien-Fleming approximation |
| Tunable conservatism (1 param) | `sfHSD` with $\gamma < 0$ | Single parameter controls shape |
| Most conservative early (1 param) | `sfExponential` with $\nu \approx 0.8$ | More conservative than sfHSD or sfPower at same late spending |
| Fit two boundary values | `sfLogistic` or `sfNormal` | General purpose; nearly identical |
| Robust to timing changes | `sfCauchy` | Flat between fitted points; spending barely changes if analysis moves |
| Conservative early + fit two | `sfExtremeValue` or `sfExtremeValue2` | Steep near fitted points; spends little early |
| Continuous Cauchy ↔ Normal | `sfTDist` with df | df=1→ Cauchy; df=∞→ Normal; intermediate df fills the gap |
| Maximum flexibility (2 param) | `sfBetaDist` | Incomplete beta; different shape family than F-inverse class |

#### Key insight from Anderson & Clark (2009)

> "Choosing other values for $s_0$ produces logistic spending functions that
> closely approximate the Hwang–Shih–DeCani and power family functions...
> Approximating the exponential function with a logistic spending function
> produces a substantially nonconvex spending function."

In other words: logistic/normal can mimic sfHSD or sfPower, but not sfExponential.
If conservative early spending is essential, use sfExponential (one-param) or
sfExtremeValue (two-param) directly.

### Worked example: boundary-driven spending for PFS {#boundary-driven-sf}

Adapted from Section 4.2 of Anderson & Clark (2009). Start from desired
Z-value boundaries, derive spending, then fit a two-parameter function.

**Setup**: PFS endpoint with 10% of total alpha (α = 0.0025). Three analyses
at 100, 400, and 900 PFS events ($t$ = 0.111, 0.444, 1.0). Control median PFS
= 6 months; design HR = 0.6.

**Step 1**: Choose desired boundary properties at each interim.

| Analysis | Events | $t$ | Desired property | $Z$ bound | HR at bound | α-spending fraction |
|----------|--------|-----|-----------------|-----------|-------------|-------------------|
| IA1 | 100 | 0.111 | Very high bar | 3.48 | 0.50 | 0.10 |
| IA2 | 400 | 0.444 | Accelerated approval | 2.87 | 0.75 | 0.90 |
| Final | 900 | 1 | — | 3.39 | 0.80 | 1.00 |

**Step 2**: HR at boundary ≈ $\exp(2Z/\sqrt{d})$ where $d$ = events.

**Step 3**: Fit a two-parameter spending function using the IA1 and IA2 points.

```r
# Alpha for PFS = 0.0025 (10% of 0.025)
# Desired: sf(0.111) = 0.0025 * 0.10, sf(0.444) = 0.0025 * 0.90
pts <- c(0.111, 0.444, 0.10, 0.90)

# Compare which family best suits the need
t <- seq(0, 1, 0.01)
plot(t, sfLogistic(0.0025, t, pts)$spend, type = "l",
     xlab = "Information fraction", ylab = "Cumulative spending",
     main = "Fitted alpha-spending for PFS")
lines(t, sfCauchy(0.0025, t, pts)$spend, lty = 2)
lines(t, sfExtremeValue(0.0025, t, pts)$spend, lty = 4)
legend("topleft", legend = c("Logistic", "Cauchy", "Extreme value"),
       lty = c(1, 2, 4))

# Use in gsDesign (logistic example)
x <- gsDesign(
  k = 3,
  test.type = 4,
  alpha = 0.0025,
  beta = 0.1,
  timing = c(100/900, 400/900),
  sfu = sfLogistic,
  sfupar = pts,           # c(t1, t2, u1, u2)
  sfl = sfHSD,
  sflpar = -2
)
gsBoundSummary(x)
```

**Choosing between families**: If the second interim might shift in timing,
prefer `sfCauchy` (flat → spending barely changes). If an unplanned early
interim is possible and minimal early spending is desired, prefer
`sfExtremeValue` (steep early → spends less if an early look occurs).

## Sequential p-values {#sequential-p-values}

Sequential p-values are the minimum significance level at which a group
sequential bound would be rejected at or before the current analysis.
Used with graphical multiplicity methods (Maurer & Bretz, 2013).

**Reference**: Liu Q, Anderson KM. On adaptive extensions of group sequential
trials for clinical investigations. *JASA* 2008; 103:1621–1630.

### Formal definition

Given a class of well-ordered boundaries $b_k(\mu)$ satisfying
$P_{\Delta=0}\{Z_1 \geq b_1(\mu) \cup \cdots \cup Z_K \geq b_K(\mu)\} = \mu$,
the sequential p-value at analysis $k$ is:

$$p_k = \sup\{\mu: \max_{1 \leq i \leq k}\{Z_i - b_i(\mu)\} \leq 0\}$$

**Interpretation**: "raise the boundary until it just touches the observed data."

### Key properties (Liu & Anderson Theorem 2)

- $p_1 \geq p_2 \geq \cdots \geq p_K$ (non-increasing as evidence accumulates)
- $P_{\Delta=0}\{p_\tau \leq \mu\} \leq \mu$ for any stopping time $\tau$ and $\mu \in (0,1)$
- $p_\tau \leq \alpha$ if and only if the significance boundary was crossed at or before $\tau$
- The final p-value adheres to the ITT principle: all available data are analyzed

### Why it works: Theorem 1 (EGS test validity)

For significance boundaries $b_k$ satisfying the Type I error equation at level
$\alpha$, the probability of rejection is $\leq \alpha$ for **any** stopping time $\tau$
— not just the classical boundary-crossing time. This means:

1. **Non-binding futility**: the trial can continue past a futility crossing
   without inflating Type I error (this is why `test.type = 4` works)
2. **Trial extension**: the trial can be extended past an efficacy boundary
   crossing (e.g., for safety or secondary endpoints) while maintaining
   error control
3. **Flexible stopping**: DSMB decisions based on the totality of data (not
   just the primary endpoint Z-statistic) are valid

### Sequential confidence intervals

From the same theory, sequential confidence bounds at analysis $k$:

$$\hat\Delta^L_k = \max_{1 \leq i \leq k}\{(Z_i - b_i) / \mathcal{I}_i^{1/2}\}, \quad
\hat\Delta^U_k = \min_{1 \leq i \leq k}\{(Z_i + b_i) / \mathcal{I}_i^{1/2}\}$$

The interval $(\hat\Delta^L_\tau, \hat\Delta^U_\tau)$ has coverage $\geq 1 - 2\alpha$
for any stopping time $\tau$.

### Using sequential p-values with multiplicity testing

Sequential p-values can be passed to graphical multiplicity procedures
(e.g., `graphicalMCP::graph_test_shortcut()`) at each analysis. This
implements the Maurer-Bretz (2013) framework for testing multiple hypotheses
in group sequential trials while controlling FWER. See the
graphicalMCP-gsDesign2 and multi-endpoint-sim skills.

**Warning** (Hung, Wang & O'Neill, 2007): Naively testing a secondary
hypothesis at level $\alpha$ whenever the primary is significant does NOT
control FWER in the group sequential setting. You must use sequential
p-values with a proper closed testing or graphical procedure.

### Computing sequential p-values with gsDesign

```r
# Design the trial
x <- gsSurv(
  k = 4, alpha = 0.025, beta = 0.1,
  timing = c(0.5, 0.65, 0.8),
  sfu = sfLDOF,
  lambdaC = log(2) / 6, hr = 0.6,
  eta = 0.01,
  gamma = c(2.5, 5, 7.5, 10), R = c(2, 2, 2, 6),
  T = 30, minfup = 18
)

# Compute sequential p-value from observed data
# n.I = observed event counts, Z = observed Z-statistics
sequentialPValue(
  gsD = x,
  n.I = c(100, 160, 190, 230),   # Observed events at each analysis
  Z = c(1.5, 2, 2.5, 3),         # Observed Z-values (positive = favorable)
  usTime = x$timing               # Use planned spending time
)

# At interim (only first 2 analyses observed)
sequentialPValue(
  gsD = x,
  n.I = c(100, 160),
  Z = c(1.5, 2)
)
```

### Important notes on sequentialPValue

- Works with `test.type = 1, 4, 6` (one-sided or non-binding futility)
- Not meaningful for `test.type = 2, 3, 5` (binding futility or symmetric)
- Z-values must be positive when the experimental treatment is favorable
- `usTime`: if NULL, spending is based on observed information fraction `n.I / max(n.I)`
- Default search interval is `c(1e-05, 0.9999)`; adjust for very small p-values
- The theoretical validity relies on Theorem 1 of Liu & Anderson (2008):
  Type I error ≤ α for any stopping time τ, regardless of why the trial stopped

## Conditional power and sample size re-estimation {#conditional-power}

### Conditional power at interim

```r
# Conditional power at the interim bound values
gsBoundCP(x)                      # At observed treatment effect (thetahat)
gsBoundCP(x, theta = x$delta)     # At the design alternative

# Conditional power at a specific z-value and analysis
gsCP(x, i = 1, z = 2.0)                     # At observed effect
gsCP(x, i = 1, z = 2.0, theta = x$delta)    # At design alternative
```

### Sample size re-estimation with ssrCP

```r
# Start with a 2-stage group sequential design
x <- gsDesign(
  k = 2, n.fix = 100, timing = 0.5,
  alpha = 0.025, beta = 0.1,
  test.type = 4,
  sfu = sfHSD, sfupar = -12,
  delta1 = 1
)

# Derive conditional power design
z <- seq(0, 4, 0.025)
xx <- ssrCP(
  x = x,
  z1 = z,
  overrun = 25,          # Patients enrolled but not yet in interim analysis
  beta = 0.1,            # Target conditional power = 90%
  cpadj = c(0.3, 0.9),  # CP range for SSR (below 0.3: stop; above 0.9: keep original)
  maxinc = 2,            # Max 2x sample size increase
  z2 = z2NC              # Inverse normal combination test
)
plot(xx)

# Compute unconditional power
Power.ssrCP(x = xx)
```

### Combination test options for ssrCP

- `z2NC`: Inverse normal combination test (Lehmacher & Wassmer, 1999)
- `z2Z`: Sufficient statistic (uses all data; only valid in specific settings)
- `z2Fisher`: Fisher combination test

## Updating bounds at analysis {#updating-bounds}

When observed timing differs from planned, update the design bounds.

```r
# Original design
x <- gsDesign(k = 3, test.type = 4, alpha = 0.025, beta = 0.1,
              n.fix = 800, sfu = sfLDOF, sfl = sfHSD, sflpar = -2)

# Actual analyses occurred at different sample sizes
y <- gsDesign(
  k = 3,
  test.type = 4,
  alpha = 0.025,
  n.fix = 800,
  n.I = c(300, 600, 860),          # Actual sample sizes
  maxn.IPlan = x$n.I[x$k],         # Planned max from original design
  sfu = sfLDOF,
  sfl = sfHSD,
  sflpar = -2
)
y

# For spending time different from information fraction
y <- gsDesign(
  k = 3,
  test.type = 1,
  alpha = 0.025,
  n.fix = 800,
  n.I = c(300, 600, 860),
  maxn.IPlan = x$n.I[x$k],
  usTime = c(0.35, 0.7, 1.0),      # Custom spending times
  sfu = sfLDOF
)
```

## Multiplicity graph with hGraph {#hgraph}

`hGraph()` creates a ggplot2-based multiplicity graph visualization.

```r
library(gsDesign)

# Simple 4-hypothesis graph
hGraph(
  nHypotheses = 4,
  nameHypotheses = c("H1: OS Sub", "H2: OS All", "H3: PFS Sub", "H4: PFS All"),
  alphaHypotheses = c(0.01, 0.01, 0.005, 0),
  m = matrix(c(
    0, 1, 0, 0,
    0, 0, 0.5, 0.5,
    0, 0, 0, 1,
    1, 0, 0, 0
  ), nrow = 4, byrow = TRUE),
  fill = c(1, 1, 2, 2),     # Color groups
  palette = c("#4472C4", "#ED7D31"),
  halfWid = 1.5,
  halfHgt = 0.5,
  trhw = 0.15,
  trhh = 0.075,
  digits = 4,
  trdigits = 3,
  size = 5,
  boxtextsize = 4,
  trprop = 0.35
)
```

## sfExtremeValue2 for futility bound calibration {#sfextremevalue2}

**(dev)** — requires gsDesign >= 3.9.0 from github.com/keaven/gsDesign.

`sfExtremeValue2` is a two-parameter spending function using the flipped
extreme value CDF $F(x) = 1 - \exp(-\exp(x))$ (see [Two-parameter spending
function families](#two-parameter-sf) for the general theory from Anderson
& Clark, 2009). It is particularly useful for calibrating futility bounds
to target a specific hazard ratio at the futility boundary.

The 4-parameter form `c(t1, t2, u1, u2)` specifies two points on the spending
curve: `sf(t1) = beta * u1` and `sf(t2) = beta * u2`, where `t1, t2` are
information fractions and `u1, u2` are proportions of total beta spent.

### Usage with gsSurvCalendar

```r
# OS design: sfExtremeValue2 futility targeting ~HR = 0.9 at bound
os_design <- gsSurvCalendar(
  test.type = 4,
  alpha = 0.02,
  beta = 0.1,
  sfu = sfLDOF,
  sfl = sfExtremeValue2,
  sflpar = c(0.37, 0.9, 0.75, 0.95),
  calendarTime = c(18, 24, 30, 42),
  spending = "information",
  lambdaC = log(2) / 14,
  hr = 0.75,
  eta = 0.001,
  gamma = c(0.25, 0.50, 0.75, 1.00),
  R = c(2, 2, 2, 12),
  minfup = 24
) |> toInteger()
```

### Parameterization

The 4 parameters `c(t1, t2, u1, u2)` define the spending curve (consistent
with the general two-parameter family parameterization):

- `t1`, `t2`: Information fractions (0 < t1 < t2 ≤ 1)
- `u1`, `u2`: Proportions of total beta spent at those fractions (0 < u1 < u2 ≤ 1)

To calibrate a futility bound to a target HR:

1. Determine the information fraction at each analysis from the design
2. Choose `t0` near the IF at IA1, `t1` near the IF at IA2
3. Adjust `p0` and `p1` to move the futility bound HR toward the target
4. Higher `p0`/`p1` → more aggressive futility spending → futility HR closer to 1.0

### Calibration example for PFS

```r
# PFS with IF ≈ 0.6 at IA1; target futility HR ~ 0.9
# Use small difference between p0 and p1 for gradual spending
pfs_power <- gsSurvPower(
  k = 3,
  test.type = 4,
  alpha = 0.004,
  sided = 1,
  sfu = sfLDOF,
  sfl = sfExtremeValue2,
  sflpar = c(0.6, 0.9, 0.15, 0.16),   # x=0.15 gives HR~0.9 at IA1
  lambdaC = log(2) / c(8, 16),
  S = 8,
  hr = 0.75,
  eta = 0.008,
  gamma = os_design$gamma,
  R = c(2, 2, 2, 12),
  ratio = 1,
  plannedCalendarTime = c(18, 24, 30)
)
```

## Selective futility bound testing with testLower {#testlower}

`testLower` is a logical vector of length `k` controlling which analyses
include a futility (lower) bound. Set `FALSE` at analyses where no
futility stopping is desired (e.g., at the final analysis or late interims
where efficacy testing dominates).

### Pattern: early futility only

```r
# 4-analysis OS design: futility at IA1 and IA2 only
os_design <- gsSurvCalendar(
  test.type = 4,
  alpha = 0.02,
  beta = 0.1,
  sfu = sfLDOF,
  sfl = sfExtremeValue2,
  sflpar = c(0.37, 0.9, 0.75, 0.95),
  testLower = c(TRUE, TRUE, FALSE, FALSE),   # Futility at IA1, IA2 only
  calendarTime = c(18, 24, 30, 42),
  spending = "information",
  lambdaC = log(2) / 14,
  hr = 0.75,
  eta = 0.001,
  gamma = c(0.25, 0.50, 0.75, 1.00),
  R = c(2, 2, 2, 12),
  minfup = 24
) |> toInteger()

# 3-analysis PFS design: futility at IA1 only
pfs_power <- gsSurvPower(
  k = 3,
  test.type = 4,
  sfl = sfExtremeValue2,
  sflpar = c(0.6, 0.9, 0.15, 0.16),
  testLower = c(TRUE, FALSE, FALSE),   # Futility at IA1 only
  ...
)
```

### Important: toInteger preserves testLower

As of gsDesign >= 3.9.0.9004 (branch `242-extremez`), `toInteger()` preserves
the `testLower` setting. Earlier versions reset suppressed futility bounds
to real values. Verify your gsDesign version:

```r
packageVersion("gsDesign")
# Must be >= 3.9.0.9004 for testLower + toInteger to work correctly
```

## Binomial test statistics with testBinomial and nBinomial {#testbinomial}

`testBinomial()` computes a Z-statistic for comparing two binomial proportions.
`nBinomial()` computes sample size or power for binomial comparisons.

### Sample size

```r
# Fixed design sample size for ORR: control 15%, experimental 25%
n_orr <- nBinomial(
  p1 = 0.15,     # Control proportion
  p2 = 0.25,     # Experimental proportion
  alpha = 0.001,
  beta = 0.1
)

# Power at a given sample size
beta_orr <- nBinomial(
  p1 = 0.15,
  p2 = 0.25,
  n = 500,
  alpha = 0.001
)
power_orr <- 1 - beta_orr
```

### Z-statistic (sign convention)

**Critical**: `testBinomial(x1, x2, n1, n2)` returns a positive Z when
`x1/n1 > x2/n2`. To get **positive Z = experimental better**, pass
experimental as `x1` and control as `x2`:

```r
# Positive Z = experimental better (experimental in x1 position)
z <- testBinomial(
  x1 = x_exp,    # Experimental successes
  x2 = x_ctrl,   # Control successes
  n1 = n_exp,     # Experimental sample size
  n2 = n_ctrl     # Control sample size
)
```

**Common mistake**: Passing `x1 = ctrl, x2 = exp` produces a *negative* Z when
the experimental arm is better. This causes downstream sign issues with
sequential p-values and graphical testing.

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
13. [Sequential p-values](#sequential-p-values)
14. [Conditional power and sample size re-estimation](#conditional-power)
15. [Updating bounds at analysis](#updating-bounds)
16. [Multiplicity graph with hGraph](#hgraph)

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

- `test.type = 7`: binding futility + binding harm
- `test.type = 8`: non-binding futility + non-binding harm (recommended)

```r
x <- gsDesign(
  k = 3,
  test.type = 8,           # Non-binding futility + non-binding harm
  alpha = 0.025,
  beta = 0.1,
  timing = c(0.5, 0.75),
  sfu = sfLDOF,            # Efficacy spending
  sfl = sfHSD,             # Futility spending
  sflpar = -2,
  sfharm = sfHSD,          # Harm bound spending
  sfharmparam = -2,        # Harm bound parameter
  astar = 0.05             # Total harm bound crossing probability under H0
)
gsBoundSummary(x)
```

### Selective bound testing

Control which analyses include each bound type with `testUpper`, `testLower`,
and `testHarm` (logical vectors of length `k`).

```r
# Efficacy at all analyses, futility only at IA1 and IA2, harm only at IA1
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

## Sequential p-values {#sequential-p-values}

Sequential p-values are the minimum alpha at which a group sequential bound
would be rejected. Used with graphical multiplicity methods (Maurer-Bretz).

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

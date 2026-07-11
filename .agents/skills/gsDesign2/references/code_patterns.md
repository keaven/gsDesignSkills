# Code Patterns for gsDesign2

**Note**: These patterns target gsDesign2 >= 1.1.8 (main branch at github.com/Merck/gsDesign2).

## Table of Contents
1. [Enrollment and failure rate setup](#rate-setup)
2. [Fixed design](#fixed-design)
3. [Group sequential design with gs_design_ahr](#gs-design-ahr)
4. [Power computation with gs_power_ahr](#gs-power-ahr)
5. [Non-proportional hazards scenarios](#nph-scenarios)
6. [Spending functions and bounds](#spending-bounds)
7. [Spending time (decoupling spending from information fraction)](#spending-time)
8. [Rate difference designs](#rate-difference)
9. [Weighted logrank and combination tests](#wlr-combo)
10. [Integer sample size rounding](#integer-rounding)
11. [Sequential p-values](#sequential-p-values)
12. [Updating bounds at analysis (gs_update_ahr)](#update-bounds)
13. [Conditional power](#conditional-power)
14. [Output and reporting](#output-reporting)
15. [info_scale options](#info-scale)
16. [Stratified designs](#stratified-designs)

---

## Enrollment and failure rate setup {#rate-setup}

gsDesign2 uses `define_enroll_rate()` and `define_fail_rate()` tibbles
to specify piecewise constant rates.

```r
library(gsDesign2)
library(gsDesign)

# Piecewise enrollment with ramp-up
enroll_rate <- define_enroll_rate(
  duration = c(2, 2, 2, 6),    # Period durations
  rate = c(8, 12, 16, 24)      # Patients per month per period
)

# Exponential failure with proportional hazards
fail_rate <- define_fail_rate(
  duration = Inf,               # Single period (exponential)
  fail_rate = log(2) / 12,     # Control median = 12 months
  hr = 0.7,                    # Hazard ratio under H1
  dropout_rate = 0.001
)

# Delayed treatment effect (non-proportional hazards)
fail_rate_nph <- define_fail_rate(
  duration = c(3, Inf),         # Changepoint at 3 months
  fail_rate = log(2) / c(8, 14),  # Piecewise control hazard
  hr = c(0.9, 0.6),            # HR = 0.9 early, 0.6 after delay
  dropout_rate = 0.001
)
```

### Scaling enrollment to target N

```r
# Set relative rates, then scale to target sample size
enroll_rate <- define_enroll_rate(
  duration = c(2, 2, 2, 18),
  rate = 1:4  # Relative rates
)
target_n <- 400
enroll_rate$rate <- enroll_rate$rate * target_n / sum(enroll_rate$duration * enroll_rate$rate)
```

### Stratified rates

```r
# Two strata with different enrollment and failure rates
enroll_rate_strat <- define_enroll_rate(
  stratum = rep(c("BM+", "BM-"), each = 4),
  duration = rep(c(2, 2, 2, 6), 2),
  rate = c(3, 5, 7, 10, 2, 3, 4, 6)  # Different rates per stratum
)

fail_rate_strat <- define_fail_rate(
  stratum = c("BM+", "BM-"),
  duration = c(Inf, Inf),
  fail_rate = log(2) / 12,
  hr = c(0.6, 0.85),           # Different HR per stratum
  dropout_rate = 0.001
)
```

## Fixed design {#fixed-design}

For single-analysis designs without interim looks.

```r
# AHR-based fixed design
fx <- fixed_design_ahr(
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  alpha = 0.025,
  power = 0.9,
  ratio = 1,
  study_duration = 36
)
fx |> summary() |> as_gt()

# Rate difference (binary endpoint)
fx_rd <- fixed_design_rd(
  alpha = 0.025,
  power = 0.9,
  p_c = 0.15,
  p_e = 0.30,
  rd0 = 0,
  n = NULL  # Solve for sample size
)
```

## Group sequential design with gs_design_ahr {#gs-design-ahr}

`gs_design_ahr()` derives sample size targeting a specified power.

### Basic design with information fraction

```r
x <- gs_design_ahr(
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  alpha = 0.025,
  beta = 0.1,                   # 90% power
  ratio = 1,
  info_frac = c(0.5, 0.75, 1), # Information fractions
  analysis_time = 36,           # Final analysis at 36 months
  info_scale = "h0_info",       # Use null info for spending
  upper = gs_spending_bound,
  upar = list(sf = sfLDOF, total_spend = 0.025),
  lower = gs_spending_bound,
  lpar = list(sf = sfHSD, total_spend = 0.1, param = -2),
  binding = FALSE               # Non-binding futility
)
x |> summary() |> as_gt()
```

### Design driven by analysis times (info_frac = NULL)

When `info_frac = NULL`, the information fraction is derived from
enrollment/failure rate assumptions and analysis timing.

```r
x <- gs_design_ahr(
  enroll_rate = enroll_rate,
  fail_rate = fail_rate_nph,
  alpha = 0.025,
  beta = 0.1,
  info_frac = NULL,                    # Derived from analysis_time
  analysis_time = c(18, 24, 30, 36),   # Calendar months
  info_scale = "h0_info",
  upper = gs_spending_bound,
  upar = list(sf = sfLDOF, total_spend = 0.025),
  lower = gs_spending_bound,
  lpar = list(sf = sfHSD, total_spend = 0.1, param = -2),
  binding = FALSE
)
```

### One-sided design (no futility)

```r
x <- gs_design_ahr(
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  alpha = 0.025,
  beta = 0.1,
  info_frac = c(0.5, 0.75, 1),
  analysis_time = 36,
  info_scale = "h0_info",
  upper = gs_spending_bound,
  upar = list(sf = sfLDOF, total_spend = 0.025),
  test_lower = FALSE              # No futility bound
)
```

### Spending functions as character strings

```r
# Functions can be passed by name (string) instead of function reference
x <- gs_design_ahr(
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  alpha = 0.025, beta = 0.1,
  info_frac = c(0.5, 0.75, 1),
  analysis_time = 36,
  upper = "gs_spending_bound",
  upar = list(sf = "sfLDOF", total_spend = 0.025),
  lower = "gs_spending_bound",
  lpar = list(sf = "sfHSD", total_spend = 0.1, param = -2),
  binding = FALSE
)
```

## Power computation with gs_power_ahr {#gs-power-ahr}

`gs_power_ahr()` computes power for given enrollment/failure rates and timing.
Does not solve for sample size.

```r
pwr <- gs_power_ahr(
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  ratio = 1,
  event = NULL,                   # IMPORTANT: must be NULL to use analysis_time
  analysis_time = c(20, 28, 36),
  info_scale = "h0_info",
  upper = gs_spending_bound,
  upar = list(sf = sfLDOF, total_spend = 0.025),
  lower = gs_spending_bound,
  lpar = list(sf = sfHSD, total_spend = 0.1, param = -2),
  test_lower = TRUE,
  binding = FALSE
)
pwr |> summary() |> as_gt()
```

### Important: event = NULL

`gs_power_ahr()` has a default `event = c(30, 40, 50)` which can cause
length mismatches. **Always set `event = NULL`** when using `analysis_time`
to drive the design.

## Non-proportional hazards scenarios {#nph-scenarios}

### Delayed treatment effect

```r
fail_rate_delayed <- define_fail_rate(
  duration = c(4, Inf),
  fail_rate = log(2) / 10,
  hr = c(1, 0.6),              # No effect first 4 months, then HR = 0.6
  dropout_rate = 0.001
)
```

### Crossing survival curves

```r
fail_rate_crossing <- define_fail_rate(
  duration = c(6, Inf),
  fail_rate = log(2) / c(9, 12),
  hr = c(1.2, 0.5),            # Experimental worse early, better late
  dropout_rate = 0.001
)
```

### Diminishing treatment effect

```r
fail_rate_diminish <- define_fail_rate(
  duration = c(6, Inf),
  fail_rate = log(2) / 10,
  hr = c(0.5, 0.8),            # Strong early effect, weaker late
  dropout_rate = 0.001
)
```

## Spending functions and bounds {#spending-bounds}

### Spending bound specification

```r
# Efficacy: Lan-DeMets O'Brien-Fleming
upper = gs_spending_bound
upar = list(sf = sfLDOF, total_spend = 0.025)

# Futility: Hwang-Shih-DeCani with gamma = -2
lower = gs_spending_bound
lpar = list(sf = sfHSD, total_spend = 0.1, param = -2)

# Fixed bounds (not spending-based)
upper = gs_b
upar = c(3, 2.5, 2)  # Fixed Z-value bounds

# No lower bound
test_lower = FALSE
```

### Common spending function choices

| Spending function | Parameter | Behavior |
|-------------------|-----------|----------|
| `sfLDOF` | None | Lan-DeMets O'Brien-Fleming (conservative early) |
| `sfHSD` | `param = -4` | Hwang-Shih-DeCani, O'Brien-Fleming-like |
| `sfHSD` | `param = -2` | Hwang-Shih-DeCani, moderate |
| `sfHSD` | `param = 1` | Hwang-Shih-DeCani, Pocock-like |
| `sfPoints` | `param = c(...)` | Custom piecewise linear spending |

## Spending time (decoupling spending from information fraction) {#spending-time}

By default, alpha spending tracks the information fraction. For delayed
treatment effects or multiple hypotheses, it may be desirable to use a
different spending schedule. Use `spending_time` in `upar`/`lpar`.

```r
# Spend more conservatively at interims using custom spending time
x <- gs_design_ahr(
  enroll_rate = enroll_rate,
  fail_rate = fail_rate_nph,
  alpha = 0.025,
  beta = 0.1,
  info_frac = NULL,
  analysis_time = c(18, 24, 36),
  info_scale = "h0_info",
  upper = gs_spending_bound,
  upar = list(
    sf = sfLDOF,
    total_spend = 0.025,
    timing = c(0.4, 0.65, 1)    # Spend slower than information fraction
  ),
  lower = gs_spending_bound,
  lpar = list(
    sf = sfHSD,
    total_spend = 0.1,
    param = -2,
    timing = c(0.4, 0.65, 1)
  ),
  binding = FALSE
)
```

## Rate difference designs {#rate-difference}

For binary endpoint designs using risk difference.

```r
# Design for risk difference
x_rd <- gs_design_rd(
  p_c = tibble::tibble(stratum = "All", rate = 0.15),
  p_e = tibble::tibble(stratum = "All", rate = 0.30),
  alpha = 0.025,
  beta = 0.1,
  ratio = 1,
  info_frac = c(0.5, 1),
  rd0 = 0,
  upper = gs_spending_bound,
  upar = list(sf = sfLDOF, total_spend = 0.025),
  lower = gs_b,
  lpar = rep(-Inf, 2)
)

# Power for rate difference
pwr_rd <- gs_power_rd(
  p_c = tibble::tibble(stratum = "All", rate = 0.15),
  p_e = tibble::tibble(stratum = "All", rate = 0.30),
  n = 200,
  ratio = 1,
  info_frac = c(0.5, 1),
  rd0 = 0,
  upper = gs_spending_bound,
  upar = list(sf = sfLDOF, total_spend = 0.025),
  lower = gs_b,
  lpar = rep(-Inf, 2)
)
```

## Weighted logrank and combination tests {#wlr-combo}

### Fleming-Harrington weighted logrank

```r
x_wlr <- gs_design_wlr(
  enroll_rate = enroll_rate,
  fail_rate = fail_rate_nph,
  alpha = 0.025,
  beta = 0.1,
  ratio = 1,
  weight = function(x, arm0, arm1) {
    wlr_weight_fh(x, arm0, arm1, rho = 0, gamma = 0.5)
  },
  info_frac = c(0.5, 0.75, 1),
  analysis_time = 36,
  upper = gs_spending_bound,
  upar = list(sf = sfLDOF, total_spend = 0.025),
  lower = gs_b,
  lpar = rep(-Inf, 3)
)
```

### MaxCombo (combination of FH tests)

```r
x_combo <- gs_design_combo(
  enroll_rate = enroll_rate,
  fail_rate = fail_rate_nph,
  alpha = 0.025,
  beta = 0.1,
  ratio = 1,
  fh_test = rbind(
    data.frame(rho = 0, gamma = 0, tau = -1, test = 1, analysis = 1:3, analysis_time = 36),
    data.frame(rho = 0, gamma = 0.5, tau = -1, test = 2, analysis = 1:3, analysis_time = 36)
  ),
  info_frac = c(0.5, 0.75, 1),
  upper = gs_spending_combo,
  upar = list(sf = sfLDOF, total_spend = 0.025),
  lower = gs_b,
  lpar = rep(-Inf, 3)
)
```

## Integer sample size rounding {#integer-rounding}

```r
# Round to integers (events rounded at IA, rounded up at FA)
x_int <- x |> to_integer()

# With unequal randomization
x_int <- x |> to_integer(ratio = 2)
```

## Sequential p-values {#sequential-p-values}

`sequential_pval()` computes the minimum alpha at which an efficacy bound
would be rejected. Used with graphical multiplicity methods (Maurer-Bretz).

```r
# Design
x <- gs_design_ahr(
  enroll_rate = define_enroll_rate(duration = c(2, 2, 2, 6), rate = c(2.5, 5, 7.5, 10)),
  fail_rate = define_fail_rate(duration = Inf, fail_rate = log(2) / 6,
                               hr = 0.6, dropout_rate = .001),
  info_frac = c(.5, .65, .8, 1),
  analysis_time = 30,
  upper = gs_spending_bound,
  upar = list(sf = "sfLDOF", total_spend = 0.025),
  lower = gs_spending_bound,
  lpar = list(sf = "sfHSD", total_spend = 0.1, param = 2),
  binding = FALSE
)

# Compute sequential p-value from observed data
sequential_pval(
  gs_design = x,
  event = c(100, 160, 190),       # Observed events at each analysis
  z = c(1.5, 2, 2.5)              # Observed Z-statistics (positive = favorable)
)

# At interim (partial data)
sequential_pval(gs_design = x, event = c(100, 160), z = c(1.5, 2))

# With custom spending time
sequential_pval(
  gs_design = x,
  event = c(100, 160, 190, 230),
  z = c(1.5, 2, 2.5, 3),
  ustime = c(0.4, 0.6, 0.8, 1)
)
```

## Updating bounds at analysis (gs_update_ahr) {#update-bounds}

When observed events differ from planned, update bounds with `gs_update_ahr()`.

```r
library(gsDesign)
library(gsDesign2)

# Original design
x <- gs_design_ahr(
  enroll_rate = define_enroll_rate(duration = c(2, 2, 10), rate = (1:3) / 3),
  fail_rate = define_fail_rate(duration = c(3, Inf), fail_rate = log(2) / 9,
                               hr = c(1, 0.6), dropout_rate = .0001),
  alpha = 0.025, beta = 0.1,
  info_frac = NULL,
  analysis_time = c(20, 36),
  info_scale = "h0_info",
  upper = gs_spending_bound,
  upar = list(sf = sfLDOF, total_spend = 0.025),
  test_lower = FALSE
)

# Update with observed events and spending time
x_upd <- gs_update_ahr(
  x = x,
  alpha = 0.025,
  ustime = c(0.8, 1),               # Spending time at observed analyses
  event_tbl = data.frame(
    analysis = c(1, 2),
    event = c(180, 250)              # Observed events
  )
)
```

### Stratified update with piecewise events

```r
# Piecewise event table (events per interval per analysis)
event_tbl <- data.frame(
  analysis = c(1, 1, 2, 2),
  event = c(30, 100, 30, 200)        # 30 events in delayed period, rest after
)

x_upd <- gs_update_ahr(
  x = x,
  alpha = 0.025,
  ustime = c(0.55, 1),
  event_tbl = event_tbl
)
```

## Conditional power {#conditional-power}

```r
# Conditional power under NPH
cp <- gs_cp_npe(
  theta = x$analysis$theta,
  info = x$analysis$info,
  info0 = x$analysis$info0,
  z = 2.0,                    # Observed Z at interim
  i = 1,                      # Which interim analysis
  x = x
)
```

## Output and reporting {#output-reporting}

```r
# Summary table
x |> summary() |> as_gt()

# Bound summary (similar to gsDesign::gsBoundSummary)
x |> gs_bound_summary()

# Text summary
x |> text_summary()

# RTF output
x |> summary() |> as_rtf(file = "design.rtf")

# Suppress footnotes
x |> summary() |> as_gt(footnote = FALSE)

# Multiple alpha levels in bound summary
x |> gs_bound_summary(alpha = c(0.025, 0.0125))
```

## info_scale options {#info-scale}

Controls how information is computed for bound calculations.

| info_scale | Description | When to use |
|------------|-------------|-------------|
| `"h0_h1_info"` | Average of H0 and H1 info (default) | General use |
| `"h0_info"` | Null hypothesis info only | Matches gsDesign convention; recommended for multiplicity |
| `"h1_info"` | Alternative hypothesis info only | Rarely used |

```r
# Match gsDesign convention
x <- gs_design_ahr(
  ...,
  info_scale = "h0_info"
)
```

## Stratified designs {#stratified-designs}

For overall population designs with subgroups, use stratified enrollment
and failure rates.

```r
# Subgroup prevalence
prevalence <- 0.4

# Subgroup enrollment from a driving design
enroll_rate_sub <- x_sub$enroll_rate  # From gs_design_ahr output

# Build overall enrollment (BM+ and BM- strata)
enroll_rate_overall <- define_enroll_rate(
  stratum = rep(c("BM+", "BM-"), each = nrow(enroll_rate_sub)),
  duration = rep(enroll_rate_sub$duration, 2),
  rate = c(
    enroll_rate_sub$rate,
    enroll_rate_sub$rate * (1 - prevalence) / prevalence
  )
)

# Stratified failure rates
fail_rate_overall <- define_fail_rate(
  stratum = c("BM+", "BM-"),
  duration = c(Inf, Inf),
  fail_rate = log(2) / 12,
  hr = c(0.6, 0.85),
  dropout_rate = 0.001
)

# Power for overall population
pwr_overall <- gs_power_ahr(
  enroll_rate = enroll_rate_overall,
  fail_rate = fail_rate_overall,
  ratio = 1,
  event = NULL,
  analysis_time = c(20, 28, 36),
  info_scale = "h0_info",
  upper = gs_spending_bound,
  upar = list(sf = sfLDOF, total_spend = 0.025),
  test_lower = FALSE,
  binding = FALSE
)
```

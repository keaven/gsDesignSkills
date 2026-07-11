# Code Patterns for simtrial

**Note**: These patterns target simtrial >= 1.0.2 (main branch at github.com/Merck/simtrial).

## Table of Contents
1. [Generating survival data with sim_pw_surv](#sim-pw-surv)
2. [Data cutting](#data-cutting)
3. [Weighted logrank test (wlr)](#wlr)
4. [MaxCombo test](#maxcombo)
5. [Milestone and RMST tests](#milestone-rmst)
6. [Multiple tests with multitest](#multitest)
7. [Fixed-sample simulation with sim_fixed_n](#sim-fixed-n)
8. [Group sequential simulation with sim_gs_n](#sim-gs-n)
9. [Integrating with gsDesign2 designs](#gsdesign2-integration)
10. [Weight functions](#weight-functions)
11. [Converting rate formats with to_sim_pw_surv](#to-sim-pw-surv)
12. [Flexible analysis timing with get_analysis_date](#analysis-date)
13. [Example NPH datasets](#example-datasets)
14. [Stratified simulations](#stratified)

---

## Generating survival data with sim_pw_surv {#sim-pw-surv}

`sim_pw_surv()` generates individual patient data with piecewise exponential
enrollment, failure, and dropout.

### Basic usage (proportional hazards)

```r
library(simtrial)

# Simulate 400 patients, 1:1 randomization
data <- sim_pw_surv(
  n = 400,
  stratum = data.frame(stratum = "All", p = 1),
  block = c(rep("control", 2), rep("experimental", 2)),
  enroll_rate = data.frame(rate = 25, duration = 16),
  fail_rate = data.frame(
    stratum = rep("All", 2),
    period = c(1, 1),
    treatment = c("control", "experimental"),
    duration = c(100, 100),
    rate = c(log(2) / 12, log(2) / 12 * 0.7)
  ),
  dropout_rate = data.frame(
    stratum = rep("All", 2),
    period = c(1, 1),
    treatment = c("control", "experimental"),
    duration = c(100, 100),
    rate = c(0.001, 0.001)
  )
)
# Returns: stratum, enroll_time, treatment, fail_time, dropout_time, cte, fail
```

### Delayed treatment effect (non-proportional hazards)

```r
data <- sim_pw_surv(
  n = 400,
  stratum = data.frame(stratum = "All", p = 1),
  block = c(rep("control", 2), rep("experimental", 2)),
  enroll_rate = data.frame(rate = c(10, 30), duration = c(4, 12)),
  fail_rate = data.frame(
    stratum = rep("All", 4),
    period = rep(1:2, 2),
    treatment = c(rep("control", 2), rep("experimental", 2)),
    duration = rep(c(3, 100), 2),
    rate = log(2) / c(9, 9, 9, 18)   # HR = 1 first 3 months, then 0.5
  ),
  dropout_rate = data.frame(
    stratum = rep("All", 2),
    period = c(1, 1),
    treatment = c("control", "experimental"),
    duration = c(100, 100),
    rate = c(0.001, 0.001)
  )
)
```

### Using to_sim_pw_surv for simpler rate specification

```r
# Define rates in the simpler format (like sim_fixed_n)
fail_rate <- data.frame(
  stratum = "All",
  duration = c(3, 100),
  fail_rate = log(2) / c(9, 18),
  hr = c(0.9, 0.6),
  dropout_rate = 0.001
)

# Convert to sim_pw_surv format
rates <- to_sim_pw_surv(fail_rate)

data <- sim_pw_surv(
  n = 400,
  enroll_rate = data.frame(rate = c(10, 30), duration = c(4, 12)),
  fail_rate = rates$fail_rate,
  dropout_rate = rates$dropout_rate
)
```

---

## Data cutting {#data-cutting}

After generating data with `sim_pw_surv()`, cut it at a specific calendar time
or event count to create an analysis dataset.

### Cut by calendar date

```r
# Cut at 36 months from study start
analysis_data <- cut_data_by_date(data, cut_date = 36)
# Returns: tte, event, stratum, treatment (class "tte_data")
```

### Cut by event count

```r
# Cut when 250 events observed
analysis_data <- cut_data_by_event(data, event = 250)
```

### Using create_cut for pipelines

`create_cut()` wraps `get_analysis_date()` arguments into a reusable function.

```r
# Cut at planned calendar time
cut1 <- create_cut(planned_calendar_time = 24)

# Cut at target events
cut2 <- create_cut(target_event_overall = 250)

# Cut at target events but no more than 6 months extension
cut3 <- create_cut(
  target_event_overall = 300,
  max_extension_for_target_event = 6
)

# Apply cutting function to simulated data
analysis_data <- cut1(data)
```

---

## Weighted logrank test (wlr) {#wlr}

`wlr()` performs a single weighted logrank test on cut data.

```r
# Standard logrank
result <- wlr(
  data = cut_data_by_event(data, 250),
  weight = fh(rho = 0, gamma = 0)
)
result$z  # Z-statistic (negative favors experimental)

# Fleming-Harrington (0, 0.5) - emphasizes late differences
result_fh <- wlr(
  data = cut_data_by_event(data, 250),
  weight = fh(rho = 0, gamma = 0.5)
)

# With variance estimate
result_var <- wlr(
  data = cut_data_by_event(data, 250),
  weight = fh(rho = 0, gamma = 0),
  return_variance = TRUE
)
```

---

## MaxCombo test {#maxcombo}

`maxcombo()` combines multiple Fleming-Harrington weighted logrank tests,
using the maximum test statistic with a correlation-adjusted p-value.

```r
# MaxCombo: logrank + FH(0, 0.5)
result <- maxcombo(
  data = cut_data_by_event(data, 250),
  rho = c(0, 0),
  gamma = c(0, 0.5)
)
result$z        # Z-statistics for each test
result$p_value  # p-value (accounts for correlation)

# Three-test MaxCombo
result3 <- maxcombo(
  data = cut_data_by_event(data, 250),
  rho = c(0, 0, 1),
  gamma = c(0, 1, 1),
  return_corr = TRUE   # Also return correlation matrix
)
```

---

## Milestone and RMST tests {#milestone-rmst}

### Milestone (survival difference at fixed time)

```r
result <- milestone(
  data = cut_data_by_event(data, 250),
  ms_time = 12,            # Compare survival at 12 months
  test_type = "naive"      # Or "log-log"
)
result$estimate  # Survival difference (experimental - control)
result$z         # Z-statistic
```

### RMST (restricted mean survival time difference)

```r
result <- rmst(
  data = cut_data_by_event(data, 250),
  tau = 24                 # RMST up to 24 months
)
```

---

## Multiple tests with multitest {#multitest}

Apply several tests to the same cut dataset in one call.

```r
# Create parameterized test functions
logrank <- create_test(wlr, weight = fh(rho = 0, gamma = 0))
fh05 <- create_test(wlr, weight = fh(rho = 0, gamma = 0.5))
rmst_test <- create_test(rmst, tau = 24)
combo <- create_test(maxcombo, rho = c(0, 0), gamma = c(0, 0.5))

# Apply all tests to one dataset
cut <- cut_data_by_event(data, 250)
results <- multitest(
  data = cut,
  logrank = logrank,
  fh05 = fh05,
  rmst = rmst_test,
  combo = combo
)
# results is a named list of test outputs
```

---

## Fixed-sample simulation with sim_fixed_n {#sim-fixed-n}

`sim_fixed_n()` is a self-contained pipeline: simulate data, cut, and test
in one call. Simpler but less flexible than `sim_gs_n()`.

```r
# Basic simulation with Fleming-Harrington tests
results <- sim_fixed_n(
  n_sim = 1000,
  sample_size = 400,
  target_event = 300,
  enroll_rate = data.frame(duration = c(4, 12), rate = c(10, 30)),
  fail_rate = data.frame(
    stratum = "All",
    duration = c(3, 100),
    fail_rate = log(2) / c(9, 18),
    hr = c(0.9, 0.6),
    dropout_rate = 0.001
  ),
  total_duration = 36,
  timing_type = 1:5,
  rho_gamma = data.frame(rho = c(0, 0), gamma = c(0, 0.5))
)

# timing_type controls the cutoff method:
# 1: Planned study duration (total_duration)
# 2: Time target_event is achieved
# 3: Minimum follow-up after enrollment complete
# 4: max(planned duration, target event time)
# 5: max(target event time, minimum follow-up)
```

### Estimating power from sim_fixed_n

```r
# Power = proportion of simulations with z < -1.96
power <- results |>
  dplyr::filter(cut == "Planned duration" & rho == 0 & gamma == 0) |>
  dplyr::summarize(power = mean(z <= qnorm(0.025)))
```

---

## Group sequential simulation with sim_gs_n {#sim-gs-n}

`sim_gs_n()` simulates a group sequential trial: generate data, cut at
multiple analysis times, and apply tests.

### Basic group sequential simulation

```r
library(simtrial)
library(gsDesign2)
library(gsDesign)

# Define cutting functions for 3 analyses
ia1_cut <- create_cut(target_event_overall = 150)
ia2_cut <- create_cut(target_event_overall = 225)
fa_cut <- create_cut(target_event_overall = 300)

# Simulate with logrank test
results <- sim_gs_n(
  n_sim = 1000,
  sample_size = 400,
  enroll_rate = data.frame(duration = c(4, 12), rate = c(10, 30)),
  fail_rate = data.frame(
    stratum = "All",
    duration = c(3, 100),
    fail_rate = log(2) / c(9, 18),
    hr = c(0.9, 0.6),
    dropout_rate = 0.001
  ),
  test = wlr,
  cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
  weight = fh(rho = 0, gamma = 0)
)
```

### MaxCombo in sim_gs_n

MaxCombo must be the only test (cannot be combined with other tests in the
same call).

```r
results_combo <- sim_gs_n(
  n_sim = 1000,
  sample_size = 400,
  enroll_rate = data.frame(duration = c(4, 12), rate = c(10, 30)),
  fail_rate = data.frame(
    stratum = "All",
    duration = c(3, 100),
    fail_rate = log(2) / c(9, 18),
    hr = c(0.9, 0.6),
    dropout_rate = 0.001
  ),
  test = maxcombo,
  cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut),
  rho = c(0, 0),
  gamma = c(0, 0.5)
)
```

### Multiple WLR tests in sim_gs_n

```r
results_multi <- sim_gs_n(
  n_sim = 1000,
  sample_size = 400,
  enroll_rate = data.frame(duration = c(4, 12), rate = c(10, 30)),
  fail_rate = data.frame(
    stratum = "All",
    duration = c(3, 100),
    fail_rate = log(2) / c(9, 18),
    hr = c(0.9, 0.6),
    dropout_rate = 0.001
  ),
  test = list(
    logrank = create_test(wlr, weight = fh(rho = 0, gamma = 0)),
    fh05 = create_test(wlr, weight = fh(rho = 0, gamma = 0.5))
  ),
  cut = list(ia1 = ia1_cut, ia2 = ia2_cut, fa = fa_cut)
)
```

---

## Integrating with gsDesign2 designs {#gsdesign2-integration}

Use a `gsDesign2` design to drive sample size, event targets, and boundary
updates in `sim_gs_n()`.

```r
library(simtrial)
library(gsDesign2)
library(gsDesign)

# Enrollment and failure rate assumptions
enroll_rate <- define_enroll_rate(
  duration = c(4, 12),
  rate = c(10, 30)
)

fail_rate <- define_fail_rate(
  duration = c(3, 100),
  fail_rate = log(2) / 9,
  hr = c(0.9, 0.6),
  dropout_rate = 0.001
)

# Group sequential design
design <- gs_design_ahr(
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  alpha = 0.025,
  beta = 0.1,
  analysis_time = c(12, 24, 36),
  info_scale = "h0_info",
  upper = gs_spending_bound,
  upar = list(sf = sfLDOF, total_spend = 0.025),
  lower = gs_b,
  lpar = rep(-Inf, 3)
)

# Create cuts from the design's planned events
cuts <- list(
  ia1 = create_cut(target_event_overall = ceiling(design$analysis$event[1])),
  ia2 = create_cut(target_event_overall = ceiling(design$analysis$event[2])),
  fa = create_cut(target_event_overall = ceiling(design$analysis$event[3]))
)

# Simulate with boundary updating from original design
results <- sim_gs_n(
  n_sim = 1000,
  sample_size = ceiling(max(design$analysis$n)),
  enroll_rate = data.frame(
    duration = design$enroll_rate$duration,
    rate = design$enroll_rate$rate
  ),
  fail_rate = data.frame(
    stratum = "All",
    duration = fail_rate$duration,
    fail_rate = fail_rate$fail_rate,
    hr = fail_rate$hr,
    dropout_rate = fail_rate$dropout_rate
  ),
  test = wlr,
  cut = cuts,
  original_design = design,
  weight = fh(rho = 0, gamma = 0)
)

# Summarize results
results |> summary() |> as_gt()
```

### Alpha spending options

```r
# How interim alpha is spent when observed events differ from planned
results <- sim_gs_n(
  ...,
  original_design = design,
  ia_alpha_spending = "min_planned_actual",  # Default: min of planned and actual info fraction
  fa_alpha_spending = "full_alpha"           # Default: spend full alpha at final analysis
)

# Alternative: spend by actual info fraction at FA (for underrunning)
results <- sim_gs_n(
  ...,
  original_design = design,
  ia_alpha_spending = "actual",
  fa_alpha_spending = "info_frac"
)
```

---

## Weight functions {#weight-functions}

### Fleming-Harrington family

```r
fh(rho = 0, gamma = 0)    # Standard logrank
fh(rho = 0, gamma = 0.5)  # Emphasizes late differences
fh(rho = 0, gamma = 1)    # Stronger late emphasis (Peto-logrank-like)
fh(rho = 1, gamma = 0)    # Emphasizes early differences (Prentice-Wilcoxon)
fh(rho = 1, gamma = 1)    # Moderate weight all around
```

### Magirr-Burman weights (for delayed effects)

```r
mb(delay = 4)              # Down-weight first 4 months
mb(delay = 4, w_max = 2)  # Cap maximum weight at 2
mb(delay = Inf, w_max = 2) # Magirr (2021) recommendation
```

### Early zero weight

```r
early_zero(early_period = 4)  # Zero weight before 4 months
```

---

## Converting rate formats with to_sim_pw_surv {#to-sim-pw-surv}

`sim_fixed_n()` and `sim_gs_n()` use a simpler rate format (control rate + HR).
`sim_pw_surv()` uses treatment-specific rates. `to_sim_pw_surv()` converts
between them.

```r
# Simple format (control rate + HR)
fail_rate_simple <- data.frame(
  stratum = "All",
  duration = c(3, 100),
  fail_rate = log(2) / c(9, 18),
  hr = c(0.9, 0.6),
  dropout_rate = 0.001
)

# Convert to sim_pw_surv format
rates <- to_sim_pw_surv(fail_rate_simple)
rates$fail_rate    # Treatment-specific failure rates
rates$dropout_rate # Treatment-specific dropout rates

# Use in sim_pw_surv
data <- sim_pw_surv(
  n = 400,
  enroll_rate = data.frame(rate = 25, duration = 16),
  fail_rate = rates$fail_rate,
  dropout_rate = rates$dropout_rate
)
```

---

## Flexible analysis timing with get_analysis_date {#analysis-date}

`get_analysis_date()` determines analysis timing from multiple criteria.
`create_cut()` wraps it for use in pipelines.

### Calendar-time-driven analysis

```r
cut <- create_cut(planned_calendar_time = 36)
```

### Event-driven analysis

```r
cut <- create_cut(target_event_overall = 300)
```

### Event-driven with maximum extension

```r
# Target 300 events, but extend at most 6 months beyond planned time
cut <- create_cut(
  target_event_overall = 300,
  max_extension_for_target_event = 42  # Absolute max calendar time
)
```

### Minimum follow-up

```r
# At least 12 months follow-up after last enrollment
cut <- create_cut(min_followup = 12)
```

### Combining criteria

```r
# Analysis at the later of target events or minimum follow-up
cut <- create_cut(
  target_event_overall = 300,
  min_followup = 12
)

# Minimum time between analyses
cut_ia2 <- create_cut(
  target_event_overall = 225,
  min_time_after_previous_analysis = 6
)
```

### Per-stratum event targets

```r
cut <- create_cut(
  target_event_per_stratum = c("BM+" = 100, "BM-" = 150)
)
```

---

## Example NPH datasets {#example-datasets}

simtrial includes pre-built datasets from the Cross-Pharma Non-Proportional
Hazards Working Group.

```r
# Delayed treatment effect
data(ex1_delayed_effect)
# Columns: id, month (tte), evntd (event indicator), trt (0=ctrl, 1=exp)

data(ex2_delayed_effect)  # Another delayed effect scenario
data(ex3_cure_with_ph)    # Cure model with proportional hazards
data(ex4_belly)           # "Belly-shaped" survival curves
data(ex5_widening)        # Widening survival curves
data(ex6_crossing)        # Crossing survival curves

# Use with wlr (need column renaming)
library(dplyr)
analysis_data <- ex1_delayed_effect |>
  rename(tte = month, event = evntd, treatment = trt) |>
  mutate(
    treatment = ifelse(treatment == 1, "experimental", "control"),
    stratum = "All"
  )
wlr(analysis_data, weight = fh(rho = 0, gamma = 0))
```

---

## Stratified simulations {#stratified}

### Generating stratified data

```r
data <- sim_pw_surv(
  n = 600,
  stratum = data.frame(stratum = c("BM+", "BM-"), p = c(0.4, 0.6)),
  block = c(rep("control", 2), rep("experimental", 2)),
  enroll_rate = data.frame(rate = 30, duration = 20),
  fail_rate = data.frame(
    stratum = rep(c("BM+", "BM+", "BM-", "BM-"), each = 1),
    period = rep(1, 4),
    treatment = rep(c("control", "experimental"), 2),
    duration = rep(100, 4),
    rate = c(log(2)/12, log(2)/12 * 0.6,   # BM+: HR = 0.6
             log(2)/12, log(2)/12 * 0.85)   # BM-: HR = 0.85
  ),
  dropout_rate = data.frame(
    stratum = rep(c("BM+", "BM-"), each = 2),
    period = rep(1, 4),
    treatment = rep(c("control", "experimental"), 2),
    duration = rep(100, 4),
    rate = rep(0.001, 4)
  )
)

# Cut and test with stratified WLR
cut <- cut_data_by_event(data, event = 350)
wlr(cut, weight = fh(rho = 0, gamma = 0))
```

### Stratified sim_fixed_n

```r
results <- sim_fixed_n(
  n_sim = 1000,
  sample_size = 600,
  target_event = 350,
  stratum = data.frame(stratum = c("BM+", "BM-"), p = c(0.4, 0.6)),
  enroll_rate = data.frame(duration = c(4, 16), rate = c(10, 30)),
  fail_rate = data.frame(
    stratum = rep(c("BM+", "BM-"), each = 2),
    duration = rep(c(3, 100), 2),
    fail_rate = log(2) / 12,
    hr = c(0.6, 0.6, 0.85, 0.85),
    dropout_rate = 0.001
  ),
  total_duration = 36,
  timing_type = 1:2,
  rho_gamma = data.frame(rho = 0, gamma = 0)
)
```

---

## Standalone wlr with illness-death model data {#wlr-illness-death}

`wlr()` can be used standalone (outside `sim_gs_n()`) with data from the
illness-death model. The key is mapping ADTTE columns to the `wlr()` format.

### Column mapping

`wlr()` expects columns: `tte`, `event`, `stratum`, `treatment`

```r
# From illness-death ADTTE data (columns: AVAL, CNSR, STRATUM, TRT)
logrank_z <- function(data) {
  d <- data.frame(
    tte = data$AVAL,
    event = 1L - data$CNSR,
    stratum = data$STRATUM,
    treatment = data$TRT
  )
  if (sum(d$event) < 2) return(NA_real_)
  wlr(d, weight = fh(rho = 0, gamma = 0))$z
}
```

### Sign convention

`wlr()` returns **positive Z when experimental is better** (more events in
control than expected under null). This matches the convention needed for
`sequentialPValue()` and graphical testing.

### Usage in simulation loops

```r
adtte_18 <- cut_illness_death(sim_data, cut_date = 18)
z_pfs_18 <- logrank_z(adtte_18[adtte_18$PARAMCD == "PFS", ])
z_os_18 <- logrank_z(adtte_18[adtte_18$PARAMCD == "OS", ])

# Collect event counts for spending time
ev_pfs_18 <- sum(1L - adtte_18$CNSR[adtte_18$PARAMCD == "PFS"])
ev_os_18 <- sum(1L - adtte_18$CNSR[adtte_18$PARAMCD == "OS"])
```

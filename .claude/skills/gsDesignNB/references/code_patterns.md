# Code Patterns for gsDesignNB

## Table of Contents
1. [Fixed sample size calculation](#sample-size)
2. [Power calculation](#power-calc)
3. [Event gaps and piecewise accrual](#event-gaps)
4. [Non-inferiority and super-superiority](#non-inferiority)
5. [Group sequential design](#group-sequential)
6. [Simulation of recurrent events](#simulation)
7. [Group sequential simulation](#gs-simulation)
8. [Data cutting and analysis timing](#data-cutting)
9. [Wald test for treatment rate ratio](#wald-test)
10. [Blinded and unblinded information](#information)
11. [Sample size re-estimation](#ssr)
12. [Completers-based interim analysis](#completers)
13. [Seasonal event simulation](#seasonal)
14. [Verification of operating characteristics](#verification)

---

## Fixed sample size calculation {#sample-size}

`sample_size_nbinom()` implements the Zhu-Lakkis Method 3 for negative binomial
outcomes, accounting for variable accrual, dropout, maximum follow-up, and event gaps.

### Basic calculation

```r
library(gsDesignNB)

design <- sample_size_nbinom(
  lambda1 = 0.5,          # Control event rate (events/time unit)
  lambda2 = 0.3,          # Experimental event rate
  dispersion = 0.5,       # Overdispersion parameter k
  power = 0.9,            # Target power
  alpha = 0.025,          # One-sided significance level
  accrual_rate = 10,      # Enrollment rate (patients/time unit)
  accrual_duration = 12,  # Total enrollment duration
  trial_duration = 24     # Total trial duration
)
design
summary(design)
```

### Piecewise accrual ramp-up

```r
design <- sample_size_nbinom(
  lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.5,
  power = 0.9, alpha = 0.025,
  accrual_rate = c(5, 10),        # Rate per period
  accrual_duration = c(3, 3),     # Duration per period (total: 6 months)
  trial_duration = 12
)
```

### With dropout and maximum follow-up

```r
design <- sample_size_nbinom(
  lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.5,
  power = 0.9, alpha = 0.025,
  accrual_rate = 10, accrual_duration = 12,
  trial_duration = 24,
  dropout_rate = 0.05,       # Annual dropout rate (scalar: same for both groups)
  max_followup = 12          # Cap follow-up at 12 months
)
```

### Group-specific parameters

```r
design <- sample_size_nbinom(
  lambda1 = 0.5, lambda2 = 0.3,
  dispersion = c(0.5, 0.3),       # Different dispersion per group
  power = 0.9, alpha = 0.025,
  accrual_rate = 10, accrual_duration = 12, trial_duration = 24,
  dropout_rate = c(0.10, 0.05),   # Different dropout per group
  max_followup = c(12, 12)        # Can differ per group
)
```

### Unequal allocation

```r
design <- sample_size_nbinom(
  lambda1 = 0.5, lambda2 = 0.3, dispersion = 0.5,
  power = 0.9, alpha = 0.025, ratio = 2,   # 2:1 experimental:control
  accrual_rate = 10, accrual_duration = 12, trial_duration = 24
)
# design$n1 = control, design$n2 = experimental (n2 = ratio * n1)
```

### Key output fields

```r
design$n_total        # Total sample size
design$n1             # Control group N
design$n2             # Experimental group N
design$power          # Achieved power
design$total_events   # Expected total events
design$events_n1      # Expected control events
design$events_n2      # Expected experimental events
design$exposure       # Average exposure per group
design$variance       # Variance of log rate ratio
```

---

## Power calculation {#power-calc}

Set `power = NULL` to compute power for a given sample size:

```r
power_result <- sample_size_nbinom(
  lambda1 = 0.5, lambda2 = 0.4,
  dispersion = 0.5,
  power = NULL,                     # Request power calculation
  alpha = 0.025,
  accrual_rate = design$accrual_rate,
  accrual_duration = design$inputs$accrual_duration,
  trial_duration = design$inputs$trial_duration,
  dropout_rate = design$inputs$dropout_rate,
  max_followup = design$inputs$max_followup
)
power_result$power
```

---

## Event gaps and piecewise accrual {#event-gaps}

Event gaps model time "off risk" after each event (e.g., recovery period).

```r
design <- sample_size_nbinom(
  lambda1 = 2.0, lambda2 = 1.0,    # Higher rates (events/year)
  dispersion = 0.5,
  power = 0.9, alpha = 0.025,
  accrual_rate = 10, accrual_duration = 12, trial_duration = 24,
  event_gap = 30 / 365.25          # 30-day gap after each event (in years)
)

# Event gaps reduce effective exposure:
# lambda_eff = lambda / (1 + lambda * gap)
# Exposure_at_risk < Exposure_total
design$exposure_at_risk_n1  # Exposure at risk (excludes gap time)
```

---

## Non-inferiority and super-superiority {#non-inferiority}

The `rr0` parameter specifies the null hypothesis rate ratio.

### Non-inferiority (rr0 > 1)

```r
# H0: RR >= 1.1, H1: RR < 1.1
ni_design <- sample_size_nbinom(
  lambda1 = 0.1, lambda2 = 0.09,   # Experimental slightly better
  rr0 = 1.1,                        # Non-inferiority margin
  dispersion = 0.5,
  power = 0.9, alpha = 0.025,
  accrual_rate = 100, accrual_duration = 1, trial_duration = 2
)
```

### Super-superiority (rr0 < 1)

```r
# H0: RR >= 0.5, H1: RR < 0.5 (must show > 50% reduction)
ss_design <- sample_size_nbinom(
  lambda1 = 0.1, lambda2 = 0.02,
  rr0 = 0.5,
  dispersion = 0.5,
  power = 0.9, alpha = 0.025,
  accrual_rate = 100, accrual_duration = 1, trial_duration = 2,
  max_followup = 1
)
```

---

## Group sequential design {#group-sequential}

`gsNBCalendar()` creates a group sequential design for negative binomial outcomes
with calendar-time analysis schedule, compatible with gsDesign spending functions.

### Basic group sequential design

```r
# Step 1: fixed design
nb_ss <- sample_size_nbinom(
  lambda1 = 1.5/12, lambda2 = 1.0/12, dispersion = 0.5,
  power = 0.9, alpha = 0.025,
  accrual_rate = 1, accrual_duration = 12, trial_duration = 24,
  max_followup = 12, event_gap = 20/30.4375
)

# Step 2: group sequential design
gs_nb <- gsNBCalendar(
  x = nb_ss,
  k = 3,                              # 3 analyses
  test.type = 4,                       # 2-sided, non-binding futility
  alpha = 0.025,
  sfu = gsDesign::sfLDOF,             # Upper spending: Lan-DeMets O'Brien-Fleming
  sfl = gsDesign::sfHSD, sflpar = -8, # Lower spending: Hwang-Shih-DeCani
  analysis_times = c(10, 18, 24)       # Calendar months
)

# Step 3: integer sample sizes
gs_nb <- toInteger(gs_nb)

# Display bounds
gsDesign::gsBoundSummary(gs_nb, deltaname = "RR", logdelta = TRUE)
```

### With custom spending times

```r
gs_nb <- gsNBCalendar(
  x = nb_ss,
  k = 3, test.type = 1,               # One-sided (efficacy only)
  sfu = gsDesign::sfLinear,
  sfupar = c(0.5, 0.5),               # Linear spending parameters
  usTime = c(0.1, 0.18, 1),           # Spending time (may differ from info fraction)
  analysis_times = c(10, 18, 24)
)
```

### Key output fields

```r
gs_nb$n.I           # Information at each analysis
gs_nb$n_total       # Total N at each analysis (vector)
gs_nb$n1            # Control N at each analysis
gs_nb$n2            # Experimental N at each analysis
gs_nb$events        # Expected total events at each analysis
gs_nb$exposure      # Expected exposure at each analysis
gs_nb$T             # Calendar times of analyses
gs_nb$upper$bound   # Upper (efficacy) Z bounds
gs_nb$lower$bound   # Lower (futility) Z bounds
```

---

## Simulation of recurrent events {#simulation}

`nb_sim()` simulates a single trial using a Gamma-Poisson mixture model.

### Basic simulation

```r
sim_data <- nb_sim(
  enroll_rate = data.frame(rate = 20, duration = 12),
  fail_rate = data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.5, 0.3),
    dispersion = c(0.5, 0.5)
  ),
  dropout_rate = data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.05, 0.05),
    duration = c(100, 100)
  ),
  max_followup = 12,
  n = 200
)

# Output: one row per event/censoring per patient
head(sim_data)
# id, treatment, enroll_time, tte, calendar_time, event
```

### With enrollment ramp-up and event gaps

```r
sim_data <- nb_sim(
  enroll_rate = data.frame(
    rate = c(5, 10, 20),
    duration = c(2, 2, 8)
  ),
  fail_rate = data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(1.5/12, 1.0/12),
    dispersion = c(0.5, 0.5)
  ),
  dropout_rate = data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.001, 0.001),
    duration = c(100, 100)
  ),
  max_followup = 12,
  event_gap = 5/365.25     # 5-day gap after each event
)
```

---

## Group sequential simulation {#gs-simulation}

`sim_gs_nbinom()` runs many group sequential trials and returns per-analysis results.

### Basic usage

```r
sim_results <- sim_gs_nbinom(
  n_sims = 1000,
  enroll_rate = data.frame(rate = 20, duration = 12),
  fail_rate = data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(1.5/12, 1.0/12),
    dispersion = c(0.5, 0.5)
  ),
  dropout_rate = data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.001, 0.001),
    duration = c(100, 100)
  ),
  max_followup = 12,
  event_gap = 20/30.4375,
  analysis_times = c(10, 18, 24),
  n_target = gs_nb$n_total[3],
  design = gs_nb
)

# Check bounds
sim_checked <- check_gs_bound(sim_results, design = gs_nb, info_scale = "blinded")

# Summarize operating characteristics
oc <- summarize_gs_sim(sim_checked)
oc$power                 # Overall power
oc$analysis_summary      # Per-analysis cumulative crossing probabilities
```

### With custom cutting criteria

```r
sim_results <- sim_gs_nbinom(
  n_sims = 500,
  enroll_rate = data.frame(rate = 20, duration = 12),
  fail_rate = fail_rate_df,
  dropout_rate = dropout_df,
  max_followup = 12,
  event_gap = 5/365.25,
  n_target = 240,
  design = gs_nb,
  cuts = list(
    # Analysis 1: calendar time 10 OR 80 events, whichever later
    list(planned_calendar = 10, target_events = 80),
    # Analysis 2: calendar time 18 OR 160 events
    list(planned_calendar = 18, target_events = 160),
    # Analysis 3: calendar time 24 (final)
    list(planned_calendar = 24)
  )
)
```

---

## Data cutting and analysis timing {#data-cutting}

### Cut at a calendar date

```r
# Aggregate to one row per subject at cut date
cut_summary <- cut_data_by_date(sim_data, cut_date = 18, event_gap = 5/365.25)
# Output: id, treatment, enroll_time, tte, tte_total, events
```

### Find date when target events reached

```r
event_date <- get_analysis_date(sim_data, planned_events = 150, event_gap = 5/365.25)
```

### Find date satisfying multiple criteria

```r
cut_date <- get_cut_date(
  data = sim_data,
  planned_calendar = 18,       # At least 18 months
  target_events = 150,          # At least 150 events
  target_completers = 100,      # At least 100 completers
  event_gap = 5/365.25
)
```

### Find date for completers

```r
comp_date <- cut_date_for_completers(sim_data, target_completers = 80)
cut_comp <- cut_completers(sim_data, cut_date = comp_date, event_gap = 5/365.25)
```

---

## Wald test for treatment rate ratio {#wald-test}

`mutze_test()` performs a Wald test based on negative binomial (or Poisson) GLM.

```r
# Cut data first
cut_summary <- cut_data_by_date(sim_data, cut_date = 24, event_gap = 5/365.25)

# Negative binomial Wald test
result <- mutze_test(cut_summary)
result$z              # Z statistic
result$p_value        # P-value
result$rate_ratio     # Estimated rate ratio (exp(estimate))
result$conf_int       # Confidence interval for rate ratio
result$dispersion     # Estimated NB dispersion (1/theta)
result$group_summary  # Events, exposure by group

# Poisson fallback
result_poisson <- mutze_test(cut_summary, method = "poisson")

# Two-sided test
result_2sided <- mutze_test(cut_summary, sided = 2, conf_level = 0.95)
```

### Extracting information from the test

```r
# Statistical information = 1 / SE^2
info <- 1 / result$se^2

# Log rate ratio estimate
log_rr <- result$estimate
```

---

## Blinded and unblinded information {#information}

### Blinded information at interim

```r
blinded_res <- calculate_blinded_info(
  data = cut_summary,
  ratio = 1,                         # 1:1 randomization
  lambda1_planning = 0.5,            # Planned control rate
  lambda2_planning = 0.3             # Planned experimental rate
)
blinded_res$blinded_info             # Estimated information
blinded_res$dispersion_blinded       # Estimated dispersion
blinded_res$lambda_blinded           # Blinded overall rate
```

### Method of moments estimation

```r
# Blinded (pooled) estimate
mom_blinded <- estimate_nb_mom(cut_summary)
mom_blinded$lambda      # Overall rate
mom_blinded$dispersion  # Dispersion k

# Unblinded (by group) estimate
mom_unblinded <- estimate_nb_mom(cut_summary, group = "treatment")
mom_unblinded$lambda    # Named vector: Control, Experimental
mom_unblinded$dispersion  # Common dispersion
```

### Information at a planning time

```r
info <- compute_info_at_time(
  analysis_time = 18,
  accrual_rate = 20,
  accrual_duration = 12,
  lambda1 = 0.5, lambda2 = 0.3,
  dispersion = 0.5,
  ratio = 1,
  dropout_rate = 0.05,
  event_gap = 5/365.25,
  max_followup = 12
)
```

---

## Sample size re-estimation {#ssr}

### Blinded SSR (Friede & Schmidli)

```r
blinded_ssr_res <- blinded_ssr(
  data = cut_summary,
  ratio = 1,
  lambda1_planning = 0.5,
  lambda2_planning = 0.3,
  power = 0.9,
  alpha = 0.025,
  accrual_rate = 20,
  accrual_duration = 12,
  trial_duration = 24,
  max_followup = 12,
  event_gap = 5/365.25
)
blinded_ssr_res$n_total_blinded    # Re-estimated total N
blinded_ssr_res$info_fraction      # Current information fraction
```

### Unblinded SSR

```r
unblinded_ssr_res <- unblinded_ssr(
  data = cut_summary,
  ratio = 1,
  lambda1_planning = 0.5,
  lambda2_planning = 0.3,
  power = 0.9,
  alpha = 0.025,
  accrual_rate = 20,
  accrual_duration = 12,
  trial_duration = 24,
  max_followup = 12,
  event_gap = 5/365.25
)
unblinded_ssr_res$n_total_unblinded
unblinded_ssr_res$lambda1_unblinded  # Observed control rate
unblinded_ssr_res$lambda2_unblinded  # Observed experimental rate
```

### Updating group sequential bounds after SSR

```r
# Observed information fraction
info_frac <- blinded_ssr_res$info_fraction

# Update bounds with observed information
gs_update <- gsDesign::gsDesign(
  k = 2,
  test.type = 1,
  alpha = 0.025,
  sfu = gsDesign::sfLDOF,
  timing = info_frac         # Observed information fraction
)
# Use gs_update$upper$bound for updated efficacy bounds
```

---

## Completers-based interim analysis {#completers}

An interim analysis triggered by the number of patients completing full follow-up.

```r
target_completers <- 80

# Find date when 80 patients have completed max_followup
comp_date <- cut_date_for_completers(sim_data, target_completers)

# Cut data at that date
interim_data <- cut_completers(sim_data, cut_date = comp_date, event_gap = 5/365.25)

# Test
result <- mutze_test(interim_data)
z_interim <- result$z
info_interim <- 1 / result$se^2

# Final analysis
final_date <- max(sim_data$calendar_time)
final_data <- cut_data_by_date(sim_data, cut_date = final_date, event_gap = 5/365.25)
result_final <- mutze_test(final_data)
info_final <- 1 / result_final$se^2

# Adaptive bounds based on observed information fraction
gs_adaptive <- gsDesign::gsDesign(
  k = 2, test.type = 1, alpha = 0.025,
  sfu = gsDesign::sfLDOF,
  timing = info_interim / info_final
)
```

---

## Seasonal event simulation {#seasonal}

`nb_sim_seasonal()` simulates events with season-dependent rates.

```r
fail_rate <- data.frame(
  treatment = rep(c("Control", "Experimental"), each = 4),
  season = rep(c("Winter", "Spring", "Summer", "Fall"), 2),
  rate = c(2.0, 0.5, 0.2, 1.5,       # Control by season
           1.4, 0.35, 0.14, 1.05),    # Experimental (30% reduction)
  dispersion = 0.5
)

sim_seasonal <- nb_sim_seasonal(
  enroll_rate = data.frame(rate = 40, duration = 0.5),
  fail_rate = fail_rate,
  dropout_rate = data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.05, 0.05),
    duration = c(100, 100)
  ),
  max_followup = 1,
  randomization_start_date = as.Date("2024-01-01"),
  n = 200
)

# Output: one row per season-interval per patient
# Columns: id, treatment, season, enroll_time, start, end, event,
#           calendar_start, calendar_end

# Cut and analyze
cut_seasonal <- cut_data_by_date(sim_seasonal, cut_date = 1.25, event_gap = 7/365.25)

# Fit NB GLM with seasonal adjustment
fit <- MASS::glm.nb(
  events ~ treatment + season + offset(log(tte)),
  data = cut_seasonal[cut_seasonal$tte > 0, ]
)
rr <- exp(coef(fit)["treatmentExperimental"])
```

---

## Verification of operating characteristics {#verification}

Compare theoretical (formula-based) predictions to simulation results.

```r
# Theoretical predictions from design
design <- sample_size_nbinom(
  lambda1 = 0.4, lambda2 = 0.3, dispersion = 0.5,
  power = 0.9, alpha = 0.025,
  accrual_rate = c(1, 2), accrual_duration = c(6, 6),
  trial_duration = 24, dropout_rate = 0.1/12,
  max_followup = 12, event_gap = 20/30.42
)

# Simulate many trials
sim_results <- sim_gs_nbinom(
  n_sims = 3600,   # SE of power estimate ~ 0.005
  enroll_rate = data.frame(
    rate = design$accrual_rate * design$n_total,
    duration = design$inputs$accrual_duration
  ),
  fail_rate = data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.4, 0.3),
    dispersion = c(0.5, 0.5)
  ),
  dropout_rate = data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.1/12, 0.1/12),
    duration = c(100, 100)
  ),
  max_followup = 12,
  event_gap = 20/30.42,
  analysis_times = design$inputs$trial_duration,
  n_target = design$n_total,
  design = design
)

# Compare theoretical vs. simulated
sim_final <- sim_results[sim_results$analysis == 1, ]
comparison <- data.frame(
  Metric = c("Events (Control)", "Events (Experimental)",
             "Exposure at Risk (Control)", "Power"),
  Theoretical = c(design$events_n1, design$events_n2,
                   design$exposure_at_risk_n1 * design$n1, design$power),
  Simulated = c(mean(sim_final$events_ctrl), mean(sim_final$events_exp),
                 mean(sim_final$exposure_at_risk_ctrl), mean(sim_final$z_stat > qnorm(1 - 0.025)))
)
```

# Code Patterns for wpgsd

**Note**: These patterns target wpgsd from github.com/Merck/wpgsd.

## Table of Contents
1. [Event count tables and correlation](#event-corr)
2. [Generating event tables from datasets](#event-from-data)
3. [Bonferroni bounds (type = 0)](#bonferroni-bounds)
4. [WPGSD bounds with overall alpha spending (type = 2)](#overall-spending)
5. [WPGSD bounds with separate alpha spending (type = 3)](#separate-spending)
6. [WPGSD bounds with fixed alpha spending (type = 1)](#fixed-spending)
7. [Closed testing with observed p-values](#closed-testing)
8. [Sequential p-values with calc_seq_p](#sequential-p)
9. [Overlapping populations example](#overlapping)
10. [Common control example](#common-control)
11. [Manual correlation for k > 2](#manual-corr)

---

## Event count tables and correlation {#event-corr}

The key input to wpgsd is an event count table that encodes the overlap
between hypothesis populations. `generate_corr()` converts this into a
correlation matrix for all test statistics across hypotheses and analyses.

### Event count table structure

The event table has 4 columns: `H1`, `H2`, `Analysis`, `Event`.

- `(H1=i, H2=i)`: events for hypothesis i alone
- `(H1=i, H2=j)` where `i != j`: events in the intersection of populations i and j

```r
library(wpgsd)
library(gsDesign)
library(tibble)

# 3 hypotheses, 2 analyses (IA + FA)
event <- tribble(
  ~H1, ~H2, ~Analysis, ~Event,
  # Individual hypothesis events at IA
  1,   1,   1,         100,     # H1 events at IA
  2,   2,   1,         110,     # H2 events at IA
  3,   3,   1,         225,     # H3 events at IA
  # Pairwise intersection events at IA
  1,   2,   1,         80,      # Events in H1 ∩ H2 population at IA
  1,   3,   1,         100,     # Events in H1 ∩ H3 population at IA
  2,   3,   1,         110,     # Events in H2 ∩ H3 population at IA
  # Individual hypothesis events at FA
  1,   1,   2,         200,
  2,   2,   2,         220,
  3,   3,   2,         450,
  # Pairwise intersection events at FA
  1,   2,   2,         160,
  1,   3,   2,         200,
  2,   3,   2,         220
)
```

### Generating correlation matrix

```r
corr <- generate_corr(event)
# Returns a named matrix with rows/columns like H1_A1, H2_A1, H3_A1, H1_A2, ...
# Dimension: (k * n_hyp) x (k * n_hyp) = 6 x 6 for 3 hypotheses, 2 analyses
round(corr, 2)
```

### Correlation formula

The correlation between test statistics Z_{ik} and Z_{i'k'} is:

```
Corr(Z_{ik}, Z_{i'k'}) = n_{i∧i', k∧k'} / sqrt(n_{ik} * n_{i'k'})
```

Where `n_{i∧i', k∧k'}` is the number of events in the intersection of
populations i and i' at the earlier of analyses k and k'.

---

## Generating event tables from datasets {#event-from-data}

`generate_event_table()` reads ADSL and ADTTE SAS datasets to compute
the event count table automatically.

```r
# Paths to analysis datasets (one per analysis)
paths <- system.file("extdata/", package = "wpgsd")

# Selection criteria for each hypothesis
h_select <- tribble(
  ~Hypothesis, ~Crit,
  1, "PARAMCD == 'OS' & TRT01P %in% c('Xanomeline High Dose', 'Placebo')",
  2, "PARAMCD == 'OS' & TRT01P %in% c('Xanomeline Low Dose', 'Placebo')"
)

result <- generate_event_table(
  paths,
  h_select,
  adsl_name = "adsl",
  adtte_name = "adtte",
  key_var = "USUBJID",
  cnsr_var = "CNSR"
)

event <- result$event   # Event count table for generate_corr/generate_bounds
dsets <- result$dsets    # Analysis datasets for each hypothesis
```

---

## Bonferroni bounds (type = 0) {#bonferroni-bounds}

Type 0 computes standard weighted Bonferroni group sequential bounds
(no correlation adjustment). This serves as the baseline for comparison.

```r
# Multiplicity graph
w <- c(0.3, 0.3, 0.4)   # Initial weights
m <- matrix(c(
  0,   3/7, 4/7,
  3/7, 0,   4/7,
  0.5, 0.5, 0
), nrow = 3, byrow = TRUE)

# Correlation matrix (still needed but not used for bound inflation)
corr <- generate_corr(event)

# Bonferroni bounds: separate spending per hypothesis
bound_bonf <- generate_bounds(
  type = 0,
  k = 2,
  w = w,
  m = m,
  corr = corr,
  alpha = 0.025,
  sf = list(sfHSD, sfHSD, sfHSD),       # One per hypothesis
  sfparm = list(-4, -4, -4),            # One per hypothesis
  t = list(c(0.5, 1), c(0.5, 1), c(0.5, 1))  # Info fractions per hypothesis
)

# Output: tibble with columns Analysis, Hypotheses, H1, H2, H3
# Each row = one intersection hypothesis at one analysis
# Values = nominal p-value bounds
bound_bonf
```

---

## WPGSD bounds with overall alpha spending (type = 2) {#overall-spending}

Type 2 uses a single spending function for the overall intersection
hypothesis, finding `alpha*` that accounts for known correlations.
This is Method 3b in Anderson et al. (2022).

```r
set.seed(1234)

bound_wpgsd <- generate_bounds(
  type = 2,
  k = 2,
  w = w,
  m = m,
  corr = corr,
  alpha = 0.025,
  sf = sfHSD,              # Single spending function
  sfparm = -4,             # Single parameter
  t = c(0.5, 1)            # Single info fraction vector
)

# Compare to Bonferroni: WPGSD bounds are inflated (less conservative)
bound_wpgsd
```

### Using minimum information fraction as spending time

```r
# Use minimum IF across hypotheses for spending time
IF_IA <- c(100/200, 110/220, 225/450)  # IFs for H1, H2, H3

bound_wpgsd <- generate_bounds(
  type = 2,
  k = 2,
  w = w,
  m = m,
  corr = corr,
  alpha = 0.025,
  sf = sfHSD,
  sfparm = -4,
  t = c(min(IF_IA), 1)
)
```

---

## WPGSD bounds with separate alpha spending (type = 3) {#separate-spending}

Type 3 uses separate spending functions per hypothesis, finding an inflation
factor xi that accounts for correlations. This is Method 3c in Anderson
et al. (2022).

```r
set.seed(1234)

# Each hypothesis can have its own spending function and info fraction
bound_wpgsd_3 <- generate_bounds(
  type = 3,
  k = 2,
  w = w,
  m = m,
  corr = corr,
  alpha = 0.025,
  sf = list(sfLDOF, sfLDOF, sfLDOF),
  sfparm = list(0, 0, 0),
  t = list(c(100/200, 1), c(110/220, 1), c(225/450, 1))
)

# Output includes an inflation factor (xi) column
bound_wpgsd_3
```

### Different spending functions per hypothesis

```r
bound_mixed <- generate_bounds(
  type = 3,
  k = 2,
  w = w,
  m = m,
  corr = corr,
  alpha = 0.025,
  sf = list(sfLDOF, sfHSD, sfHSD),
  sfparm = list(0, -4, -2),
  t = list(c(0.5, 1), c(0.5, 1), c(0.5, 1))
)
```

---

## WPGSD bounds with fixed alpha spending (type = 1) {#fixed-spending}

Type 1 uses pre-specified cumulative alpha spending at each analysis.
This is Method 3a in Anderson et al. (2022).

```r
# Specify cumulative alpha spent at each analysis
cum_alpha <- c(0.001, 0.025)  # Spend 0.001 at IA, full 0.025 at FA

bound_fixed <- generate_bounds(
  type = 1,
  k = 2,
  w = w,
  m = m,
  corr = corr,
  alpha = 0.025,
  cum_alpha = cum_alpha
)
bound_fixed
```

---

## Closed testing with observed p-values {#closed-testing}

After computing bounds, use `closed_test()` to test observed p-values
against the bounds using the closed testing principle.

```r
# Observed p-values at each analysis
p_obs <- tribble(
  ~Analysis, ~H1,    ~H2,     ~H3,
  1,         0.01,   0.0004,  0.03,
  2,         0.05,   0.002,   0.015
)

# Closed testing
result <- closed_test(bound_wpgsd, p_obs)
result
# Output: matrix showing which hypotheses are rejected
# A hypothesis is rejected only if ALL intersection hypotheses containing it
# are rejected (closed testing principle)
```

### Interpreting closed_test output

The output matrix has rows for each elementary hypothesis and columns
indicating rejection status. A hypothesis H_i is rejected if and only if
every intersection hypothesis containing H_i has its observed p-value
below the corresponding bound.

---

## Sequential p-values with calc_seq_p {#sequential-p}

`calc_seq_p()` computes sequential p-values for each intersection hypothesis
at a given analysis. These can be used with graphical testing procedures.

### WPGSD sequential p-values (type = 2)

```r
# Observed p-values (must include all analyses, even future)
p_obs <- tibble(
  analysis = 1:2,
  H1 = c(0.02, 0.015),
  H2 = c(0.01, 0.012),
  H3 = c(0.012, 0.010)
)

# Sequential p-value for full intersection {H1, H2, H3} at IA
seq_p <- calc_seq_p(
  test_analysis = 1,
  test_hypothesis = "H1, H2, H3",
  p_obs = p_obs,
  alpha_spending_type = 2,           # Overall spending (type 2)
  n_analysis = 2,
  initial_weight = w,
  transition_mat = m,
  z_corr = corr,
  spending_fun = sfHSD,
  spending_fun_par = -4,
  info_frac = c(0.5, 1),
  interval = c(1e-4, 0.2)
)

# Sequential p-value for subset {H1, H2} at IA
seq_p_12 <- calc_seq_p(
  test_analysis = 1,
  test_hypothesis = "H1, H2",
  p_obs = p_obs,
  alpha_spending_type = 2,
  n_analysis = 2,
  initial_weight = w,
  transition_mat = m,
  z_corr = corr,
  spending_fun = sfHSD,
  spending_fun_par = -4,
  info_frac = c(0.5, 1),
  interval = c(1e-4, 0.2)
)

# Sequential p-value for elementary {H1} at IA
seq_p_1 <- calc_seq_p(
  test_analysis = 1,
  test_hypothesis = "H1",
  p_obs = p_obs,
  alpha_spending_type = 2,
  n_analysis = 2,
  initial_weight = w,
  transition_mat = m,
  z_corr = corr,
  spending_fun = sfHSD,
  spending_fun_par = -4,
  info_frac = c(0.5, 1),
  interval = c(1e-4, 0.2)
)
```

### Bonferroni sequential p-values (type = 0)

```r
# For Bonferroni, spending functions are lists (one per hypothesis)
seq_p_bonf <- calc_seq_p(
  test_analysis = 1,
  test_hypothesis = "H1, H2, H3",
  p_obs = p_obs,
  alpha_spending_type = 0,
  n_analysis = 2,
  initial_weight = w,
  transition_mat = m,
  z_corr = corr,
  spending_fun = list(sfHSD, sfHSD, sfHSD),
  spending_fun_par = list(-4, -4, -4),
  info_frac = list(c(0.5, 1), c(0.5, 1), c(0.5, 1)),
  interval = c(1e-4, 0.3)
)
```

### Sequential p-values at final analysis

```r
seq_p_fa <- calc_seq_p(
  test_analysis = 2,                 # Final analysis

  test_hypothesis = "H1, H2, H3",
  p_obs = p_obs,
  alpha_spending_type = 2,
  n_analysis = 2,
  initial_weight = w,
  transition_mat = m,
  z_corr = corr,
  spending_fun = sfHSD,
  spending_fun_par = -4,
  info_frac = c(0.5, 1),
  interval = c(1e-4, 0.2)
)
```

---

## Overlapping populations example {#overlapping}

Three hypotheses from overlapping biomarker populations (e.g., biomarker
A+, biomarker B+, overall). Correlation arises from shared patients.

```r
library(wpgsd)
library(gsDesign)
library(tibble)
library(dplyr)

# Multiplicity graph
w <- c(0.3, 0.3, 0.4)  # H1: BM-A+, H2: BM-B+, H3: Overall
m <- matrix(c(
  0,   3/7, 4/7,
  3/7, 0,   4/7,
  0.5, 0.5, 0
), nrow = 3, byrow = TRUE)

# Event counts at IA and FA
# Key: (i, j) intersection means patients in BOTH population i and j
# For nested: H1 ⊂ H3, so n_{1∧3} = n_1 (all H1 patients are in H3)
event <- tribble(
  ~H1, ~H2, ~Analysis, ~Event,
  # IA
  1,   1,   1,         100,    # BM-A+ events
  2,   2,   1,         110,    # BM-B+ events
  3,   3,   1,         225,    # Overall events
  1,   2,   1,         80,     # BM-A+ ∩ BM-B+ (AB positive)
  1,   3,   1,         100,    # BM-A+ ∩ Overall = BM-A+
  2,   3,   1,         110,    # BM-B+ ∩ Overall = BM-B+
  # FA
  1,   1,   2,         200,
  2,   2,   2,         220,
  3,   3,   2,         450,
  1,   2,   2,         160,
  1,   3,   2,         200,
  2,   3,   2,         220
)

corr <- generate_corr(event)

# Information fractions
IF_IA <- c(100/200, 110/220, 225/450)

# WPGSD bounds (overall spending)
set.seed(1234)
bound <- generate_bounds(
  type = 2, k = 2, w = w, m = m,
  corr = corr, alpha = 0.025,
  sf = sfHSD, sfparm = -4,
  t = c(min(IF_IA), 1)
)

# Observed p-values
p_obs <- tribble(
  ~Analysis, ~H1,   ~H2,    ~H3,
  1,         0.01,  0.0004, 0.03,
  2,         0.05,  0.002,  0.015
)

# Closed testing
result <- closed_test(bound, p_obs)
result
```

---

## Common control example {#common-control}

Three experimental arms vs a common control. Correlation arises from
shared control arm events.

```r
# Holm-like equal-weight graph
w <- c(1/3, 1/3, 1/3)
m <- matrix(c(
  0,   0.5, 0.5,
  0.5, 0,   0.5,
  0.5, 0.5, 0
), nrow = 3, byrow = TRUE)

# Events per arm at IA and FA
# H_i: Experimental i vs Control
# n_{i,i}: events in Exp_i + Control
# n_{i,j}: events in Control only (shared)
event <- tribble(
  ~H1, ~H2, ~Analysis, ~Event,
  # IA
  1,   1,   1,         155,    # Exp1 + Control events
  2,   2,   1,         160,    # Exp2 + Control events
  3,   3,   1,         165,    # Exp3 + Control events
  1,   2,   1,         85,     # Control events (shared by H1, H2)
  1,   3,   1,         85,     # Control events (shared by H1, H3)
  2,   3,   1,         85,     # Control events (shared by H2, H3)
  # FA
  1,   1,   2,         305,
  2,   2,   2,         320,
  3,   3,   2,         335,
  1,   2,   2,         170,
  1,   3,   2,         170,
  2,   3,   2,         170
)

corr <- generate_corr(event)

IF_IA <- c(155/305, 160/320, 165/335)

# WPGSD bounds (separate spending, method 3c)
set.seed(1234)
bound <- generate_bounds(
  type = 3, k = 2, w = w, m = m,
  corr = corr, alpha = 0.025,
  sf = list(sfLDOF, sfLDOF, sfLDOF),
  sfparm = list(0, 0, 0),
  t = list(c(IF_IA[1], 1), c(IF_IA[2], 1), c(IF_IA[3], 1))
)

bound
```

---

## Manual correlation for k > 2 {#manual-corr}

`generate_corr()` has a known bug for k > 2 analyses: it incorrectly
computes within-hypothesis cross-analysis entries for non-adjacent analyses.
For k > 2, build the correlation matrix manually.

### Formula

```
Corr(Z_{ik}, Z_{i'k'}) = n_{i∧i', k∧k'} / sqrt(n_{ik} * n_{i'k'})
```

Where `k∧k' = min(k, k')` means events at the earlier analysis.

### Manual construction for 2 hypotheses, 3 analyses

```r
# Event counts
# n[i, k]: events for hypothesis i at analysis k
# n_int[i, j, k]: intersection events at analysis k
n <- matrix(c(
  100, 200, 300,    # H1 at analyses 1, 2, 3
  120, 240, 360     # H2 at analyses 1, 2, 3
), nrow = 2, byrow = TRUE)

n_int <- array(0, dim = c(2, 2, 3))
# Intersection H1∩H2 at each analysis
n_int[1, 2, ] <- n_int[2, 1, ] <- c(60, 120, 180)
# Diagonal = individual counts
n_int[1, 1, ] <- n[1, ]
n_int[2, 2, ] <- n[2, ]

# Build 6x6 correlation matrix (2 hyps x 3 analyses)
n_hyp <- 2
n_anal <- 3
dim <- n_hyp * n_anal
D <- matrix(0, dim, dim)

for (i in 1:n_hyp) {
  for (ip in 1:n_hyp) {
    for (k in 1:n_anal) {
      for (kp in 1:n_anal) {
        row <- (k - 1) * n_hyp + i
        col <- (kp - 1) * n_hyp + ip
        kmin <- min(k, kp)
        D[row, col] <- n_int[i, ip, kmin]
      }
    }
  }
}

# Convert to correlation
corr_manual <- diag(1 / sqrt(diag(D))) %*% D %*% diag(1 / sqrt(diag(D)))

# Add names
labs <- paste0("H", rep(1:n_hyp, n_anal), "_A", rep(1:n_anal, each = n_hyp))
rownames(corr_manual) <- colnames(corr_manual) <- labs
round(corr_manual, 3)
```

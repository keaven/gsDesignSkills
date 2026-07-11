# Code Patterns for the Illness-Death Model

## Table of Contents
1. [Building transition rates](#build-rates)
2. [Simulating a trial](#simulate)
3. [Analysis cut dates](#cut-dates)
4. [Cutting data to ADTTE format](#cut-data)
5. [Nominal p-values (logrank and risk difference)](#p-values)
6. [Theoretical survival curves](#theoretical-curves)
7. [Overall population curves (prevalence-weighted)](#overall-curves)
8. [Full workflow: design assumptions to analysis](#full-workflow)

---

## Building transition rates {#build-rates}

`build_transition_rates()` calibrates illness-death transition rates from
clinically interpretable parameters. The model has states:
- 0: Alive, no response, no progression
- 1: Responded, no progression
- 2: Progressed
- 3: Dead (absorbing)

### Basic usage

```r
source("inst/simulation/sim_illness_death.R")
source("inst/simulation/cut_illness_death.R")

transition_rate <- build_transition_rates(
  strata = c("BM+", "BM-"),
  treatments = c("control", "experimental"),
  median_pfs = c("BM+" = 5, "BM-" = 5),       # Control median PFS (months)
  median_os = c("BM+" = 12, "BM-" = 12),       # Control median OS (months)
  orr = list(
    "BM+" = c(control = 0.15, experimental = 0.35),
    "BM-" = c(control = 0.15, experimental = 0.25)
  ),
  hr_pfs = c("BM+" = 0.65, "BM-" = 0.85),     # Experimental/control HR
  hr_os = c("BM+" = 0.65, "BM-" = 0.85)
)
```

### Output format

The result is a data frame with columns: `stratum`, `treatment`, `transition`,
`rate`, `duration`. Six transitions per (stratum, treatment) combination:

| Transition | From | To | Description |
|-----------|------|-----|-------------|
| `response` | 0 | 1 | Response (before progression/death) |
| `prog_0` | 0 | 2 | Progression without prior response |
| `death_0` | 0 | 3 | Death without prior progression |
| `prog_1` | 1 | 2 | Progression after response |
| `death_1` | 1 | 3 | Death after response (without progression) |
| `death_2` | 2 | 3 | Death after progression |

### Calibration details

Control arm rates are calibrated via `uniroot()`:
1. `lambda_0` (total exit rate from state 0) is found such that the marginal PFS CDF equals 0.5 at the specified median PFS.
2. `lambda_death2` (post-progression death rate) is found such that the marginal OS CDF equals 0.5 at the specified median OS.
3. Response rate `lambda_resp = ORR * lambda_0` ensures the competing risks probability of response matches the specified ORR.

Experimental arm rates:
1. Progression rates scaled by `hr_pfs`: `lambda_prog0_exp = lambda_prog0_ctrl * hr_pfs`
2. Death rates scaled by `hr_os`: `lambda_death0_exp = lambda_death0_ctrl * hr_os`
3. Response rate solved to match experimental ORR exactly.

### Advanced parameters

```r
transition_rate <- build_transition_rates(
  strata = "All",
  treatments = c("control", "experimental"),
  median_pfs = c(All = 5),
  median_os = c(All = 12),
  orr = list(All = c(control = 0.30, experimental = 0.45)),
  hr_pfs = c(All = 0.70),
  hr_os = c(All = 0.65),
  death_wo_prog_rate = 0.02,     # Rate of death without progression (default 0.02)
  responder_prog_ratio = 0.5     # Responders progress at 50% of non-responder rate
)
```

### Weaker simulation effects (for testing scenarios)

To simulate under effects weaker than design assumptions (e.g., to illustrate
a scenario where not all hypotheses are rejected):

```r
hr_os_sim_bm_pos <- 0.70    # Design assumed 0.65
hr_pfs_sim_bm_pos <- 0.65   # Design assumed 0.65
hr_os_sim_bm_neg <- 1.1     # Design assumed 0.85 (reversed effect in BM-)
hr_pfs_sim_bm_neg <- 1.2    # Design assumed 0.85 (reversed effect in BM-)
orr_exp_sim_bm_pos <- 0.30  # Design assumed 0.35
orr_exp_sim_bm_neg <- 0.12  # Design assumed 0.25

transition_rate <- build_transition_rates(
  strata = c("BM+", "BM-"),
  treatments = c("control", "experimental"),
  median_pfs = c("BM+" = 5, "BM-" = 5),
  median_os = c("BM+" = 12, "BM-" = 12),
  orr = list(
    "BM+" = c(control = 0.15, experimental = orr_exp_sim_bm_pos),
    "BM-" = c(control = 0.15, experimental = orr_exp_sim_bm_neg)
  ),
  hr_pfs = c("BM+" = hr_pfs_sim_bm_pos, "BM-" = hr_pfs_sim_bm_neg),
  hr_os = c("BM+" = hr_os_sim_bm_pos, "BM-" = hr_os_sim_bm_neg)
)
```

---

## Simulating a trial {#simulate}

`sim_illness_death()` generates patient-level data with enrollment, randomization,
and multi-state outcomes.

```r
library(simtrial)

set.seed(54321)
n_total <- 500
steady_rate <- n_total / 15  # 18-month ramp: 0.25*2 + 0.5*2 + 0.75*2 + 1*12 = 15 rate-months

sim_data <- sim_illness_death(
  n = n_total,
  stratum = data.frame(stratum = c("BM+", "BM-"), p = c(0.50, 0.50)),
  block = c("control", "control", "experimental", "experimental"),
  enroll_rate = data.frame(
    rate = steady_rate * c(0.25, 0.50, 0.75, 1.00),
    duration = c(2, 2, 2, 12)
  ),
  transition_rate = transition_rate,
  dropout_rate = 0.001
)

fpe <- max(sim_data$ENRLTIME)
cat("Total N:", nrow(sim_data), "\n")
cat("N BM+:", sum(sim_data$STRATUM == "BM+"),
    "N BM-:", sum(sim_data$STRATUM == "BM-"), "\n")
cat("FPE:", round(fpe, 1), "months\n")
```

### Output columns

| Column | Description |
|--------|-------------|
| `USUBJID` | Patient ID |
| `STRATUM` | Stratum assignment |
| `TRT` | Treatment ("control" or "experimental") |
| `ENRLTIME` | Enrollment time (months from study start) |
| `OS_TIME` | Time from randomization to death |
| `PFS_TIME` | Time from randomization to progression or death |
| `TTR` | Time to response (Inf if no response) |
| `ORR` | 1 if responded, 0 otherwise |
| `DROPOUT_TIME` | Time to dropout (independent exponential) |
| `CTE_OS` | Calendar time of OS event (`ENRLTIME + min(OS_TIME, DROPOUT_TIME)`) |
| `CTE_PFS` | Calendar time of PFS event |

### Single-stratum simulation

```r
sim_data <- sim_illness_death(
  n = 400,
  stratum = data.frame(stratum = "All", p = 1),
  block = c("control", "control", "experimental", "experimental"),
  enroll_rate = data.frame(rate = 30, duration = 14),
  transition_rate = transition_rate
)
```

---

## Analysis cut dates {#cut-dates}

`get_analysis_dates()` computes calendar cut dates from timing rules that combine
minimum follow-up, event-driven triggers, and maximum extension caps.

```r
# Target event counts from group sequential designs
target_pfs_events_bm <- ceiling(max(pfssub$analysis$event))
target_os_events_bm <- ceiling(max(ossub$analysis$event))

analyses <- list(
  # IA1: 6 months after FPE (calendar-only, no event trigger)
  list(min_followup = 6, endpoint = NULL, event_target = NULL,
       target_stratum = NULL, max_followup = NULL),
  # IA2: 14 months after FPE OR targeted PFS events in BM+, capped at 17 months
  list(min_followup = 14, endpoint = "PFS", event_target = target_pfs_events_bm,
       target_stratum = "BM+", max_followup = 17),
  # Final: 24 months after FPE OR targeted OS events in BM+, capped at 30 months
  list(min_followup = 24, endpoint = "OS", event_target = target_os_events_bm,
       target_stratum = "BM+", max_followup = 30)
)

cut_dates <- get_analysis_dates(sim_data, analyses)

for (i in seq_along(cut_dates)) {
  cat(sprintf("Analysis %d: %.1f months (FPE + %.1f)\n",
              i, cut_dates[i], cut_dates[i] - fpe))
}
```

### Timing logic

For each analysis:
1. `calendar_min = FPE + min_followup`
2. If event-driven: find the calendar time when `event_target` events accumulate in `target_stratum`
3. `cut_date = max(calendar_min, event_date)` — both conditions must be met
4. If `max_followup` set: `cut_date = min(cut_date, FPE + max_followup)` — cap the delay

---

## Cutting data to ADTTE format {#cut-data}

`cut_illness_death()` applies administrative censoring and produces long-format
ADTTE datasets.

```r
adtte <- lapply(seq_along(cut_dates), function(i) {
  d <- cut_illness_death(sim_data, cut_dates[i])
  d$ANALYSIS <- i
  d
})
names(adtte) <- paste0("IA", seq_along(cut_dates))
```

### Output format

Each row represents one patient-endpoint combination:

| Column | Description |
|--------|-------------|
| `USUBJID` | Patient ID |
| `STRATUM` | Stratum |
| `TRT` | Treatment |
| `PARAMCD` | Endpoint code: "OS", "PFS", "TTR", "ORR" |
| `PARAM` | Endpoint name |
| `AVAL` | Analysis value (months for TTE; 0/1 for ORR) |
| `AVALU` | Units |
| `CNSR` | Censoring indicator (0 = event, 1 = censored; NA for ORR) |
| `EVNTDESC` | Event description |

### Selecting endpoints

```r
# Only OS and PFS (no ORR/TTR)
adtte_tte <- cut_illness_death(sim_data, cut_date = 36, paramcd = c("OS", "PFS"))

# Only ORR
adtte_orr <- cut_illness_death(sim_data, cut_date = 24, paramcd = "ORR")
```

### Event counts summary

```r
event_summary <- do.call(rbind, lapply(seq_along(adtte), function(i) {
  endpoints <- if (i == 1) c("OS", "PFS", "ORR") else if (i == 2) c("OS", "PFS") else "OS"
  adtte[[i]] %>%
    filter(PARAMCD %in% endpoints) %>%
    group_by(PARAMCD, STRATUM) %>%
    summarise(
      count = if (first(PARAMCD) == "ORR") n() else sum(CNSR == 0),
      .groups = "drop"
    ) %>%
    mutate(Analysis = i)
}))
```

---

## Nominal p-values (logrank and risk difference) {#p-values}

### One-sided logrank test (OS, PFS)

```r
logrank_pval <- function(data, stratified = FALSE) {
  data$event <- 1L - data$CNSR
  if (stratified) {
    fit <- survdiff(Surv(AVAL, event) ~ TRT + strata(STRATUM), data = data)
  } else {
    fit <- survdiff(Surv(AVAL, event) ~ TRT, data = data)
  }
  # IMPORTANT: with strata(), obs and exp are matrices (groups x strata)
  # Must sum across strata columns
  obs_ctrl <- if (is.matrix(fit$obs)) sum(fit$obs[1, ]) else fit$obs[1]
  exp_ctrl <- if (is.matrix(fit$exp)) sum(fit$exp[1, ]) else fit$exp[1]
  z_sign <- sign(obs_ctrl - exp_ctrl)
  z <- z_sign * sqrt(fit$chisq)
  pnorm(-z)  # 1-sided: small when experimental is better
}
```

**Key pitfall**: `survdiff()` with `strata()` returns `obs` and `exp` as
matrices (treatment groups x strata). Using `fit$obs[1]` extracts just the
first stratum's control group value. Always check `is.matrix()` and sum.

### Stratified risk difference test (ORR)

```r
rd_pval <- function(data, stratified = FALSE) {
  if (!stratified) {
    n_exp <- sum(data$TRT == "experimental")
    n_ctrl <- sum(data$TRT == "control")
    resp_exp <- sum(data$AVAL[data$TRT == "experimental"])
    resp_ctrl <- sum(data$AVAL[data$TRT == "control"])
    p_exp <- resp_exp / n_exp
    p_ctrl <- resp_ctrl / n_ctrl
    p_pool <- (resp_exp + resp_ctrl) / (n_exp + n_ctrl)
    se <- sqrt(p_pool * (1 - p_pool) * (1 / n_exp + 1 / n_ctrl))
    if (se == 0) return(0.5)
    z <- (p_exp - p_ctrl) / se
    return(pnorm(-z))
  }

  # Stratified: combine stratum-specific RDs with sample size weights
  strata <- unique(data$STRATUM)
  n_s <- rd_s <- var_s <- numeric(length(strata))
  for (k in seq_along(strata)) {
    ds <- data[data$STRATUM == strata[k], ]
    n_exp <- sum(ds$TRT == "experimental")
    n_ctrl <- sum(ds$TRT == "control")
    n_s[k] <- n_exp + n_ctrl
    p_exp <- sum(ds$AVAL[ds$TRT == "experimental"]) / n_exp
    p_ctrl <- sum(ds$AVAL[ds$TRT == "control"]) / n_ctrl
    rd_s[k] <- p_exp - p_ctrl
    var_s[k] <- p_exp * (1 - p_exp) / n_exp + p_ctrl * (1 - p_ctrl) / n_ctrl
  }
  w <- n_s / sum(n_s)
  rd_combined <- sum(w * rd_s)
  se_combined <- sqrt(sum(w^2 * var_s))
  if (se_combined == 0) return(0.5)
  z <- rd_combined / se_combined
  pnorm(-z)
}
```

### Computing p-values across analyses

```r
pvals <- list()
for (i in seq_along(cut_dates)) {
  d <- adtte[[i]]
  # BM+ subgroup (unstratified within subgroup)
  pvals[[paste0("H1_", i)]] <- logrank_pval(d %>% filter(PARAMCD == "OS", STRATUM == "BM+"))
  # Overall population (stratified by BM+/BM-)
  pvals[[paste0("H2_", i)]] <- logrank_pval(d %>% filter(PARAMCD == "OS"), stratified = TRUE)
  if (i <= 2) {
    pvals[[paste0("H3_", i)]] <- logrank_pval(d %>% filter(PARAMCD == "PFS", STRATUM == "BM+"))
    pvals[[paste0("H4_", i)]] <- logrank_pval(d %>% filter(PARAMCD == "PFS"), stratified = TRUE)
  }
  if (i == 1) {
    pvals[[paste0("H5_", i)]] <- rd_pval(d %>% filter(PARAMCD == "ORR", STRATUM == "BM+"))
    pvals[[paste0("H6_", i)]] <- rd_pval(d %>% filter(PARAMCD == "ORR"), stratified = TRUE)
  }
}
```

---

## Theoretical survival curves {#theoretical-curves}

The internal CDF helpers (`.pfs_cdf()`, `.os_cdf()`) compute marginal PFS and
OS survival curves analytically from the transition rates, without simulation.

```r
get_theoretical_surv <- function(transition_rate, stratum_val, trt_val, t_grid) {
  tr <- transition_rate[transition_rate$stratum == stratum_val &
                          transition_rate$treatment == trt_val, ]
  get_rate <- function(trans) tr$rate[tr$transition == trans]

  lambda_resp <- get_rate("response")
  lambda_prog0 <- get_rate("prog_0")
  lambda_death0 <- get_rate("death_0")
  lambda_prog1 <- get_rate("prog_1")
  lambda_death1 <- get_rate("death_1")
  lambda_death2 <- get_rate("death_2")

  lambda_0 <- lambda_resp + lambda_prog0 + lambda_death0
  orr_param <- lambda_resp / lambda_0
  lambda_1 <- lambda_prog1 + lambda_death1

  pfs_surv <- 1 - sapply(t_grid, function(t)
    .pfs_cdf(t, lambda_0, orr_param, lambda_1))
  os_surv <- 1 - sapply(t_grid, function(t)
    .os_cdf(t, lambda_0, orr_param, lambda_1,
            lambda_death0, lambda_prog0, lambda_death1,
            lambda_prog1, lambda_death2))

  tibble(time = rep(t_grid, 2), surv = c(pfs_surv, os_surv),
         endpoint = rep(c("PFS", "OS"), each = length(t_grid)),
         stratum = stratum_val, treatment = trt_val)
}

# Compute curves for each stratum and treatment
t_grid <- seq(0, 42, by = 0.25)
curves <- bind_rows(
  get_theoretical_surv(transition_rate, "BM+", "control", t_grid),
  get_theoretical_surv(transition_rate, "BM+", "experimental", t_grid),
  get_theoretical_surv(transition_rate, "BM-", "control", t_grid),
  get_theoretical_surv(transition_rate, "BM-", "experimental", t_grid)
)

ggplot(curves, aes(x = time, y = surv, color = treatment, linetype = stratum)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ endpoint) +
  labs(y = "Survival probability", x = "Time (months)")
```

---

## Overall population curves (prevalence-weighted) {#overall-curves}

The overall population survival curve is a prevalence-weighted mixture of
stratum-specific curves:

```r
bm_prev <- 0.50  # BM+ prevalence

bm_pos_curves <- bind_rows(
  get_theoretical_surv(transition_rate, "BM+", "control", t_grid),
  get_theoretical_surv(transition_rate, "BM+", "experimental", t_grid)
)
bm_neg_curves <- bind_rows(
  get_theoretical_surv(transition_rate, "BM-", "control", t_grid),
  get_theoretical_surv(transition_rate, "BM-", "experimental", t_grid)
)

overall_curves <- bm_pos_curves %>%
  select(time, endpoint, treatment, surv_pos = surv) %>%
  left_join(
    bm_neg_curves %>% select(time, endpoint, treatment, surv_neg = surv),
    by = c("time", "endpoint", "treatment")
  ) %>%
  mutate(
    surv = bm_prev * surv_pos + (1 - bm_prev) * surv_neg,
    stratum = "Overall"
  ) %>%
  select(time, surv, endpoint, stratum, treatment)

# Combine by-stratum and overall for a 2x2 faceted plot
all_curves <- bind_rows(
  curves %>% mutate(population = "By stratum"),
  overall_curves %>% mutate(population = "Overall")
)

ggplot(all_curves, aes(x = time, y = surv,
                        color = treatment, linetype = stratum)) +
  geom_line(linewidth = 0.8) +
  facet_grid(population ~ endpoint) +
  labs(y = "Survival probability", x = "Time (months)")
```

---

## Full workflow: design assumptions to analysis {#full-workflow}

This shows the complete pipeline from design assumptions through p-value computation.

```r
library(simtrial)
library(survival)
library(dplyr)
library(tibble)

source("inst/simulation/sim_illness_death.R")
source("inst/simulation/cut_illness_death.R")

# 1. Design assumptions
prevalence_bm <- 0.50
osmedian <- 12
pfsmedian <- 5

# 2. Build transition rates
transition_rate <- build_transition_rates(
  strata = c("BM+", "BM-"),
  treatments = c("control", "experimental"),
  median_pfs = c("BM+" = pfsmedian, "BM-" = pfsmedian),
  median_os = c("BM+" = osmedian, "BM-" = osmedian),
  orr = list(
    "BM+" = c(control = 0.15, experimental = 0.35),
    "BM-" = c(control = 0.15, experimental = 0.25)
  ),
  hr_pfs = c("BM+" = 0.65, "BM-" = 0.85),
  hr_os = c("BM+" = 0.65, "BM-" = 0.85)
)

# 3. Simulate trial
set.seed(12345)
n_total <- 500
sim_data <- sim_illness_death(
  n = n_total,
  stratum = data.frame(stratum = c("BM+", "BM-"), p = c(0.50, 0.50)),
  block = c("control", "control", "experimental", "experimental"),
  enroll_rate = data.frame(
    rate = (n_total / 15) * c(0.25, 0.50, 0.75, 1.00),
    duration = c(2, 2, 2, 12)
  ),
  transition_rate = transition_rate,
  dropout_rate = 0.001
)

# 4. Determine cut dates
analyses <- list(
  list(min_followup = 6, endpoint = NULL, event_target = NULL,
       target_stratum = NULL, max_followup = NULL),
  list(min_followup = 14, endpoint = "PFS", event_target = 250,
       target_stratum = "BM+", max_followup = 17),
  list(min_followup = 24, endpoint = "OS", event_target = 200,
       target_stratum = "BM+", max_followup = 30)
)
cut_dates <- get_analysis_dates(sim_data, analyses)

# 5. Cut data
adtte <- lapply(seq_along(cut_dates), function(i) {
  d <- cut_illness_death(sim_data, cut_dates[i])
  d$ANALYSIS <- i
  d
})

# 6. Compute p-values
# H1: OS BM+ (logrank, unstratified)
p_H1 <- logrank_pval(adtte[[3]] %>% filter(PARAMCD == "OS", STRATUM == "BM+"))
# H2: OS All (logrank, stratified)
p_H2 <- logrank_pval(adtte[[3]] %>% filter(PARAMCD == "OS"), stratified = TRUE)
# H5: ORR BM+ (risk difference, unstratified)
p_H5 <- rd_pval(adtte[[1]] %>% filter(PARAMCD == "ORR", STRATUM == "BM+"))
```

---

## Piecewise transition rate modification {#piecewise-rates}

`build_transition_rates()` produces constant rates. To implement piecewise
hazards (e.g., progression rate that halves after 8 months, or response rate
that varies over time), modify the output data frame by splitting rows into
multiple time periods.

### Piecewise progression rates

Model PFS with median 8 months for $t < 8$, median 16 months for $t \geq 8$:

```r
transition_rate <- build_transition_rates(
  strata = "All",
  treatments = c("control", "experimental"),
  median_pfs = c(All = 8),     # Calibrate with first-period median

  median_os = c(All = 14),
  orr = list(All = c(control = 0.15, experimental = 0.25)),
  hr_pfs = c(All = 0.75),
  hr_os = c(All = 0.75)
)

# Split prog_0 and prog_1 into two periods
prog_transitions <- c("prog_0", "prog_1")
pw_rows <- list()
for (i in seq_len(nrow(transition_rate))) {
  row <- transition_rate[i, ]
  if (row$transition %in% prog_transitions) {
    # Period 1: original rate, duration = changepoint
    row1 <- row; row1$duration <- 8
    pw_rows[[length(pw_rows) + 1]] <- row1
    # Period 2: half the rate, duration = Inf
    row2 <- row; row2$rate <- row$rate / 2; row2$duration <- Inf
    pw_rows[[length(pw_rows) + 1]] <- row2
  } else {
    pw_rows[[length(pw_rows) + 1]] <- row
  }
}
transition_rate_pw <- do.call(rbind, pw_rows)
```

### Piecewise response rates

Model early high response rate (first 6 months) with low late response:

```r
orr_changepoint <- 6
orr_rate_mult <- 2.25     # Multiply base rate for t < 6
orr_rate_late <- 0.1      # 10% of base rate for t >= 6

# Add to the piecewise loop above
if (row$transition == "response") {
  row1 <- row
  row1$rate <- row$rate * orr_rate_mult
  row1$duration <- orr_changepoint
  pw_rows[[length(pw_rows) + 1]] <- row1
  row2 <- row
  row2$rate <- row$rate * orr_rate_late
  row2$duration <- Inf
  pw_rows[[length(pw_rows) + 1]] <- row2
}
```

### Separate design vs simulation transition rates

Use different HRs for design (optimistic) and simulation (weaker effect):

```r
# Design rates: HR = 0.75 (for sample size and bounds)
transition_rate_design <- build_transition_rates(
  ..., hr_pfs = c(All = 0.75), hr_os = c(All = 0.75)
)

# Simulation rates: HR = 0.80 (weaker, for operating characteristics)
transition_rate_sim <- build_transition_rates(
  ..., hr_pfs = c(All = 0.80), hr_os = c(All = 0.80)
)
# Apply the same piecewise modifications to both
```

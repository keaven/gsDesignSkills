---
name: illness-death
description: >
  Guide users through simulating clinical trials using the illness-death model
  with response. Use this skill when the user asks about: multi-state models for
  oncology trials, simulating correlated OS/PFS/ORR endpoints, transition rates
  between disease states, illness-death model calibration, building ADTTE datasets
  from simulation, analysis cut date determination, or theoretical survival curves
  from transition rates.
---

# Illness-Death Model with Response

This skill covers the illness-death model implementation in `gMCPLite/inst/simulation/`,
which simulates correlated OS, PFS, and ORR endpoints for oncology trial design.

## Source files

The model lives in `inst/simulation/` within the gMCPLite package:

- `sim_illness_death.R` — `build_transition_rates()`, `sim_illness_death()`, and internal CDF helpers (`.pfs_cdf()`, `.os_cdf()`)
- `cut_illness_death.R` — `get_analysis_dates()`, `cut_illness_death()`

Load via `source()`:
```r
source(file.path(rprojroot::find_root("DESCRIPTION"), "inst", "simulation", "sim_illness_death.R"))
source(file.path(rprojroot::find_root("DESCRIPTION"), "inst", "simulation", "cut_illness_death.R"))
```

## Model overview

The illness-death model has four states with piecewise exponential transition times:

```
  State 0 ─── response ───> State 1 (responded)
  (alive,     prog_0 ────> State 2 (progressed)
   no resp)   death_0 ───> State 3 (dead)

  State 1 ─── prog_1 ────> State 2
  (responded)  death_1 ───> State 3

  State 2 ─── death_2 ───> State 3
  (progressed)
```

### Derived endpoints
- **PFS**: Time from randomization to progression or death (whichever first)
- **OS**: Time from randomization to death
- **ORR**: Binary indicator of response before progression/death
- **TTR**: Time to response

### Key property: endpoint correlation
Because PFS and OS share the same underlying state transitions, they are
automatically correlated — a patient who progresses early tends to die earlier.
This makes the model more realistic than independently simulating PFS and OS.

## Key functions

### `build_transition_rates()`
Calibrates transition rates from clinically interpretable inputs (median PFS/OS,
ORR, hazard ratios). Uses numerical root-finding so that control arm marginal
PFS and OS medians match the inputs exactly.

### `sim_illness_death()`
Simulates enrollment, randomization, and multi-state outcomes. Uses
`simtrial::rpwexp()` for piecewise exponential transition times and
`simtrial::rpwexp_enroll()` for enrollment.

### `get_analysis_dates()`
Computes calendar cut dates from timing rules: minimum follow-up after FPE,
event-driven triggers with caps, and maximum extensions.

### `cut_illness_death()`
Applies administrative censoring at a cut date and produces ADTTE-format datasets
with OS, PFS, TTR, and ORR endpoints.

### Internal CDF helpers
- `.pfs_cdf(t, ...)` — Marginal PFS CDF from transition rates
- `.os_cdf(t, ...)` — Marginal OS CDF from transition rates
- `.conv2_cdf(t, a, b)` — CDF of sum of two independent exponentials
- `.conv3_cdf(t, a, b, c)` — CDF of sum of three independent exponentials

These allow computing theoretical survival curves without simulation.

## Workflow patterns

For detailed code templates, read `references/code_patterns.md`.

Topics covered:
- Building transition rates from clinical assumptions
- Simulating a single trial
- Analysis cut date determination (calendar + event-driven)
- Cutting data to ADTTE format
- Computing nominal p-values (logrank, stratified logrank, risk difference)
- Piecewise transition rate modification (progression, response)
- Separate design vs simulation effect sizes
- Theoretical survival curves from transition rates
- Prevalence-weighted overall population curves
- Multiple strata with different treatment effects

## Important design considerations

- **Piecewise rates**: `build_transition_rates()` produces constant rates. To model piecewise PFS or response, modify the output by splitting rows into multiple time periods with different rates and durations.
- **Design vs simulation effects**: Build separate transition rate tables with different HRs for design (sample size) and simulation (operating characteristics). Weaker simulation HRs reveal realistic power under conservative assumptions.
- **Response rate timing**: Multiply the base response rate for early periods (e.g., 2.25× for $t < 6$ months) and reduce it later (e.g., 10% of base) to model response windows.

- **Calibration**: `build_transition_rates()` numerically calibrates control arm rates to match specified median PFS and OS. Experimental arm rates are derived by scaling with HRs, then re-calibrating response rate for the target ORR.
- **`death_wo_prog_rate`**: Baseline rate of death without prior progression (default 0.02/month). This controls the fraction of patients who die without progressing.
- **`responder_prog_ratio`**: Ratio of progression rate for responders vs non-responders (default 0.5). Lower values mean responders have longer PFS.
- **Piecewise rates**: `transition_rate` supports multiple rows per (stratum, treatment, transition) with different durations, enabling time-varying hazards.
- **Enrollment**: Uses `simtrial::rpwexp_enroll()` — specify ramp-up as increasing rates over duration periods.
- **FPE-relative timing**: Analysis cut dates are computed relative to the final patient enrolled (FPE), not study start.
- **Event-driven timing with caps**: `get_analysis_dates()` supports both calendar-based and event-driven triggers with maximum extension caps.
- **Stratified analyses**: For overall population p-values, use `strata()` in `survdiff()` (logrank) or sample-size-weighted stratum-specific risk differences (ORR). Beware that `survdiff()` with `strata()` returns `obs`/`exp` as matrices — sum across strata columns.

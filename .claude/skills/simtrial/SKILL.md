---
name: simtrial
description: >
  Guide users through clinical trial simulation using the simtrial R package.
  Use this skill when the user asks about: simulating survival trials, simfix,
  sim_pw_surv, cutting data at calendar or event times, weighted logrank tests,
  MaxCombo tests, or simulation-based power.
---

# Clinical Trial Simulation with simtrial

## API reference

For full function documentation (arguments, return values, examples), read `references/llms.txt`.
Source: https://gsDesign.ai/packages/simtrial/llms.txt

## Key functions

### Simulation
- `sim_pw_surv()` - Simulate piecewise exponential survival data
- `sim_fixed_n()` - Fixed-sample simulation with analysis pipeline
- `sim_gs_n()` - Group sequential simulation
- `to_sim_pw_surv()` - Convert enrollment/failure rates to sim_pw_surv format

### Data cutting
- `cut_data_by_date()` - Cut simulated data at a calendar date
- `cut_data_by_event()` - Cut simulated data at a target event count
- `get_cut_date_by_event()` - Find calendar date for a target event count
- `get_analysis_date()` - Get analysis date from multiple criteria
- `create_cut()` - Create a cutting function for use in pipelines

### Statistical tests
- `wlr()` - Weighted logrank test
- `maxcombo()` - MaxCombo test (multiple weight functions)
- `rmst()` / `rmst_two_arm()` / `rmst_single_arm()` - Restricted mean survival time
- `milestone()` - Milestone analysis (survival at fixed time)
- `multitest()` - Apply multiple tests to one dataset
- `create_test()` - Create a test function for use in pipelines

### Weight functions
- `fh()` - Fleming-Harrington weights
- `mb()` / `mb_delayed_effect()` - Magirr-Burman weights
- `early_zero()` - Early zero weight function
- `wlr_weight()` - General WLR weight specification

### Utilities
- `counting_process()` - At-risk/event tables from survival data
- `rpwexp()` - Random piecewise exponential generation
- `rpwexp_enroll()` - Random piecewise enrollment times
- `fit_pwexp()` - Fit piecewise exponential model
- `randomize_by_fixed_block()` - Block randomization

### Example datasets
- `ex1_delayed_effect` through `ex6_crossing` - Pre-built scenarios for NPH

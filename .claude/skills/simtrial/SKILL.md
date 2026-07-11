---
name: simtrial
description: >
  Guide users through clinical trial simulation using the simtrial R package.
  Use this skill when the user asks about: simulating survival trials, simfix,
  sim_pw_surv, cutting data at calendar or event times, weighted logrank tests,
  MaxCombo tests, or simulation-based power.
---

# Clinical Trial Simulation with simtrial

**Note**: This skill targets simtrial >= 1.0.2 (main branch at github.com/Merck/simtrial).

## API reference

- Vendored function docs: `references/llms.txt`
- Full function docs (local v1.0.2): `references/llms_local.txt`
- Workflow patterns: `references/code_patterns.md`

## Key functions

### Simulation
- `sim_pw_surv()` - Simulate piecewise exponential survival data (individual patient data)
- `sim_fixed_n()` - Fixed-sample simulation with analysis pipeline (simpler, less flexible)
- `sim_gs_n()` - Group sequential simulation (integrates with gsDesign2 designs)
- `to_sim_pw_surv()` - Convert simple rate format (control rate + HR) to sim_pw_surv format

### Data cutting
- `cut_data_by_date()` - Cut simulated data at a calendar date
- `cut_data_by_event()` - Cut simulated data at a target event count
- `get_cut_date_by_event()` - Find calendar date for a target event count
- `get_analysis_date()` - Get analysis date from multiple criteria (events, calendar time, follow-up)
- `create_cut()` - Create a cutting function for use in sim_gs_n pipelines

### Statistical tests
- `wlr()` - Weighted logrank test (single weight function)
- `maxcombo()` - MaxCombo test (multiple FH weight functions, correlation-adjusted p-value)
- `rmst()` / `rmst_two_arm()` / `rmst_single_arm()` - Restricted mean survival time
- `milestone()` - Milestone analysis (survival difference at fixed time)
- `multitest()` - Apply multiple tests to one dataset
- `create_test()` - Create a parameterized test function for use in pipelines

### Weight functions
- `fh()` - Fleming-Harrington weights (rho, gamma)
- `mb()` - Magirr-Burman weights (delay period for NPH)
- `early_zero()` - Early zero weight function
- `wlr_weight()` - General WLR weight specification

### Utilities
- `counting_process()` - At-risk/event tables from survival data (for custom analyses)
- `rpwexp()` - Random piecewise exponential generation
- `rpwexp_enroll()` - Random piecewise enrollment times
- `fit_pwexp()` - Fit piecewise exponential model
- `randomize_by_fixed_block()` - Block randomization

### Output
- `summary()` - Summarize sim_gs_n results (power, events, timing)
- `as_gt()` - Convert summary to gt table

### Example datasets
- `ex1_delayed_effect` through `ex6_crossing` - Pre-built NPH scenarios from Cross-Pharma Working Group

## Workflow patterns

For detailed code templates, read `references/code_patterns.md`.

Topics covered:
- Generating survival data with `sim_pw_surv()` (PH and NPH)
- Data cutting (calendar, event, flexible criteria)
- Weighted logrank tests (`wlr()`) and MaxCombo (`maxcombo()`)
- Milestone and RMST tests
- Multiple tests with `multitest()` and `create_test()`
- Fixed-sample simulation with `sim_fixed_n()` (timing_type options)
- Group sequential simulation with `sim_gs_n()` (multiple tests, boundary updating)
- Integration with gsDesign2 designs (events from design, `original_design` parameter)
- Weight functions (Fleming-Harrington, Magirr-Burman, early zero)
- Rate format conversion with `to_sim_pw_surv()`
- Flexible analysis timing with `get_analysis_date()` / `create_cut()`
- Stratified simulations
- Standalone `wlr()` with illness-death model ADTTE data
- Example NPH datasets

## Important design considerations

- **`maxcombo` must be used alone** in `sim_gs_n()`: it cannot be combined with other tests in the same `test` list
- **`to_sim_pw_surv()`**: Use this to convert the simpler rate format (control rate + HR) to the treatment-specific format needed by `sim_pw_surv()`
- **`sim_gs_n()` + `original_design`**: Pass a gsDesign2 design object to get boundaries updated based on actual vs planned information fraction
- **`ia_alpha_spending`**: Controls how alpha is spent when observed events differ from planned ("min_planned_actual" is conservative default)
- **`fa_alpha_spending = "full_alpha"`**: Spends full alpha at final analysis (default); use "info_frac" for event underrunning scenarios
- **`create_cut()` and `create_test()`**: These factory functions are essential for building `sim_gs_n()` pipelines
- **Example datasets**: Columns are `id`, `month`, `evntd`, `trt` — need renaming to `tte`, `event`, `treatment` for use with `wlr()` / `maxcombo()`
- **Standalone `wlr()`**: Can be used outside `sim_gs_n()` with any data having `tte`, `event`, `stratum`, `treatment` columns. Returns positive Z when experimental is better. See `references/code_patterns.md` for illness-death model integration.

---
name: gsDesign2
description: >
  Guide users through next-generation group sequential design using the gsDesign2
  R package. Use this skill when the user asks about: gs_design_ahr, gs_power_ahr,
  gs_update_ahr, sequential_pval, average hazard ratio designs, non-proportional
  hazards, piecewise enrollment/failure rates, spending time, or information
  fraction computation.
---

# Group Sequential Design with gsDesign2

**Note**: This skill targets gsDesign2 >= 1.1.9 (main branch at github.com/Merck/gsDesign2).
The `llms.txt` from gsDesign.ai may lag behind; local API docs are in `llms_local.txt`.

## API reference

- Full function docs (CRAN): `references/llms.txt` (source: https://gsDesign.ai)
- Full function docs (local v1.1.9+): `references/llms_local.txt`
- Workflow patterns: `references/code_patterns.md`

## Key functions

### Design (solve for sample size)
- `gs_design_ahr()` - AHR-based group sequential design (primary workhorse)
- `gs_design_wlr()` - Weighted logrank design
- `gs_design_rd()` - Rate difference design (now with minimum risk weighting in 1.1.9)
- `gs_design_combo()` - Combination test design
- `gs_design_npe()` - Non-proportional effect design (general)
- `fixed_design()` - Fixed (non-sequential) designs

### Power (given fixed assumptions)
- `gs_power_ahr()` - Power for AHR designs
- `gs_power_wlr()` - Power for weighted logrank designs
- `gs_power_rd()` - Power for rate difference designs
- `gs_power_combo()` - Power for combination tests
- `gs_power_npe()` - Power for general NPE designs

### Information and bounds
- `gs_info_ahr()` / `gs_info_wlr()` / `gs_info_rd()` / `gs_info_combo()` - Statistical information
- `gs_spending_bound()` - Spending function bounds
- `gs_spending_combo()` - Spending bounds for combination tests
- `gs_b()` - Fixed boundary values
- `gs_bound_summary()` - Formatted bound summary (supports multiple alpha levels)

### Analysis
- `gs_update_ahr()` - Update bounds with observed data (supports stratified piecewise events)
- `gs_cp_npe()` - Conditional power under NPH
- `sequential_pval()` - Sequential p-value for AHR designs (new in 1.1.9)

### Important: `sequential_pval()` vs `sequentialPValue()`
- **`gsDesign2::sequential_pval()`**: Use with gsDesign2 objects (output of `gs_design_ahr()`, `gs_power_ahr()`)
- **`gsDesign::sequentialPValue()`**: Use with gsDesign objects (output of `gsDesign()`, `gsSurv()`, `gsSurvCalendar()`, `gsSurvPower()`)
- Do NOT mix: each function expects the design object from its own package

### Enrollment and failure rates
- `define_enroll_rate()` - Piecewise enrollment rates (supports strata)
- `define_fail_rate()` - Piecewise failure rates with HR (supports strata)

### Utilities
- `ahr()` / `ahr_blinded()` - Average hazard ratio computation
- `expected_accrual()` / `expected_event()` / `expected_time()` - Expected quantities
- `pw_info()` - Piecewise information
- `to_integer()` - Integer sample size rounding
- `ppwe()` / `s2pwe()` - Piecewise exponential utilities
- `wlr_weight()` - Weight functions for weighted logrank

### Output
- `as_gt()` / `as_rtf()` - Table output (footnotes can be suppressed)
- `summary()` / `text_summary()` - Text summaries
- `gs_bound_summary()` - Bound summary table

## Workflow patterns

For detailed code templates, read `references/code_patterns.md`.

Topics covered:
- Enrollment and failure rate setup (piecewise, stratified, scaling to target N)
- Fixed designs and group sequential designs with `gs_design_ahr()`
- Power computation with `gs_power_ahr()` (event = NULL caveat)
- Non-proportional hazards scenarios (delayed effect, crossing, diminishing)
- Spending functions and bound specification
- Spending time (decoupling spending from information fraction)
- Rate difference designs for binary endpoints
- Weighted logrank and MaxCombo combination tests
- Integer rounding with `to_integer()`
- Sequential p-values with `sequential_pval()` for multiplicity
- Updating bounds with `gs_update_ahr()` (stratified piecewise events)
- Conditional power with `gs_cp_npe()`
- Output and reporting (gt, RTF, bound summary)
- `info_scale` options (h0_info, h1_info, h0_h1_info)
- Stratified designs (subgroup + complement populations)

## Important design considerations

- **`info_frac = NULL`**: Use with `analysis_time` to let timing drive the design; gsDesign2 derives the information fraction
- **`event = NULL` in gs_power_ahr**: Always set this when using `analysis_time`; the default `c(30, 40, 50)` causes length mismatches
- **`info_scale = "h0_info"`**: Matches gsDesign convention; recommended for multiplicity workflows
- **Spending time**: Use `timing` in `upar`/`lpar` to decouple spending from information fraction (critical for delayed effects and multiple hypotheses)
- **Stratified designs**: Use `stratum` column in `define_enroll_rate()` and `define_fail_rate()`; scale complement enrollment from subgroup using prevalence
- **`sequential_pval()` (v1.1.9+)**: Works directly with `gs_design_ahr()` output. Use this for gsDesign2-native workflows. For gsDesign objects (from `gsSurv()`, `gsSurvCalendar()`, `gsSurvPower()`), use `gsDesign::sequentialPValue()` instead.
- **`gs_design_ahr()` spending time output (v1.1.8+)**: Can now output spending time directly for transparency
- **`gs_update_ahr()` event_tbl**: Supports piecewise event tables for delayed-effect designs where events per interval are tracked
- **Non-binding futility**: Use `binding = FALSE` so efficacy bounds are computed ignoring the futility bound
- **Spending functions as strings**: `sf = "sfLDOF"` works in addition to `sf = sfLDOF`

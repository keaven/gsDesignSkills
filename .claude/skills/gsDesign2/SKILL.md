---
name: gsDesign2
description: >
  Guide users through next-generation group sequential design using the gsDesign2
  R package. Use this skill when the user asks about: gs_design_ahr, gs_power_ahr,
  gs_update_ahr, sequential_pval, average hazard ratio designs, non-proportional
  hazards, piecewise enrollment/failure rates, or information fraction computation.
---

# Group Sequential Design with gsDesign2

## API reference

For full function documentation (arguments, return values, examples), read `references/llms.txt`.
Source: https://gsDesign.ai/packages/gsdesign2/llms.txt

## Key functions

### Design
- `gs_design_ahr()` - AHR-based group sequential design (primary workhorse)
- `gs_design_wlr()` - Weighted logrank design
- `gs_design_rd()` - Rate difference design
- `gs_design_combo()` - Combination test design
- `gs_design_npe()` - Non-proportional effect design (general)
- `fixed_design()` - Fixed (non-sequential) designs

### Power
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
- `gs_bound_summary()` - Formatted bound summary

### Analysis
- `gs_update_ahr()` - Update bounds with observed data
- `gs_cp_npe()` - Conditional power

### Enrollment and failure rates
- `define_enroll_rate()` - Piecewise enrollment rates (supports strata)
- `define_fail_rate()` - Piecewise failure rates with HR (supports strata)

### Utilities
- `ahr()` / `ahr_blinded()` - Average hazard ratio computation
- `expected_accrual()` / `expected_event()` / `expected_time()` - Expected quantities
- `pw_info()` - Piecewise information
- `to_integer()` - Integer sample size rounding
- `ppwe()` / `s2pwe()` - Piecewise exponential utilities

### Output
- `as_gt()` / `as_rtf()` - Table output
- `summary()` / `text_summary()` - Text summaries

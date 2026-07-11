---
name: multi-endpoint-sim
description: >
  Guide users through multi-endpoint group sequential trial simulation with
  multiplicity-controlled testing. Use this skill when the user asks about:
  simulating trials with OS, PFS, and ORR endpoints, illness-death model simulation
  with gsDesign bounds, sequential p-values in simulation loops, combining
  graphicalMCP with gsDesign for simulation-based operating characteristics,
  cumulative rejection probabilities, or building a full pipeline from design
  through simulation to multiplicity-adjusted testing.
---

# Multi-Endpoint Group Sequential Trial Simulation

This cross-package skill covers the full pipeline for simulating multi-endpoint
trials with group sequential bounds and graphical multiplicity control.

## Required packages

```r
library(gsDesign)       # Group sequential design and sequential p-values
library(gMCPLite)       # Illness-death model (sim_illness_death, cut_illness_death)
library(graphicalMCP)   # Multiplicity graph and graph_test_shortcut
library(simtrial)       # wlr() for logrank Z-statistics
library(parallel)       # mclapply for simulation
```

## When to use this skill

- Designing a trial with multiple endpoints (OS, PFS, ORR) tested at different analyses
- Simulating correlated endpoints via the illness-death model
- Computing simulation-based operating characteristics (power, rejection probabilities)
- Applying graphical multiplicity testing per simulated trial
- Using sequential p-values from `gsDesign::sequentialPValue()` in simulation loops

## Pipeline overview

1. **Design**: Use `gsDesign::gsSurvCalendar()` for the sample-size-driving endpoint (OS),
   `gsDesign::gsSurvPower()` for secondary time-to-event (PFS), and
   `gsDesign::nBinomial()` for binary endpoints (ORR).
2. **Multiplicity graph**: Build with `graphicalMCP::graph_create()` allocating alpha across hypotheses.
3. **Transition rates**: Calibrate with `build_transition_rates()`, modify for piecewise hazards.
4. **Simulation**: Use `sim_illness_death()` + `cut_illness_death()` to generate ADTTE data at each analysis time.
5. **Test statistics**: `simtrial::wlr()` for TTE endpoints, `gsDesign::testBinomial()` for binary.
6. **Sequential testing**: Compute sequential p-values with `gsDesign::sequentialPValue()`, then
   test with `graphicalMCP::graph_test_shortcut()` at each analysis.
7. **Operating characteristics**: Track first rejection analysis per hypothesis, compute cumulative rejection probabilities.

## Key code patterns

For detailed code templates, read `references/code_patterns.md`.

Topics covered:
- Test statistic functions (logrank via wlr, binomial via testBinomial)
- Sign conventions for Z-statistics across packages
- Simulation loop structure (parallel processing)
- Spending time computation from actual vs planned events
- Sequential p-value computation at each analysis
- Per-trial graphical testing loop
- Cumulative rejection probability table
- Correlation matrix of test statistics

## Related skills

- **gsDesign**: Design and spending functions (`gsSurvCalendar`, `sfExtremeValue2`, `testLower`, `sequentialPValue`). Sequential p-value theory: Liu & Anderson (2008).
- **illness-death**: Illness-death model calibration and piecewise rate modification
- **simtrial**: `wlr()` for logrank Z-statistics
- **graphicalMCP**: `graph_create()` and `graph_test_shortcut()` for multiplicity testing
- **graphicalMCP-gsDesign2**: Similar workflow using gsDesign2 instead of gsDesign (observed data, not simulation). Includes theoretical basis for the Maurer-Bretz (2013) framework.

## Important design considerations

- **Theoretical basis**: The per-trial sequential testing loop implements Algorithm 1 of Maurer & Bretz (2013). Sequential p-values from `sequentialPValue()` (Liu & Anderson, 2008) are passed to `graph_test_shortcut()` at each analysis to control FWER.

- **Sign conventions**: `wlr()` returns positive Z when experimental is better. `testBinomial(x1=exp, x2=ctrl)` returns positive Z when experimental is better. Both conventions must match for `sequentialPValue()` (which expects positive Z = favorable).
- **Spending time in simulation**: Use `pmin(planned_events, actual_events) / planned_final_events` at interim analyses and spending time = 1 at final analyses. This prevents over-spending when simulated events exceed planned.
- **ORR is not group sequential**: ORR is tested at a single analysis (e.g., IA2). Use the nominal p-value `pnorm(-z_orr)` directly rather than a sequential p-value.
- **Separate design vs simulation effects**: Use the design HR for sample size and bounds (e.g., 0.75), but weaker simulation HR (e.g., 0.80) for realistic operating characteristics.
- **Graph carries forward**: Use the same `graph_test_shortcut()` call at each analysis with updated sequential p-values. The graph handles alpha reallocation from rejected hypotheses.
- **Track first rejection**: Record the first analysis at which each hypothesis is rejected. Use this to compute cumulative rejection probabilities by analysis.

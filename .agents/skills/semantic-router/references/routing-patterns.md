# Routing Patterns

## First-pass intent fields

Use these fields when converting a user request into semantic intent:

- `endpoint_family`: normal, binary, time_to_event, recurrent_event,
  negative_binomial_count, multi_endpoint, graphical_multiplicity.
- `estimand`: difference, ratio, log_ratio, hazard_ratio, rate_ratio,
  odds_ratio, average_hazard_ratio.
- `hypothesis_type`: superiority, non_inferiority, equivalence,
  super_superiority.
- `design_family`: fixed_sample_size, group_sequential, adaptive,
  sample_size_reestimation, graphical_mcp, simulation.
- `timing_basis`: calendar_time, event_count, information_fraction,
  spending_time, completer_count, accrual_cut.
- `information_source`: planned, observed, blinded, unblinded.
- `boundaries`: efficacy, futility, binding_futility, nonbinding_futility,
  alpha_spending, beta_spending.
- `multiplicity`: none, bonferroni, graphical, gatekeeping, closed_testing,
  parametric.
- `adaptation`: none, blinded_ssr, unblinded_ssr, conditional_power,
  graph_update.
- `simulation`: none, operating_characteristics, type_i_error, power,
  expected_sample_size.

## Package routing cues

Use the smallest matching skill set.

- `gsDesignNB`: recurrent events, negative binomial counts, rate ratios,
  dispersion, event gaps, calendar analyses for count endpoints, blinded or
  unblinded SSR for recurrent events.
- `gsDesign`: classical group sequential boundaries, spending functions,
  survival designs with established gsDesign workflows.
- `gsDesign2`: average hazard ratio, non-proportional hazards, spending time,
  sequential p-values, next-generation survival group sequential designs.
- `simtrial`: survival trial simulation, weighted logrank tests, MaxCombo,
  event-driven data cuts, simulation-based power.
- `graphicalMCP`: multiplicity graphs, graph nodes, transition weights,
  Bonferroni-based closed testing, graphical alpha reallocation.
- `graphicalMCP-gsDesign2`: group sequential designs with graphical MCP,
  Maurer-Bretz style workflows, sequential p-values with multiplicity graphs.
- `rpact`: adaptive designs, inverse normal or Fisher combination tests,
  conditional power, confirmatory adaptive workflows.
- `wpgsd`: correlated test statistics across hypotheses, nested populations,
  weighted parametric multiplicity with group sequential designs.
- `illness-death`: multi-state oncology simulation with correlated OS, PFS,
  and response endpoints.
- `multi-endpoint-sim`: full operating-characteristic simulation across
  multiple endpoints with sequential testing and multiplicity adjustment.

## Clarify when terms are ambiguous

- "Event" can mean a survival event, recurrent event, or count endpoint.
- "Information" can mean planned information, observed information, blinded
  information, or information fraction.
- "Interim" can be calendar based, event-count based, completer-count based,
  or information-fraction based.
- "Multiple endpoints" can require a simple alpha split, a graphical MCP,
  gatekeeping, or a simulation workflow.
- "Adaptive" can mean SSR, conditional power, combination testing, or graph
  updates.

## gsDesignNB workflow matches

Use `../../../crosswalks/gsDesignNB.yaml` for the canonical mapping. Common
matches:

- Fixed recurrent-event sample size: `sample_size_nbinom()`.
- Calendar-time group sequential design: `sample_size_nbinom()`,
  `gsNBCalendar()`, optionally `sim_gs_nbinom()` and `summarize_gs_sim()`.
- Recurrent-event simulation: `nb_sim()` or `nb_sim_seasonal()`.
- Analysis of simulated recurrent-event data: `mutze_test()`.
- Blinded information monitoring: `calculate_blinded_info()`.
- Blinded SSR: `blinded_ssr()` or `sim_ssr_nbinom()` depending on whether the
  user wants one calculation or operating characteristics.
- Unblinded SSR: `unblinded_ssr()` or `sim_ssr_nbinom()`.

## Recommended response shape

For routing-only answers, return:

```yaml
semantic_parse:
  endpoint_family:
  estimand:
  design_family:
  timing_basis:
  information_source:
selected_skills:
  - semantic-router
  - package-skill
recommended_workflow:
  - function_or_step
open_inputs:
  - only inputs that block the workflow
```

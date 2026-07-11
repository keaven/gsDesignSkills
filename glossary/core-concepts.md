# Trial Design Core Concepts

This glossary defines shared semantic fields for routing clinical trial design
requests across the gsDesign ecosystem. Package documentation remains the
source of truth for executable methods and argument details.

## Normalized Intent Record

Use these fields when a request spans multiple concepts or packages:

```yaml
endpoint_family:
estimand:
hypothesis_type:
design_family:
timing_basis:
information_source:
boundaries:
multiplicity:
adaptation:
simulation:
```

Not every request needs every field. Fill only fields that are relevant to the
decision being made.

## Endpoint Family

- `normal`: continuous endpoint analyzed on an approximately normal scale.
- `binary`: response, event/non-event, or proportion endpoint.
- `time_to_event`: survival endpoint such as OS or PFS.
- `recurrent_event`: repeated events per participant over follow-up.
- `negative_binomial_count`: recurrent-event or count endpoint modeled with
  negative binomial overdispersion.
- `multi_endpoint`: multiple clinical endpoints analyzed together.
- `graphical_multiplicity`: a family of hypotheses controlled through a graph.

## Estimand and Effect Scale

- `difference`: treatment minus control on the endpoint scale.
- `ratio`: generic experimental over control ratio.
- `log_ratio`: log transform of a ratio estimand.
- `hazard_ratio`: relative hazard for time-to-event endpoints.
- `average_hazard_ratio`: weighted or average hazard ratio used in
  non-proportional hazards workflows.
- `rate_ratio`: experimental event rate divided by control event rate.
- `odds_ratio`: odds ratio for binary endpoints.

## Hypothesis Type

- `superiority`: show benefit beyond no effect.
- `non_inferiority`: show the treatment is not worse than a margin.
- `equivalence`: show effect is inside an equivalence margin.
- `super_superiority`: show benefit beyond an elevated superiority threshold.

## Design Family

- `fixed_sample_size`: no planned interim decision rule.
- `group_sequential`: interim analyses with efficacy and/or futility
  boundaries.
- `adaptive`: design can change using pre-specified adaptation rules.
- `sample_size_reestimation`: sample size changes after interim information.
- `graphical_mcp`: multiple testing through graph weights and transitions.
- `simulation`: operating characteristics estimated by repeated trial
  simulation.

## Timing and Information

- `calendar_time`: analysis occurs at specified calendar times.
- `event_count`: analysis occurs when a target event count is reached.
- `information_fraction`: analysis uses planned or observed information
  fraction.
- `spending_time`: time scale used by alpha or beta spending.
- `completer_count`: analysis occurs when enough participants complete
  follow-up.
- `accrual_cut`: analysis data are cut by enrollment or calendar rule.
- `planned`: information fixed by the design.
- `observed`: information calculated from observed data.
- `blinded`: information calculated without treatment labels.
- `unblinded`: information calculated using treatment labels.

## Boundaries and Spending

- `efficacy`: upper boundary for evidence of benefit.
- `futility`: lower boundary for stopping due to insufficient benefit.
- `binding_futility`: futility rule included in Type I error control.
- `nonbinding_futility`: futility rule not required for Type I error control.
- `alpha_spending`: controls cumulative Type I error over analyses.
- `beta_spending`: controls cumulative Type II error or futility behavior.

## Multiplicity

- `none`: one primary hypothesis or no family-wise adjustment.
- `bonferroni`: weighted or unweighted Bonferroni allocation.
- `graphical`: graph-based alpha allocation and reallocation.
- `gatekeeping`: ordered or family-based testing.
- `closed_testing`: closed testing principle.
- `parametric`: correlation-based multiplicity adjustment.

## Adaptation

- `none`: no planned adaptation.
- `blinded_ssr`: sample size re-estimation without treatment labels.
- `unblinded_ssr`: sample size re-estimation using treatment labels.
- `conditional_power`: adaptation based on interim conditional power.
- `graph_update`: adaptive change to a multiplicity graph.

## Simulation

- `operating_characteristics`: repeated simulation to estimate trial behavior.
- `power`: probability of rejecting under an alternative.
- `type_i_error`: probability of rejecting under a null.
- `expected_sample_size`: average enrolled or analyzed sample size.
- `early_stopping`: probability and timing of stopping before final analysis.

# Trial Design Semantic Layer Strategy

## Summary

This proposal outlines a semantic layer for clinical trial design across the
`gsDesign` ecosystem of packages, including `gsDesign`, `gsDesign2`,
`simtrial`, `graphicalMCP`, `rpact`, `gsDesignNB`, and future or adjacent
packages such as `gsDesignMB`.

The goal is to make user intent, statistical design concepts, and package
implementations easier to connect. The semantic layer should define shared
trial-design concepts once, then map those concepts into package-specific
functions, arguments, workflows, examples, and validation checks.

`gsDesignSkills` is the right home for this cross-package layer. Individual R
packages should remain the source of truth for executable statistical methods,
but this repository can serve as the hub for vocabulary, routing, adapters,
skills, examples, and evaluation tests.

## Motivation

The current skills already encode substantial design knowledge, but most of it
is package-specific. For example, a user might ask:

> Design a two-arm group sequential trial with a recurrent event endpoint,
> blinded information monitoring, and possible sample size re-estimation.

That request includes concepts that span several semantic levels:

- Endpoint family: recurrent event counts modeled with a negative binomial
  distribution.
- Estimand: treatment effect as a rate ratio.
- Timing: calendar time, follow-up, accrual, event gaps, interim analyses.
- Information: planned, blinded, unblinded, and observed information.
- Design structure: fixed or group sequential design with boundaries.
- Adaptation: blinded or unblinded sample size re-estimation.
- Implementation: `gsDesignNB::sample_size_nbinom()`,
  `gsDesignNB::gsNBCalendar()`, `gsDesignNB::sim_gs_nbinom()`, or
  `gsDesignNB::sim_ssr_nbinom()`.

At present, these mappings are distributed across package documentation,
skills, examples, vignettes, and code patterns. A semantic layer would make the
mapping explicit and reusable.

## Proposed Principle

Each package should be a semantic provider. `gsDesignSkills` should be the
semantic hub.

Package repositories should maintain:

- Statistically validated functions.
- Package-specific documentation.
- Function argument definitions.
- Return classes and object structures.
- Vignettes and examples tied to package behavior.

`gsDesignSkills` should maintain:

- Cross-package vocabulary.
- Package crosswalks.
- Natural-language-to-workflow routing.
- Agent skills and examples.
- Evaluation cases.
- Compatibility notes across packages.

This avoids pushing fast-changing assistant and routing logic into CRAN
packages while still allowing each package to expose stable semantic anchors.

## Proposed Repository Structure

```text
schema/
  trial-design.yaml
  trial-design.schema.json

glossary/
  endpoints.md
  estimands.md
  hypotheses.md
  information.md
  timing.md
  boundaries.md
  multiplicity.md
  adaptation.md
  simulation.md

crosswalks/
  gsDesign.yaml
  gsDesign2.yaml
  simtrial.yaml
  graphicalMCP.yaml
  rpact.yaml
  gsDesignNB.yaml
  gsDesignMB.yaml

.claude/skills/
  semantic-router/
    SKILL.md
    references/
      routing-patterns.md
      concept-map.md

examples/
  negative-binomial-group-sequential.yaml
  survival-group-sequential.yaml
  graphical-mcp-with-gs.yaml
  multi-endpoint-simulation.yaml

evals/
  natural-language-requests.yaml
  expected-routing.yaml
```

This structure separates vocabulary, package mappings, skills, examples, and
test cases while keeping the current package-specific skills intact.

## Core Semantic Concepts

The first version should focus on concepts that already recur across the
existing packages.

### Trial Structure

- Trial phase or setting.
- Number of arms.
- Treatment labels.
- Allocation ratio.
- Accrual model.
- Follow-up model.
- Analysis population.
- Calendar duration.

### Endpoint Family

- Normal endpoint.
- Binary endpoint.
- Time-to-event endpoint.
- Recurrent event or negative binomial endpoint.
- Multiple endpoints.
- Graphical multiple comparison procedure.

### Estimand and Effect Scale

- Difference.
- Ratio.
- Log ratio.
- Hazard ratio.
- Rate ratio.
- Odds ratio.
- Restricted mean or average hazard ratio where applicable.
- Superiority, non-inferiority, equivalence, or super-superiority target.

### Hypotheses and Multiplicity

- Null hypothesis.
- Alternative hypothesis.
- One-sided or two-sided alpha.
- Non-inferiority margin.
- Family-wise error rate.
- Graph nodes, weights, and transitions.
- Gatekeeping or graphical testing strategy.

### Information and Timing

- Planned information.
- Observed information.
- Blinded information.
- Unblinded information.
- Information fraction.
- Spending time.
- Calendar time.
- Event count.
- Completer count.
- Accrual cut date.
- Analysis date.

### Boundaries and Spending

- Upper efficacy boundary.
- Lower futility boundary.
- Binding or non-binding futility.
- Alpha spending function.
- Beta spending function.
- Calendar-based spending.
- Information-based spending.

### Adaptation

- Sample size re-estimation.
- Blinded re-estimation.
- Unblinded re-estimation.
- Conditional power.
- Promising zone.
- Adaptive graph update, if applicable.

### Simulation

- Data-generating model.
- Analysis model.
- Simulation replicate count.
- Random seed.
- Operating characteristics.
- Power.
- Type I error.
- Expected sample size.
- Probability of early stopping.

## Package Crosswalk Pattern

Each package crosswalk should answer the same basic questions:

- What endpoint families does the package support?
- What design families does it support?
- What are the core user-facing functions?
- What high-level concepts map to which arguments?
- What object classes or return structures are produced?
- What concepts are unsupported or only partially supported?
- Which workflows are recommended?
- Which examples are canonical?

The crosswalk files should be simple YAML at first. They can become JSON Schema
or another formal representation later if needed.

Example outline:

```yaml
package: gsDesignNB
role: negative-binomial recurrent event design provider

endpoint_families:
  - recurrent_event
  - negative_binomial_count

design_families:
  - fixed_sample_size
  - group_sequential
  - sample_size_reestimation
  - simulation

concepts:
  control_event_rate:
    maps_to:
      sample_size_nbinom: lambda1
      nb_sim: fail_rate
  experimental_event_rate:
    maps_to:
      sample_size_nbinom: lambda2
      nb_sim: fail_rate
  rate_ratio:
    definition: experimental event rate divided by control event rate
    maps_to:
      sample_size_nbinom: lambda2 / lambda1
      mutze_test: ratio estimate
  dispersion:
    maps_to:
      sample_size_nbinom: dispersion
      nb_sim: dispersion
      calculate_blinded_info: dispersion estimate
  event_gap:
    maps_to:
      sample_size_nbinom: event_gap
      nb_sim: event_gap
  calendar_analysis_times:
    maps_to:
      gsNBCalendar: analysis_times
      sim_gs_nbinom: analysis_times

canonical_workflows:
  fixed_design:
    - sample_size_nbinom
  group_sequential_design:
    - sample_size_nbinom
    - gsNBCalendar
    - sim_gs_nbinom
    - summarize_gs_sim
  sample_size_reestimation:
    - sample_size_nbinom
    - sim_ssr_nbinom
    - summarize_ssr_sim
```

## Initial gsDesignNB Crosswalk

`gsDesignNB` is a good first implementation because it exposes several hard
semantic-layer issues in one package:

- Endpoint-specific terminology.
- Event rates and rate ratios.
- Negative binomial dispersion.
- Minimum event gaps and effective exposure.
- Calendar-time group sequential analyses.
- Blinded versus unblinded information.
- Sample size re-estimation.
- Simulation versus analytical design.

The first `gsDesignNB` crosswalk should cover:

- `sample_size_nbinom()`
- `compute_info_at_time()`
- `gsNBCalendar()`
- `nb_sim()`
- `nb_sim_seasonal()`
- `sim_gs_nbinom()`
- `sim_ssr_nbinom()`
- `blinded_ssr()`
- `unblinded_ssr()`
- `calculate_blinded_info()`
- `mutze_test()`
- `summarize_gs_sim()`
- `summarize_ssr_sim()`

It should also document semantic anchors for return objects such as
`sample_size_nbinom_result` and `gsNB`.

## Relationship to Existing Skills

The existing package skills should continue to work as they do now. The
semantic layer should add a routing layer above them.

For example:

1. The semantic router parses a request into endpoint, estimand, design family,
   timing, multiplicity, and simulation requirements.
2. The router selects one or more package skills.
3. The selected package skill provides package-specific code patterns and
   implementation details.
4. The semantic layer validates that the selected workflow matches the user
   intent.

This keeps each package skill focused while allowing cross-package designs to
be composed more reliably.

## Non-Goals

The first version should not attempt to:

- Replace package documentation.
- Reimplement statistical calculations.
- Force all packages into one common R API.
- Encode every possible trial design method.
- Add CRAN-facing dependencies or assistant logic to R packages.
- Build a full ontology before proving usefulness with examples.

The first useful version should be small, inspectable, and practical.

## Implementation Phases

### Phase 1: Minimum Viable Semantic Layer

- Add `glossary/` pages for core trial design concepts.
- Add `crosswalks/gsDesignNB.yaml`.
- Add a `semantic-router` skill that uses the glossary and crosswalk.
- Add 10 to 20 natural-language examples with expected routing.
- Validate the examples manually against existing package skills.

### Phase 2: Cross-Package Coverage

- Add crosswalks for `gsDesign`, `gsDesign2`, `simtrial`, `graphicalMCP`, and
  `rpact`.
- Add examples that require multiple packages.
- Add routing rules for survival designs, graphical MCP, and simulation.
- Identify where terminology differs across packages.

### Phase 3: Evaluation and Regression Tests

- Add machine-readable evaluation cases.
- Track expected package, function, and workflow selection.
- Include ambiguous prompts and expected clarification questions.
- Add examples for unsupported requests and safe refusal or redirection.

### Phase 4: Documentation Site Integration

- Publish the glossary and crosswalks through the Quarto site.
- Add a semantic-layer overview page.
- Link package skill pages to their crosswalks.
- Add a maintainer guide for adding new package adapters.

## Acceptance Criteria for Phase 1

Phase 1 is complete when:

- A reviewer can inspect the glossary and understand the core vocabulary.
- `gsDesignNB` has a package crosswalk covering its main design and simulation
  workflows.
- The semantic-router skill can route negative binomial design requests to the
  correct `gsDesignNB` workflows.
- At least 10 natural-language examples have expected routing outputs.
- The proposal does not require changes to the `gsDesignNB` CRAN package.

## Open Questions

- Should crosswalk files be YAML only, or should they be validated by JSON
  Schema from the start?
- Should `gsDesignSkills` eventually expose a small website page for each
  package crosswalk?
- Should package repositories include a minimal semantic manifest, or should
  all crosswalks live only in `gsDesignSkills`?
- How should versioning work when package APIs change?
- Should `gsDesignMB` be added immediately as a placeholder, or only when its
  scope is more clearly defined?

## Proposed GitHub Issue

Title:

```text
Create a cross-package trial design semantic layer
```

Body:

```markdown
## Goal

Create a semantic layer in `gsDesignSkills` that maps high-level clinical trial
design concepts to package-specific workflows across `gsDesign`, `gsDesign2`,
`simtrial`, `graphicalMCP`, `rpact`, `gsDesignNB`, and future adjacent
packages.

## Rationale

The existing skills encode useful package-specific knowledge, but there is no
shared vocabulary or crosswalk for concepts such as endpoint family, estimand,
information, timing, boundaries, multiplicity, adaptation, and simulation. A
semantic layer would make package routing and multi-package workflows more
reliable.

## Initial Scope

- Add a small glossary of core trial-design concepts.
- Add a package crosswalk format.
- Implement the first crosswalk for `gsDesignNB`.
- Add a semantic-router skill.
- Add natural-language routing examples and expected outputs.

## Deliverables

- `glossary/*.md`
- `crosswalks/gsDesignNB.yaml`
- `.claude/skills/semantic-router/SKILL.md`
- `.claude/skills/semantic-router/references/routing-patterns.md`
- `evals/natural-language-requests.yaml`
- Links from the existing skills index or README as appropriate.

## Acceptance Criteria

- `gsDesignNB` negative binomial design requests can be routed to the correct
  package workflows.
- The design is inspectable and does not require changes to CRAN packages.
- The structure can be extended to `gsDesign`, `gsDesign2`, `simtrial`,
  `graphicalMCP`, and `rpact`.
```

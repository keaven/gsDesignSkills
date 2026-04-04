---
name: wpgsd
description: >
  Guide users through weighted parametric group sequential design using the
  wpgsd R package. Use this skill when the user asks about: correlated test
  statistics across hypotheses, generate_bounds, closed_test, correlation
  matrices for nested populations, or parametric multiplicity adjustment
  with group sequential designs.
---

# Weighted Parametric Group Sequential Design with wpgsd

**Note**: This skill targets wpgsd from github.com/Merck/wpgsd.

## API reference

- Full function docs: `references/llms.txt` (source: local Rd man pages)
- Workflow patterns: `references/code_patterns.md`

## Key functions

### Bounds and testing
- `generate_bounds()` - Compute group sequential bounds accounting for correlations across hypotheses
- `closed_test()` - Closed testing procedure using weighted parametric tests
- `calc_seq_p()` - Calculate sequential p-values

### Correlation structure
- `generate_corr()` - Generate correlation matrix from event counts (for nested populations)
- `generate_event_table()` - Build event count table for correlation computation

### Internal utilities
- `find_astar()` - Find adjusted alpha for spending function
- `find_xi()` - Find xi parameter for bounds

## Workflow patterns

For detailed code templates, read `references/code_patterns.md`.

Topics covered:
- Event count tables and correlation matrix generation
- Generating event tables from ADSL/ADTTE datasets
- Bonferroni bounds (type = 0, baseline comparison)
- WPGSD bounds with overall alpha spending (type = 2, Method 3b)
- WPGSD bounds with separate alpha spending (type = 3, Method 3c)
- WPGSD bounds with fixed alpha spending (type = 1, Method 3a)
- Closed testing with observed p-values
- Sequential p-values with `calc_seq_p()` (WPGSD and Bonferroni)
- Overlapping populations example (biomarker subgroups + overall)
- Common control example (multiple experimental arms vs control)
- Manual correlation construction for k > 2 analyses

## Important design considerations

- **`generate_bounds()` type parameter**: 0 = Bonferroni (baseline), 1 = fixed spending, 2 = overall spending (single SF), 3 = separate spending (per-hypothesis SFs with inflation factor xi)
- **Spending functions**: For type 0 and 3, `sf`/`sfparm`/`t` are lists (one per hypothesis); for type 2, they are scalars (single overall SF)
- **Correlation sources**: Overlapping populations (shared patients), nested populations (subgroup ⊂ overall), or common control arm (shared control events)
- **Event table (H1, H2) pairs**: `(i, i)` = events for hypothesis i alone; `(i, j)` = events in the intersection of populations i and j
- **Nested populations**: If H1 ⊂ H3 (e.g., subgroup within overall), then `n_{1∧3} = n_1` (all H1 events are also H3 events)
- **Common control**: If H1, H2 share a control arm, then `n_{1∧2}` = control arm events only
- **Closed testing is required**: WPGSD does not guarantee consonance, so full closed testing via `closed_test()` is needed (not just shortcut testing)
- **Power at design stage**: Either (a) design with Bonferroni and use WPGSD only at analysis, or (b) use the minimum p-value bound across intersection hypotheses to power each elementary hypothesis

## Known limitations

- **`generate_corr()` bug for k > 2**: Incorrectly computes within-hypothesis cross-analysis entries for non-adjacent analyses. For k > 2, build the event-count matrix D manually and compute `corr = diag(1/sqrt(diag(D))) %*% D %*% diag(1/sqrt(diag(D)))`. See `code_patterns.md` section on manual correlation.

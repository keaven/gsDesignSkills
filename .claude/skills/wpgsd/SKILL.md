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

## API reference

For full function documentation (arguments, return values, examples), read `references/llms.txt`.
Source: https://merck.github.io/wpgsd/

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

## Known limitations

- `generate_corr()` has a bug for k > 2 analyses: it incorrectly computes within-hypothesis cross-analysis entries for non-adjacent analyses. For k > 2, build the event-count matrix D manually and compute `corr = diag(1/sqrt(diag(D))) %*% D %*% diag(1/sqrt(diag(D)))`.

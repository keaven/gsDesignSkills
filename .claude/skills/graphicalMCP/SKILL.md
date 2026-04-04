---
name: graphicalMCP
description: >
  Guide users through graphical multiple comparison procedures using the
  graphicalMCP R package. Use this skill when the user asks about: multiplicity
  graphs, Bonferroni-based procedures, graph_create, graph_test_shortcut,
  graph_update, transition matrices, alpha reallocation, or closed testing
  with graphs.
---

# Graphical Multiple Comparison Procedures with graphicalMCP

## API reference

- Full function docs: `references/llms.txt` (source: https://openpharma.github.io/graphicalMCP/)
- Workflow patterns: `references/code_patterns.md`

## Key functions

### Graph creation and manipulation
- `graph_create()` - Create a multiplicity graph (hypotheses, weights, transitions)
- `graph_update()` - Update graph by deleting rejected hypotheses
- `as_graph()` - Convert from gMCP, igraph, or matrix objects

### Testing
- `graph_test_shortcut()` - Shortcut (Bonferroni-based) graphical testing
- `graph_test_closure()` - Full closure-based testing (supports Simes, parametric)
- `graph_rejection_orderings()` - Enumerate all valid rejection orderings

### Adjusted p-values and weights
- `adjust_p()` - Adjusted p-values (Bonferroni, Simes, parametric)
- `adjust_weights()` - Adjusted significance levels (Bonferroni, Simes, parametric)
- `graph_generate_weights()` - Generate weights for all intersection hypotheses

### Power
- `graph_calculate_power()` - Power simulation via multivariate normal

### Plotting
- `plot.initial_graph()` - Plot multiplicity graph
- `plot.updated_graph()` - Plot updated graph sequence

### Example graphs
- `example_graphs()` - Pre-built example graphs (Bonferroni-Holm, fixed sequence, etc.)

## Workflow patterns

For detailed code templates, read `references/code_patterns.md`.

Topics covered:
- Creating graphs with `graph_create()` (weights, transitions)
- Common graph structures (fixed sequence, Bonferroni-Holm, fallback, successive)
- Shortcut (Bonferroni) testing with `graph_test_shortcut()`
- Closure testing with `graph_test_closure()` (Simes, parametric, Hochberg, mixed)
- Updating graphs manually with `graph_update()`
- Generating intersection weights with `graph_generate_weights()`
- Rejection orderings with `graph_rejection_orderings()`
- Power simulation with `graph_calculate_power()` (marginal power, correlation)
- Custom success criteria for power (e.g., "reject primary + secondary")
- Plotting graphs (layout, edge curves, epsilon edges, colors)
- Built-in example graphs (20+ pre-built procedures)
- Converting between graphicalMCP, gMCPLite, and igraph formats
- Parametric tests with known correlation (nested populations, shared controls)
- Multi-population multi-endpoint designs (oncology PFS + OS in subgroup + overall)

## Important design considerations

- **`graph_test_shortcut()` vs `graph_test_closure()`**: Shortcut is Bonferroni-only and faster; closure supports Simes, parametric, and Hochberg tests for more power
- **Simes test**: Valid when test statistics are positively correlated (PRDS condition); more powerful than Bonferroni for correlated hypotheses
- **Parametric test**: Requires known `test_corr` between test statistics; most powerful when correlation is high
- **`test_groups` and `test_types`**: Allow different test types for different groups of hypotheses (e.g., parametric for co-primary, Simes for secondary)
- **`test_corr`**: Must match `test_types` in length; use `NA` for non-parametric tests (Bonferroni, Simes, Hochberg)
- **`sim_corr` vs `test_corr`**: `sim_corr` is for generating correlated p-values in power simulation; `test_corr` is the known correlation used in the parametric test itself
- **`power_marginal`**: Marginal power for each hypothesis at full alpha (not at allocated alpha); higher values = stronger signal
- **Transition matrix rows must sum to ≤ 1**: Each row represents how alpha propagates when that hypothesis is rejected
- **Sequential p-values**: Use `sequential_pval()` from gsDesign2 to get p-values for group sequential hypotheses, then pass to `graph_test_shortcut()` for multiplicity control

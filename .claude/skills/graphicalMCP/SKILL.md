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

For full function documentation (arguments, return values, examples), read `references/llms.txt`.
Source: https://openpharma.github.io/graphicalMCP/

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

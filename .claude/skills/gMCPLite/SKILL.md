---
name: gMCPLite
description: >
  Guide users through graphical MCP procedures using the gMCPLite R package
  (legacy). Use this skill when the user asks about: hGraph for multiplicity
  graph visualization, gMCP for closed testing, or legacy graphical MCP
  workflows. For new projects, prefer graphicalMCP.
---

# Graphical MCP with gMCPLite (Legacy)

For new projects, prefer the `graphicalMCP` package which has a cleaner API.

## API reference

For full function documentation (arguments, return values, examples), read `references/llms.txt`.
Built from local man pages.

## Key functions

### Graph creation
- `matrix2graph()` - Create graphMCP object from transition matrix
- `graphMCP` class - Core graph representation (hypotheses, weights, transitions)
- `joinGraphs()` - Combine multiple graphs
- `subgraph()` - Extract subgraph

### Testing
- `gMCP()` - Graphical MCP testing procedure
- `gMCP.extended()` - Extended testing with parametric tests
- `graphTest()` - Test hypotheses on a graph

### Visualization
- `hGraph()` - Create multiplicity graph visualization (ggplot2-based)
- `placeNodes()` - Compute node positions for graph layout

### Test functions
- `bonferroni.test()` - Bonferroni test
- `bonferroni.trimmed.simes.test()` - Bonferroni-trimmed Simes test
- `parametric.test()` - Parametric test using correlation
- `simes.test()` - Simes test
- `simes.on.subsets.test()` - Simes test on subsets

### Utilities
- `generateWeights()` - Generate weights for intersection hypotheses
- `generatePvals()` - Generate p-values for simulation
- `simConfint()` - Simultaneous confidence intervals
- `rejectNode()` - Reject a hypothesis and update graph
- `exampleGraphs()` - Pre-built example graphs
- `checkCorrelation()` - Validate correlation matrix

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

- Full function docs: `references/llms.txt` (built from local man pages)
- Workflow patterns: `references/code_patterns.md`

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

## Workflow patterns

For detailed code templates, read `references/code_patterns.md`.

Topics covered:
- Creating multiplicity graphs with `hGraph()` (basic and custom)
- Customizing hGraph layout (positions, colors, legends, sizing, radianStart)
- Creating graphMCP objects with `matrix2graph()`
- Bonferroni, Simes, and parametric testing with `gMCP()`
- Extended testing with `gMCP.extended()` and custom test functions
- Generating intersection weights with `generateWeights()`
- Updating graphs after rejection with `rejectNode()`
- Simultaneous confidence intervals with `simConfint()`
- Built-in example graphs (BonferroniHolm, fixedSequence, fallback, etc.)
- Integration with gsDesign sequential p-values (`sequentialPValue()`)
- Complex oncology trial template (6 hypotheses: OS/PFS/ORR x Subgroup/All)
- Combining and subsetting graphs (`joinGraphs()`, `subgraph()`)

## Important design considerations

- **For new projects, prefer `graphicalMCP`**: It has a cleaner S3 API (`graph_create`, `graph_test_shortcut`, `graph_test_closure`) and is actively maintained
- **`hGraph()` remains widely used**: Even with graphicalMCP for testing, `hGraph()` from gMCPLite is commonly used for visualization in publications and presentations
- **Sequential p-values workflow**: Use `gsDesign::sequentialPValue()` to convert nominal p-values from group sequential analyses into sequential p-values, then pass to `gMCP()` for multiplicity control
- **`upscale = TRUE`**: Required for parametric tests (Bretz et al. 2011) to rescale subgraph weights to sum to 1
- **`correlation` with NA**: gMCPLite supports partially specified correlation matrices (NA for unknown entries)
- **Time travel for alpha**: When a hypothesis is rejected at a later analysis, previously tested hypotheses can be re-tested at updated alpha levels — this controls Type I error but requires careful bound re-derivation
- **`gMCP()` returns `gMCPResult`**: Access `@rejected` (logical), `@adjPValues` (adjusted p-values), and `@graphs` (sequence of updated graphs)

# Code Patterns for gMCPLite

**Note**: For new projects, prefer `graphicalMCP` which has a cleaner API.
gMCPLite is the legacy package and remains useful for `hGraph()` visualization
and integration with existing gsDesign workflows.

## Table of Contents
1. [Creating multiplicity graphs with hGraph](#hgraph)
2. [Customizing hGraph layout and appearance](#hgraph-custom)
3. [Creating graphMCP objects](#graph-create)
4. [Testing with gMCP](#gmcp-test)
5. [Extended testing with gMCP.extended](#gmcp-extended)
6. [Generating intersection weights](#generate-weights)
7. [Updating graphs (rejecting hypotheses)](#reject-node)
8. [Simultaneous confidence intervals](#sim-confint)
9. [Example graphs](#example-graphs)
10. [Integration with gsDesign (sequential p-values)](#gsdesign-integration)
11. [Complex oncology trial template](#oncology-template)
12. [Combining and subsetting graphs](#join-subgraph)

---

## Creating multiplicity graphs with hGraph {#hgraph}

`hGraph()` creates ggplot2-based multiplicity graph visualizations. This is
the primary function most users need from gMCPLite.

### Basic graph (equal weights)

```r
library(gMCPLite)

# Default: 4 hypotheses, equal alpha split, equal transitions
hGraph()
```

### Custom hypotheses and weights

```r
hGraph(
  nHypotheses = 3,
  nameHypotheses = c("H1: OS", "H2: PFS", "H3: ORR"),
  alphaHypotheses = c(0.01, 0.01, 0.005),
  m = matrix(c(
    0, 1, 0,
    0, 0, 1,
    0.5, 0.5, 0
  ), nrow = 3, byrow = TRUE)
)
```

### Using weights that sum to 1 (instead of alpha)

```r
hGraph(
  nHypotheses = 3,
  nameHypotheses = c("H1", "H2", "H3"),
  alphaHypotheses = c(0.5, 0.3, 0.2),
  wchar = "w"   # Display as w1, w2, w3 instead of alpha
)
```

### Two-hypothesis graph

```r
hGraph(
  nHypotheses = 2,
  nameHypotheses = c("H1: Primary", "H2: Secondary"),
  alphaHypotheses = c(0.02, 0.005),
  m = matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
)
```

---

## Customizing hGraph layout and appearance {#hgraph-custom}

### Custom node positions (x, y)

```r
# 6-hypothesis oncology trial layout (3 rows x 2 columns)
hGraph(
  nHypotheses = 6,
  nameHypotheses = c(
    "H1: OS\nSubgroup", "H2: OS\nAll",
    "H3: PFS\nSubgroup", "H4: PFS\nAll",
    "H5: ORR\nSubgroup", "H6: ORR\nAll"
  ),
  alphaHypotheses = c(0.01, 0.01, 0.004, 0, 0.0005, 0.0005),
  m = matrix(c(
    0, 1, 0,   0,   0,   0,
    0, 0, 0.5, 0.5, 0,   0,
    0, 0, 0,   1,   0,   0,
    0, 0, 0,   0,   0.5, 0.5,
    0, 0, 0,   0,   0,   1,
    0.5, 0.5, 0, 0, 0,   0
  ), nrow = 6, byrow = TRUE),
  x = c(-1.25, 1.25, -2.5, 2.5, -1.25, 1.25),
  y = c(2, 2, 1, 1, 0, 0),
  halfWid = 1,
  halfHgt = 0.35,
  trhw = 0.15,
  trprop = 0.4,
  offset = 0
)
```

### Controlling radianStart (circular layout)

```r
# Place first hypothesis at top center
hGraph(
  nHypotheses = 3,
  radianStart = pi / 2,
  offset = pi / 20
)

# Default formula: pi * (1/2 + 1/n) for odd n, pi * (1 + 2/n)/2 for even n
```

### Ellipse and text sizing

```r
hGraph(
  nHypotheses = 4,
  size = 5,           # Hypothesis text size (default 6)
  boxtextsize = 3,    # Transition weight text size (default 4)
  halfWid = 1,        # Ellipse half-width (default 0.5)
  halfHgt = 0.35,     # Ellipse half-height (default 0.5)
  trhw = 0.15,        # Transition box half-width (default 0.1)
  trhh = 0.1,         # Transition box half-height (default 0.075)
  digits = 4,         # Digits for alpha display (default 5)
  trdigits = 3,       # Digits for transition weights (default 2)
  arrowsize = 0.03,   # Arrow head size (default 0.02)
  xradius = 2.5,      # X-radius of circular layout (default 2)
  yradius = 1.5       # Y-radius of circular layout (default = xradius)
)
```

### Colors and legends

```r
# Color by endpoint group
cbPalette <- c("#999999", "#E69F00", "#56B4E9")

hGraph(
  nHypotheses = 6,
  nameHypotheses = c("H1: OS Sub", "H2: OS All",
                     "H3: PFS Sub", "H4: PFS All",
                     "H5: ORR Sub", "H6: ORR All"),
  alphaHypotheses = c(0.01, 0.01, 0.004, 0, 0.0005, 0.0005),
  fill = as.character(c(1, 1, 2, 2, 3, 3)),   # Group by endpoint
  palette = cbPalette,
  labels = c("OS", "PFS", "ORR"),              # Legend labels
  legend.name = "Endpoint",
  legend.position = "right"
)
```

### Transition line placement

```r
hGraph(
  nHypotheses = 4,
  trprop = 0.15,   # Slide boxes closer to source (default 1/3)
  offset = 0       # No offset for single-direction transitions
)
```

---

## Creating graphMCP objects {#graph-create}

`graphMCP` objects are needed for `gMCP()` testing. Create from a
transition matrix with `matrix2graph()`.

### From transition matrix

```r
m <- matrix(c(
  0, 1, 0,
  0, 0, 1,
  0.5, 0.5, 0
), nrow = 3, byrow = TRUE)

w <- c(0.5, 0.3, 0.2)

graph <- matrix2graph(m, weights = w)
graph
```

### Accessing graph properties

```r
graph@m         # Transition matrix
graph@weights   # Hypothesis weights

# Or use accessor methods
getMatrix(graph)
getWeights(graph)
```

### Converting back to matrix

```r
m_back <- graph2matrix(graph)
```

---

## Testing with gMCP {#gmcp-test}

`gMCP()` performs the Maurer-Bretz graphical testing procedure.

### Bonferroni-based testing (default)

```r
m <- matrix(c(
  0, 0.5, 0.5,
  0.5, 0, 0.5,
  0.5, 0.5, 0
), nrow = 3, byrow = TRUE)

graph <- matrix2graph(m, weights = c(1/3, 1/3, 1/3))

pvalues <- c(0.01, 0.007, 0.1)

result <- gMCP(graph, pvalues = pvalues, alpha = 0.025)
result

# Access results
result@rejected       # Logical vector of rejections
result@adjPValues     # Adjusted p-values
result@graphs         # Sequence of updated graphs
```

### With verbose output

```r
result <- gMCP(graph, pvalues = pvalues, alpha = 0.025, verbose = TRUE)
# Shows step-by-step rejection sequence
```

### Testing with Simes test

```r
result <- gMCP(graph, pvalues = pvalues, alpha = 0.025, test = "Simes")
```

### Testing with parametric test (known correlation)

```r
corr <- matrix(c(1, 0.5, 0.3,
                 0.5, 1, 0.4,
                 0.3, 0.4, 1), nrow = 3)

result <- gMCP(graph, pvalues = pvalues, alpha = 0.025,
               test = "parametric", correlation = corr)
```

### Using upscale

```r
# upscale = TRUE: rescale weights in subgraphs to sum to 1
# Required for parametric tests (Bretz et al. 2011)
result <- gMCP(graph, pvalues = pvalues, alpha = 0.025,
               test = "parametric", correlation = corr, upscale = TRUE)
```

---

## Extended testing with gMCP.extended {#gmcp-extended}

`gMCP.extended()` allows specifying custom test functions.

```r
# Using Simes test
result <- gMCP.extended(
  graph,
  pvalues = pvalues,
  test = simes.test,
  alpha = 0.025
)
result@rejected

# Using parametric test with correlation
result <- gMCP.extended(
  graph,
  pvalues = pvalues,
  test = parametric.test,
  alpha = 0.025,
  correlation = corr
)

# Using bonferroni.trimmed.simes.test
result <- gMCP.extended(
  graph,
  pvalues = pvalues,
  test = bonferroni.trimmed.simes.test,
  alpha = 0.025
)
```

---

## Generating intersection weights {#generate-weights}

`generateWeights()` computes the weighting strategy for all 2^m - 1
intersection hypotheses.

```r
m <- matrix(c(
  0, 0.5, 0.5,
  0.5, 0, 0.5,
  0.5, 0.5, 0
), nrow = 3, byrow = TRUE)

w <- c(1/3, 1/3, 1/3)

weights <- generateWeights(m, w)
weights
# Matrix: first m columns = membership (1/0), last m columns = allocated weights
# Row per intersection hypothesis (2^m - 1 = 7 rows for 3 hypotheses)
```

---

## Updating graphs (rejecting hypotheses) {#reject-node}

`rejectNode()` updates a graph after rejecting a specific hypothesis,
reallocating its alpha according to transition weights.

```r
graph <- matrix2graph(m, weights = c(0.5, 0.3, 0.2))

# Reject H1
updated <- rejectNode(graph, "H1")
getWeights(updated)  # H1 weight redistributed to H2, H3

# Reject H1 then H2
updated2 <- rejectNode(updated, "H2")
getWeights(updated2)
```

---

## Simultaneous confidence intervals {#sim-confint}

`simConfint()` computes simultaneous confidence intervals consistent with
the graphical testing procedure.

```r
graph <- matrix2graph(m, weights = c(0.5, 0.3, 0.2))
pvalues <- c(0.01, 0.02, 0.05)

# Normal-based simultaneous CIs
ci <- simConfint(
  object = graph,
  pvalues = pvalues,
  confint = "normal",
  alternative = "greater",
  estimates = c(0.5, 0.3, 0.1),   # Point estimates
  alpha = 0.025
)
ci  # Matrix: lower bound, estimate, upper bound

# t-based CIs (specify degrees of freedom)
ci_t <- simConfint(
  object = graph,
  pvalues = pvalues,
  confint = "t",
  alternative = "greater",
  estimates = c(0.5, 0.3, 0.1),
  df = c(100, 100, 100),
  alpha = 0.025
)
```

---

## Example graphs {#example-graphs}

gMCPLite provides many pre-built example graphs as `graphMCP` objects.

```r
# Bonferroni-Holm
g <- BonferroniHolm(3)

# Fixed sequence
g <- fixedSequence(3)

# Fallback
g <- fallback(weights = c(0.5, 0.3, 0.2))

# Bretz et al. 2011 (two primary + two secondary)
g <- BretzEtAl2011()

# Cycle graph
g <- cycleGraph(nodes = c("H1", "H2", "H3"), weights = c(1/3, 1/3, 1/3))

# Parallel gatekeeping
g <- parallelGatekeeping()
g <- improvedParallelGatekeeping()

# Truncated Holm
g <- truncatedHolm(gamma = 0.5)

# Other examples
g <- BauerEtAl2001()
g <- HommelEtAl2007()
g <- MaurerEtAl1995()
g <- HungEtWang2010()
g <- improvedFallbackI()
g <- improvedFallbackII()
```

---

## Integration with gsDesign (sequential p-values) {#gsdesign-integration}

The main use case for gMCPLite with group sequential designs: compute
sequential p-values from gsDesign, then test with gMCP.

### Workflow overview

1. Design each hypothesis with `gsDesign::gsSurv()` at its allocated alpha
2. At analysis, compute nominal p-values for each hypothesis
3. Compute sequential p-values using `gsDesign::sequentialPValue()`
4. Pass sequential p-values to `gMCP()` for multiplicity testing
5. If hypotheses are rejected, update alpha and re-test with updated bounds

### Computing sequential p-values

```r
library(gsDesign)
library(gMCPLite)

# Design for a hypothesis (e.g., OS in subgroup)
design <- gsSurv(
  k = 3, test.type = 1, alpha = 0.01, beta = 0.1,
  hr = 0.65, timing = c(0.6, 0.8),
  sfu = sfLDOF,
  lambdaC = log(2) / 12,
  eta = 0.001,
  gamma = c(2.5, 5, 7.5, 10),
  R = c(2, 2, 2, 12),
  T = 42, minfup = 24
)

# At analysis: observed events and nominal p-values
# (from logrank test or similar)
observed_events <- c(150, 200, 260)
nominal_p <- c(0.015, 0.005, 0.001)

# Compute sequential p-value at analysis k
# Sequential p-value = minimum alpha at which the boundary would be crossed
seq_p <- sequentialPValue(
  gsD = design,
  n.I = observed_events,
  Z = -qnorm(nominal_p),
  usTime = observed_events / max(design$n.I)  # Spending time
)
```

### Full testing workflow

```r
# Multiplicity graph
m <- matrix(c(
  0, 1,
  1, 0
), nrow = 2, byrow = TRUE)
graph <- matrix2graph(m, weights = c(0.5, 0.5))

# Sequential p-values for each hypothesis
seq_p <- c(seq_p_H1, seq_p_H2)  # From sequentialPValue() calls

# Test with gMCP
result <- gMCP(graph, pvalues = seq_p, alpha = 0.025)
result@rejected

# If H1 rejected, update design for H2 with increased alpha
if (result@rejected[1]) {
  # H2 now tested at full alpha (0.025)
  # Re-derive bounds for H2 at alpha = 0.025
  design_H2_updated <- gsDesign(
    k = length(observed_events_H2),
    test.type = 1,
    alpha = 0.025,
    sfu = sfLDOF,
    n.I = observed_events_H2,
    maxn.IPlan = max(design_H2$n.I)
  )
}
```

---

## Complex oncology trial template {#oncology-template}

Full template for a multi-endpoint, multi-population oncology trial
(6 hypotheses: OS/PFS/ORR x Subgroup/All).

### Define the multiplicity graph

```r
library(gsDesign)
library(gMCPLite)

nameHypotheses <- c(
  "H1: OS\nSubgroup", "H2: OS\nAll",
  "H3: PFS\nSubgroup", "H4: PFS\nAll",
  "H5: ORR\nSubgroup", "H6: ORR\nAll"
)

alphaHypotheses <- c(0.01, 0.01, 0.004, 0, 0.0005, 0.0005)

m <- matrix(c(
  0, 1, 0,   0,   0,   0,
  0, 0, 0.5, 0.5, 0,   0,
  0, 0, 0,   1,   0,   0,
  0, 0, 0,   0,   0.5, 0.5,
  0, 0, 0,   0,   0,   1,
  0.5, 0.5, 0, 0, 0,   0
), nrow = 6, byrow = TRUE)

# Visualize
cbPalette <- c("#999999", "#E69F00", "#56B4E9")
g <- hGraph(6,
  alphaHypotheses = alphaHypotheses, m = m,
  nameHypotheses = nameHypotheses,
  palette = cbPalette,
  halfWid = 1, halfHgt = 0.35, xradius = 2.5, yradius = 1,
  offset = 0, trhw = 0.15,
  x = c(-1.25, 1.25, -2.5, 2.5, -1.25, 1.25),
  y = c(2, 2, 1, 1, 0, 0),
  trprop = 0.4,
  fill = as.character(c(2, 2, 4, 4, 3, 3))
)
print(g)
```

### Design each hypothesis

```r
# OS Subgroup (H1): 3-analysis GSD
ossub <- gsSurv(
  k = 3, test.type = 1, alpha = alphaHypotheses[1], beta = 0.1,
  hr = 0.65, timing = c(0.61, 0.82), sfu = sfLDOF,
  lambdaC = log(2)/12, eta = 0.001,
  gamma = c(2.5, 5, 7.5, 10), R = c(2, 2, 2, 12),
  T = 42, minfup = 24
)

# PFS Subgroup (H3): 2-analysis GSD with futility
pfssub <- gsSurv(
  k = 2, test.type = 4, alpha = alphaHypotheses[3], beta = 0.1,
  hr = 0.65, timing = 0.7, sfu = sfLDOF, sfl = sfHSD, sflpar = -2,
  lambdaC = log(2)/5, eta = 0.001,
  gamma = c(2.5, 5, 7.5, 10), R = c(2, 2, 2, 12),
  T = 32, minfup = 14
)

# ORR Subgroup (H5): single analysis (fixed design)
```

### Test at analysis time

```r
# Sequential p-values for each hypothesis
# (computed from nominal p-values and gsDesign objects)

graph_mcp <- matrix2graph(m, weights = alphaHypotheses / sum(alphaHypotheses))

# Test
result <- gMCP(graph_mcp, pvalues = seq_p_all, alpha = sum(alphaHypotheses))
result@rejected

# Visualize result sequence
for (i in seq_along(result@graphs)) {
  g_i <- result@graphs[[i]]
  print(hGraph(
    nHypotheses = 6,
    alphaHypotheses = getWeights(g_i) * sum(alphaHypotheses),
    m = getMatrix(g_i),
    nameHypotheses = nameHypotheses
  ))
}
```

---

## Combining and subsetting graphs {#join-subgraph}

### Joining two graphs

```r
g1 <- BonferroniHolm(2)
g2 <- fixedSequence(2)

# Combine (weights rescaled if needed)
combined <- joinGraphs(g1, g2)
```

### Extracting subgraph

```r
g <- BonferroniHolm(4)

# Keep only H1 and H3
sub <- subgraph(g, c("H1", "H3"))
getWeights(sub)
getMatrix(sub)
```

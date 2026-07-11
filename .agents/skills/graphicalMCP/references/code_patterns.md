# Code Patterns for graphicalMCP

**Note**: These patterns target graphicalMCP (CRAN). See https://openpharma.github.io/graphicalMCP/ for full documentation.

## Table of Contents
1. [Creating graphs with graph_create](#graph-create)
2. [Common graph structures](#common-graphs)
3. [Shortcut (Bonferroni) testing](#shortcut-testing)
4. [Closure testing (Simes, parametric, Hochberg)](#closure-testing)
5. [Updating graphs manually](#graph-update)
6. [Generating intersection weights](#generate-weights)
7. [Rejection orderings](#rejection-orderings)
8. [Power simulation](#power)
9. [Custom success criteria](#success-criteria)
10. [Plotting graphs](#plotting)
11. [Example graphs (built-in)](#example-graphs)
12. [Converting between graph formats](#conversion)
13. [Parametric tests with correlation](#parametric)
14. [Multi-population multi-endpoint designs](#multi-pop)

---

## Creating graphs with graph_create {#graph-create}

`graph_create()` defines a multiplicity graph from hypothesis weights and
a transition matrix.

### Two-hypothesis graph

```r
library(graphicalMCP)

# Simple Bonferroni split: equal weights, pass all alpha on rejection
g <- graph_create(
  hypotheses = c(H1 = 0.5, H2 = 0.5),
  transitions = matrix(c(
    0, 1,
    1, 0
  ), nrow = 2, byrow = TRUE)
)
g
plot(g)
```

### Four-hypothesis graph (two endpoints, two populations)

```r
g <- graph_create(
  hypotheses = c(H1 = 0.25, H2 = 0.25, H3 = 0.25, H4 = 0.25),
  transitions = matrix(c(
    0,   0.5, 0.5, 0,
    0.5, 0,   0,   0.5,
    0,   0.5, 0,   0.5,
    0.5, 0,   0.5, 0
  ), nrow = 4, byrow = TRUE)
)
```

### Weighted graph with asymmetric transitions

```r
# Primary endpoint gets more alpha, secondary gets remainder
g <- graph_create(
  hypotheses = c(H1 = 0.5, H2 = 0.5, H3 = 0, H4 = 0),
  transitions = matrix(c(
    0,   0,   1,   0,
    0,   0,   0,   1,
    0,   1,   0,   0,
    1,   0,   0,   0
  ), nrow = 4, byrow = TRUE)
)
```

---

## Common graph structures {#common-graphs}

### Fixed sequence (gatekeeping)

```r
# H1 must be rejected before testing H2, etc.
g <- graph_create(
  hypotheses = c(H1 = 1, H2 = 0, H3 = 0),
  transitions = matrix(c(
    0, 1, 0,
    0, 0, 1,
    0, 0, 0
  ), nrow = 3, byrow = TRUE)
)
# Or use the built-in:
g <- fixed_sequence(3)
```

### Bonferroni-Holm

```r
# Equal weights, full propagation
g <- bonferroni_holm(3)
# Equivalent to:
g <- graph_create(
  hypotheses = c(H1 = 1/3, H2 = 1/3, H3 = 1/3),
  transitions = matrix(c(
    0,   0.5, 0.5,
    0.5, 0,   0.5,
    0.5, 0.5, 0
  ), nrow = 3, byrow = TRUE)
)
```

### Fallback procedure

```r
# Ordered but with propagation back
g <- fallback(hypotheses = c(H1 = 0.5, H2 = 0.5))
```

### Simple successive (two families)

```r
# Two primary, two secondary with gatekeeping
g <- simple_successive_1()
plot(g)
```

### Two doses, two primary, two secondary

```r
g <- two_doses_two_primary_two_secondary()
plot(g)
```

---

## Shortcut (Bonferroni) testing {#shortcut-testing}

`graph_test_shortcut()` performs the sequentially rejective Bonferroni-based
graphical MCP. This is the standard Maurer-Bretz procedure.

### Basic testing

```r
g <- graph_create(
  hypotheses = c(H1 = 0.5, H2 = 0.5),
  transitions = matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
)

p <- c(H1 = 0.012, H2 = 0.030)

result <- graph_test_shortcut(g, p = p, alpha = 0.025)
result

# Access results
result$outputs$rejected    # Logical vector of rejections
result$outputs$adjusted_p  # Adjusted p-values
```

### With verbose output (intermediate graphs)

```r
result <- graph_test_shortcut(g, p = p, alpha = 0.025, verbose = TRUE)

# See each step of the sequential procedure
result$details
```

### With test values (adjusted significance levels)

```r
result <- graph_test_shortcut(g, p = p, alpha = 0.025, test_values = TRUE)

# Adjusted significance levels at each step
result$test_values
```

---

## Closure testing (Simes, parametric, Hochberg) {#closure-testing}

`graph_test_closure()` performs the full closure principle, allowing more
powerful tests than Bonferroni for correlated endpoints.

### Simes test (for positively correlated test statistics)

```r
g <- graph_create(
  hypotheses = c(H1 = 0.5, H2 = 0.5),
  transitions = matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
)
p <- c(H1 = 0.012, H2 = 0.030)

result <- graph_test_closure(
  g, p = p, alpha = 0.025,
  test_groups = list(1:2),
  test_types = c("simes")
)
result$outputs$rejected
```

### Parametric test (known correlation)

```r
# Correlation between test statistics
corr <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

result <- graph_test_closure(
  g, p = p, alpha = 0.025,
  test_groups = list(1:2),
  test_types = c("parametric"),
  test_corr = list(corr)
)
result$outputs$rejected
```

### Mixed test types (different tests for different groups)

```r
# 4 hypotheses: parametric for H1-H2, Simes for H3-H4
g <- graph_create(
  hypotheses = c(H1 = 0.25, H2 = 0.25, H3 = 0.25, H4 = 0.25),
  transitions = matrix(c(
    0,   0.5, 0.5, 0,
    0.5, 0,   0,   0.5,
    0,   0.5, 0,   0.5,
    0.5, 0,   0.5, 0
  ), nrow = 4, byrow = TRUE)
)

p <- c(H1 = 0.01, H2 = 0.02, H3 = 0.03, H4 = 0.04)

# Correlation for the parametric group
corr12 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

result <- graph_test_closure(
  g, p = p, alpha = 0.025,
  test_groups = list(1:2, 3:4),
  test_types = c("parametric", "simes"),
  test_corr = list(corr12, NA)   # NA for non-parametric tests
)
result$outputs$rejected
```

### Hochberg test

```r
result <- graph_test_closure(
  g, p = p, alpha = 0.025,
  test_groups = list(1:4),
  test_types = c("hochberg")
)
```

### Verbose closure output

```r
result <- graph_test_closure(
  g, p = p, alpha = 0.025,
  test_groups = list(1:4),
  test_types = c("bonferroni"),
  verbose = TRUE
)

# Adjusted p-values for every intersection hypothesis
result$details$adjusted_p
```

---

## Updating graphs manually {#graph-update}

`graph_update()` shows what the graph looks like after rejecting specific
hypotheses. Useful for visualizing alpha propagation.

```r
g <- graph_create(
  hypotheses = c(H1 = 0.5, H2 = 0.5, H3 = 0, H4 = 0),
  transitions = matrix(c(
    0, 0, 1, 0,
    0, 0, 0, 1,
    0, 1, 0, 0,
    1, 0, 0, 0
  ), nrow = 4, byrow = TRUE)
)

# Update after rejecting H1
updated <- graph_update(g, delete = 1)
updated$updated_graph

# Update after rejecting H1 and H3
updated <- graph_update(g, delete = c(1, 3))
updated$updated_graph

# Using logical indexing
updated <- graph_update(g, delete = c(TRUE, FALSE, TRUE, FALSE))

# Plot the update sequence
plot(updated)
```

---

## Generating intersection weights {#generate-weights}

`graph_generate_weights()` generates the weighting strategy for all 2^m - 1
intersection hypotheses in the closure.

```r
g <- graph_create(
  hypotheses = c(H1 = 0.5, H2 = 0.5),
  transitions = matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
)

weights <- graph_generate_weights(g)
weights
# Returns matrix: first m columns = membership, last m columns = weights
#   H1 H2 H1    H2
# 1  1  1  0.50  0.50    # {H1, H2}
# 2  1  0  1.00  0.00    # {H1}
# 3  0  1  0.00  1.00    # {H2}
```

---

## Rejection orderings {#rejection-orderings}

`graph_rejection_orderings()` enumerates all valid orderings in which
hypotheses can be rejected, given the same p-values and graph.

```r
g <- bonferroni_holm(3)
p <- c(H1 = 0.01, H2 = 0.005, H3 = 0.03)

result <- graph_test_shortcut(g, p = p, alpha = 0.025)
orderings <- graph_rejection_orderings(result)
orderings
# Shows all valid sequences of rejections
```

---

## Power simulation {#power}

`graph_calculate_power()` simulates power using multivariate normal
test statistics. Requires marginal power for each hypothesis and a
correlation matrix.

### Basic power simulation

```r
g <- graph_create(
  hypotheses = c(H1 = 0.5, H2 = 0.5),
  transitions = matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
)

# Marginal power (what each test achieves individually at full alpha)
power_marginal <- c(H1 = 0.9, H2 = 0.8)

# Correlation between test statistics
sim_corr <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

pwr <- graph_calculate_power(
  graph = g,
  alpha = 0.025,
  power_marginal = power_marginal,
  sim_corr = sim_corr,
  sim_n = 1e5
)
pwr

# Key outputs
pwr$power$power_local        # Local power for each hypothesis
pwr$power$power_at_least_1   # Power to reject at least one
pwr$power$power_all          # Power to reject all
pwr$power$rejection_expected # Expected number of rejections
```

### Power with Simes test

```r
pwr <- graph_calculate_power(
  graph = g,
  alpha = 0.025,
  power_marginal = power_marginal,
  test_groups = list(1:2),
  test_types = c("simes"),
  sim_corr = sim_corr,
  sim_n = 1e5
)
```

### Power with parametric test

```r
# Test correlation (for the parametric test itself)
test_corr <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

pwr <- graph_calculate_power(
  graph = g,
  alpha = 0.025,
  power_marginal = power_marginal,
  test_groups = list(1:2),
  test_types = c("parametric"),
  test_corr = list(test_corr),
  sim_corr = sim_corr,
  sim_n = 1e5
)
```

### Power with mixed test types

```r
g4 <- graph_create(
  hypotheses = c(H1 = 0.25, H2 = 0.25, H3 = 0.25, H4 = 0.25),
  transitions = matrix(c(
    0, 0.5, 0.5, 0,
    0.5, 0, 0, 0.5,
    0, 0.5, 0, 0.5,
    0.5, 0, 0.5, 0
  ), nrow = 4, byrow = TRUE)
)

power_marginal <- c(0.9, 0.8, 0.7, 0.6)
sim_corr <- matrix(0.3, 4, 4)
diag(sim_corr) <- 1
test_corr_12 <- matrix(c(1, 0.5, 0.5, 1), 2, 2)

pwr <- graph_calculate_power(
  graph = g4,
  alpha = 0.025,
  power_marginal = power_marginal,
  test_groups = list(1:2, 3:4),
  test_types = c("parametric", "simes"),
  test_corr = list(test_corr_12, NA),
  sim_corr = sim_corr,
  sim_n = 1e5
)
```

---

## Custom success criteria {#success-criteria}

Define what "success" means beyond rejecting individual hypotheses.

```r
g <- graph_create(
  hypotheses = c(H1 = 0.5, H2 = 0.5, H3 = 0, H4 = 0),
  transitions = matrix(c(
    0, 0, 1, 0,
    0, 0, 0, 1,
    0, 1, 0, 0,
    1, 0, 0, 0
  ), nrow = 4, byrow = TRUE)
)

# Success = reject at least one primary AND at least one secondary
success_criteria <- list(
  "At least 1 primary + 1 secondary" = function(x) {
    (x[1] | x[2]) & (x[3] | x[4])
  },
  "Both primaries" = function(x) {
    x[1] & x[2]
  }
)

pwr <- graph_calculate_power(
  graph = g,
  alpha = 0.025,
  power_marginal = c(0.9, 0.8, 0.7, 0.6),
  sim_corr = diag(4),
  sim_n = 1e5,
  sim_success = success_criteria
)
pwr$power$power_success
```

---

## Plotting graphs {#plotting}

### Basic plot

```r
g <- simple_successive_1()
plot(g)
```

### Customizing layout

```r
# Grid layout with specified dimensions
plot(g, layout = "grid", nrow = 2, ncol = 2)

# Custom vertex positions (igraph-style layout matrix)
layout_matrix <- matrix(c(
  0, 1,   # H1 position (x, y)
  1, 1,   # H2
  0, 0,   # H3
  1, 0    # H4
), ncol = 2, byrow = TRUE)
plot(g, layout = layout_matrix)
```

### Customizing appearance

```r
plot(g,
  v_palette = c("#6baed6", "#cccccc"),  # Active, deleted colors
  precision = 4,                          # Decimal places
  background_color = "white",
  margins = c(1, 1, 1, 1)
)
```

### Controlling edge curvature

```r
# Curve all paired edges
plot(g, edge_curves = c("pairs" = 0.8))

# Curve specific edges
plot(g, edge_curves = c("H1|H3" = 0.5, "H3|H1" = 0.5))
```

### Plotting updated graphs

```r
updated <- graph_update(g, delete = 1)
plot(updated)  # Shows the update sequence
```

### Epsilon edges

```r
# Display very small transition weights as epsilon
g_eps <- graph_create(
  hypotheses = c(H1 = 0.5, H2 = 0.5),
  transitions = matrix(c(0, 1e-5, 1, 0), nrow = 2, byrow = TRUE)
)
plot(g_eps, eps = 1e-4)  # Edges with weight < eps shown as ε
```

---

## Example graphs (built-in) {#example-graphs}

graphicalMCP provides many pre-built example graphs.

```r
# Simple procedures
bonferroni(3)                     # Equal-weight Bonferroni
bonferroni_holm(4)                # Bonferroni-Holm step-down
fixed_sequence(3)                 # Fixed sequence (gatekeeping)
fallback(c(0.5, 0.5))            # Fallback procedure

# Weighted variants
bonferroni_weighted(c(0.6, 0.4))
bonferroni_holm_weighted(c(0.6, 0.3, 0.1))

# Simes-based
hochberg(3)                       # Hochberg step-up
hommel(3)                         # Hommel procedure
sidak(3)                          # Sidak procedure

# Dunnett-type
dunnett_single_step(3)
dunnett_single_step_weighted(c(0.5, 0.3, 0.2))
dunnett_closure_weighted(c(0.5, 0.3, 0.2))

# Complex clinical trial graphs
simple_successive_1()             # 2 primary + 2 secondary
simple_successive_2()             # 2 primary + 2 secondary (variant)
huque_etal()                      # Huque et al. procedure
two_doses_two_primary_two_secondary()
three_doses_two_primary_two_secondary()

# Improved fallback
fallback_improved_1(c(0.5, 0.5))
fallback_improved_2(c(0.5, 0.5), epsilon = 1e-4)

# Random (for testing)
random_graph(4)
```

---

## Converting between graph formats {#conversion}

```r
# graphicalMCP → gMCPLite (legacy)
g <- bonferroni_holm(3)
g_gmcp <- as_graphMCP(g)

# gMCPLite → graphicalMCP
g_back <- as_initial_graph(g_gmcp)

# graphicalMCP → igraph (for network analysis)
g_ig <- as_igraph(g)

# igraph → graphicalMCP
g_back2 <- as_initial_graph(g_ig)
```

---

## Parametric tests with correlation {#parametric}

When test statistics have known correlation (e.g., nested populations, shared
control arm), parametric tests can be more powerful than Bonferroni.

### Two endpoints with known correlation

```r
g <- graph_create(
  hypotheses = c(H1 = 0.5, H2 = 0.5),
  transitions = matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
)

p <- c(H1 = 0.015, H2 = 0.020)

# Correlation from shared patients (e.g., co-primary endpoints)
test_corr <- matrix(c(1, 0.7, 0.7, 1), nrow = 2)

# Parametric closure test
result <- graph_test_closure(
  g, p = p, alpha = 0.025,
  test_groups = list(1:2),
  test_types = c("parametric"),
  test_corr = list(test_corr)
)
result$outputs$rejected
```

### Nested populations (subgroup + overall)

```r
# H1 = subgroup PFS, H2 = overall PFS
# H3 = subgroup OS, H4 = overall OS
g <- graph_create(
  hypotheses = c(H1 = 0.25, H2 = 0.25, H3 = 0.25, H4 = 0.25),
  transitions = matrix(c(
    0,   0.5, 0.5, 0,
    0.5, 0,   0,   0.5,
    0,   0.5, 0,   0.5,
    0.5, 0,   0.5, 0
  ), nrow = 4, byrow = TRUE)
)

p <- c(H1 = 0.008, H2 = 0.012, H3 = 0.02, H4 = 0.04)

# Correlation structure (nested populations are correlated)
prevalence <- 0.4
corr_sub_full <- sqrt(prevalence)  # ~0.632

test_corr <- matrix(c(
  1,              corr_sub_full, 0.5,            0.3,
  corr_sub_full,  1,             0.3,            0.5,
  0.5,            0.3,           1,              corr_sub_full,
  0.3,            0.5,           corr_sub_full,  1
), nrow = 4, byrow = TRUE)

result <- graph_test_closure(
  g, p = p, alpha = 0.025,
  test_groups = list(1:4),
  test_types = c("parametric"),
  test_corr = list(test_corr)
)
result$outputs$rejected
```

---

## Multi-population multi-endpoint designs {#multi-pop}

### Oncology trial: PFS + OS in subgroup + overall

```r
# Typical oncology multiplicity graph
# H1: PFS in BM+, H2: PFS in ITT
# H3: OS in BM+, H4: OS in ITT

g <- graph_create(
  hypotheses = c(H1 = 0.5, H2 = 0.5, H3 = 0, H4 = 0),
  transitions = matrix(c(
    0,   0,   1,   0,
    0,   0,   0,   1,
    0,   0.5, 0,   0.5,
    0.5, 0,   0.5, 0
  ), nrow = 4, byrow = TRUE)
)

plot(g, layout = matrix(c(0, 1, 1, 1, 0, 0, 1, 0), ncol = 2, byrow = TRUE))
```

### Testing with sequential p-values from group sequential designs

```r
# Sequential p-values from gsDesign2 (one per hypothesis)
p <- c(H1 = 0.008, H2 = 0.015, H3 = 0.020, H4 = 0.045)

result <- graph_test_shortcut(g, p = p, alpha = 0.025)
result$outputs$rejected
result$outputs$adjusted_p
```

### Power for the full design

```r
# Marginal power for each hypothesis
power_marginal <- c(H1 = 0.90, H2 = 0.85, H3 = 0.60, H4 = 0.55)

# Simulation correlation (from information fraction and overlap)
sim_corr <- matrix(c(
  1.0,  0.6,  0.3, 0.2,
  0.6,  1.0,  0.2, 0.3,
  0.3,  0.2,  1.0, 0.6,
  0.2,  0.3,  0.6, 1.0
), nrow = 4, byrow = TRUE)

pwr <- graph_calculate_power(
  graph = g,
  alpha = 0.025,
  power_marginal = power_marginal,
  sim_corr = sim_corr,
  sim_n = 1e5,
  sim_success = list(
    "Reject H1 or H2" = function(x) x[1] | x[2],
    "Reject any OS" = function(x) x[3] | x[4],
    "Reject H1 and H3" = function(x) x[1] & x[3]
  )
)
pwr$power$power_local
pwr$power$power_success
```

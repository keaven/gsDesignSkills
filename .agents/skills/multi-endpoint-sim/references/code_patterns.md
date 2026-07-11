# Code Patterns for Multi-Endpoint Group Sequential Trial Simulation

## Table of Contents
1. [Test statistic functions](#test-functions)
2. [Simulation loop](#simulation-loop)
3. [Sequential testing loop](#sequential-testing)
4. [Cumulative rejection probabilities](#rejection-summary)
5. [Spending time computation](#spending-time)
6. [Correlation matrix](#correlation)

---

## Test statistic functions {#test-functions}

### Logrank Z via simtrial::wlr

Map illness-death ADTTE columns to `wlr()` format. Positive Z = experimental better.

```r
logrank_z <- function(data) {
  d <- data.frame(
    tte = data$AVAL,
    event = 1L - data$CNSR,
    stratum = data$STRATUM,
    treatment = data$TRT
  )
  if (sum(d$event) < 2) return(NA_real_)
  wlr(d, weight = fh(rho = 0, gamma = 0))$z
}
```

### Binomial Z via gsDesign::testBinomial

**Important**: Pass experimental as `x1` and control as `x2` so positive Z = experimental better.

```r
rd_z <- function(data) {
  n_ctrl <- sum(data$TRT == "control")
  n_exp <- sum(data$TRT == "experimental")
  if (n_exp == 0 || n_ctrl == 0) return(NA_real_)
  x_ctrl <- sum(data$AVAL[data$TRT == "control"])
  x_exp <- sum(data$AVAL[data$TRT == "experimental"])
  testBinomial(x1 = x_exp, x2 = x_ctrl, n1 = n_exp, n2 = n_ctrl)
}
```

### Sign convention summary

| Function | Positive Z means | Notes |
|----------|-----------------|-------|
| `simtrial::wlr()` | Experimental better | Standard logrank convention |
| `gsDesign::testBinomial(x1=exp, x2=ctrl)` | Experimental better | x1/x2 order matters |
| `gsDesign::sequentialPValue()` | Expects positive Z = favorable | Matches both above |

---

## Simulation loop {#simulation-loop}

### Structure

Each simulation generates one trial, cuts data at each analysis time, and returns
Z-statistics plus event counts.

```r
run_one_sim <- function(sim_id) {
  tryCatch({
    sim_data <- sim_illness_death(
      n = n_total,
      stratum = data.frame(stratum = "All", p = 1),
      block = c("control", "control", "experimental", "experimental"),
      enroll_rate = enroll_rate_sim,
      transition_rate = transition_rate_pw_sim,
      dropout_rate = dropout_os
    )

    # Cut at each analysis time
    adtte_18 <- cut_illness_death(sim_data, cut_date = 18)
    adtte_24 <- cut_illness_death(sim_data, cut_date = 24)
    adtte_30 <- cut_illness_death(sim_data, cut_date = 30)
    adtte_42 <- cut_illness_death(sim_data, cut_date = 42)

    # Z statistics (8 values: PFS at 3 times, OS at 4 times, ORR at 1 time)
    z_vals <- c(
      logrank_z(adtte_18[adtte_18$PARAMCD == "PFS", ]),
      logrank_z(adtte_18[adtte_18$PARAMCD == "OS", ]),
      rd_z(adtte_24[adtte_24$PARAMCD == "ORR", ]),
      logrank_z(adtte_24[adtte_24$PARAMCD == "PFS", ]),
      logrank_z(adtte_24[adtte_24$PARAMCD == "OS", ]),
      logrank_z(adtte_30[adtte_30$PARAMCD == "PFS", ]),
      logrank_z(adtte_30[adtte_30$PARAMCD == "OS", ]),
      logrank_z(adtte_42[adtte_42$PARAMCD == "OS", ])
    )

    # Event counts for spending time computation (7 values)
    ev_vals <- c(
      sum(1L - adtte_18$CNSR[adtte_18$PARAMCD == "PFS"]),
      sum(1L - adtte_24$CNSR[adtte_24$PARAMCD == "PFS"]),
      sum(1L - adtte_30$CNSR[adtte_30$PARAMCD == "PFS"]),
      sum(1L - adtte_18$CNSR[adtte_18$PARAMCD == "OS"]),
      sum(1L - adtte_24$CNSR[adtte_24$PARAMCD == "OS"]),
      sum(1L - adtte_30$CNSR[adtte_30$PARAMCD == "OS"]),
      sum(1L - adtte_42$CNSR[adtte_42$PARAMCD == "OS"])
    )

    c(z_vals, ev_vals)
  }, error = function(e) rep(NA_real_, 15))
}

# Run in parallel
z_list <- mclapply(seq_len(n_sim), run_one_sim,
                   mc.cores = n_workers, mc.set.seed = TRUE)
sim_results <- do.call(rbind, z_list)

# Split into Z-matrix and events-matrix
z_matrix <- sim_results[, 1:8]
events_matrix <- sim_results[, 9:15]
```

### Key points

- Return a fixed-length numeric vector from each simulation (Z-stats + event counts)
- Use `tryCatch()` to return `NA` for failed simulations
- Collect event counts alongside Z-statistics for spending time computation
- Use `mclapply()` with `mc.set.seed = TRUE` for reproducibility

---

## Sequential testing loop {#sequential-testing}

For each simulated trial, compute sequential p-values at each analysis and
test with `graph_test_shortcut()`.

```r
# Planned events from designs
planned_os_ev <- os_design$n.I    # 4 values
planned_pfs_ev <- pfs_power$n.I   # 3 values

# Track first rejection analysis
first_rej <- matrix(NA_integer_, nrow = n_sim, ncol = 3,
                    dimnames = list(NULL, c("OS", "PFS", "ORR")))

for (i in seq_len(n_sim)) {
  if (any(is.na(z_matrix[i, ]))) next

  z_pfs <- as.numeric(z_matrix[i, c("Z_PFS_18", "Z_PFS_24", "Z_PFS_30")])
  z_os <- as.numeric(z_matrix[i, c("Z_OS_18", "Z_OS_24", "Z_OS_30", "Z_OS_42")])
  z_orr <- as.numeric(z_matrix[i, "Z_ORR_24"])

  ev_pfs <- as.numeric(events_matrix[i, c("PFS_18", "PFS_24", "PFS_30")])
  ev_os <- as.numeric(events_matrix[i, c("OS_18", "OS_24", "OS_30", "OS_42")])

  # ORR nominal p-value
  p_orr_nom <- pnorm(-z_orr)

  rejected <- c(OS = FALSE, PFS = FALSE, ORR = FALSE)
  p_pfs_cur <- 1
  p_orr_cur <- 1

  for (k in 1:4) {
    # OS sequential p-value (all 4 analyses)
    st_os <- pmin(planned_os_ev[1:k], ev_os[1:k]) / planned_os_ev[4]
    if (k == 4) st_os[4] <- 1
    p_os_cur <- tryCatch(
      sequentialPValue(os_design, n.I = ev_os[1:k],
                       Z = z_os[1:k], usTime = st_os),
      error = function(e) 1
    )

    # PFS sequential p-value (3 analyses)
    if (k <= 3) {
      st_pfs <- pmin(planned_pfs_ev[1:k], ev_pfs[1:k]) / planned_pfs_ev[3]
      if (k == 3) st_pfs[3] <- 1
      p_pfs_cur <- tryCatch(
        sequentialPValue(pfs_power, n.I = ev_pfs[1:k],
                         Z = z_pfs[1:k], usTime = st_pfs),
        error = function(e) 1
      )
    }
    # else: carry forward p_pfs_cur from k=3

    # ORR: nominal p-value from IA2 onwards
    if (k >= 2) p_orr_cur <- p_orr_nom

    # Graphical multiplicity testing
    result <- graph_test_shortcut(g,
      p = c(p_os_cur, p_pfs_cur, p_orr_cur),
      alpha = alpha_total
    )

    # Record first rejection
    rej_now <- result$outputs$rejected
    newly <- rej_now & !rejected
    if (any(newly)) first_rej[i, which(newly)] <- k
    rejected <- rejected | rej_now

    if (all(rejected)) break
  }
}
```

### Key patterns

- **Sequential p-values decrease** as evidence accumulates; `graph_test_shortcut()` applies cumulatively
- **Carry forward**: PFS p-value from analysis 3 is carried to analysis 4 (PFS has only 3 analyses)
- **ORR starts at IA2**: Before IA2, ORR p-value = 1 (not tested yet)
- **Early exit**: Break out of the loop once all hypotheses are rejected
- **Error handling**: `tryCatch()` around `sequentialPValue()` to handle edge cases

---

## Spending time computation {#spending-time}

Spending time controls how alpha is allocated across analyses. The formula:

$$\text{spending time}_k = \frac{\min(\text{planned events}_k, \text{actual events}_k)}{\text{planned final events}}$$

At the final analysis, spending time is forced to 1.

```r
# OS spending time (4 analyses)
st_os <- pmin(planned_os_ev[1:k], ev_os[1:k]) / planned_os_ev[4]
if (k == 4) st_os[4] <- 1  # Force 1 at final

# PFS spending time (3 analyses)
st_pfs <- pmin(planned_pfs_ev[1:k], ev_pfs[1:k]) / planned_pfs_ev[3]
if (k == 3) st_pfs[3] <- 1  # Force 1 at final
```

### Why use pmin?

Using `pmin` prevents over-spending when actual events exceed planned. If the
simulation produces more events than the design assumed, we cap spending at
the planned information fraction to maintain Type I error control.

---

## Cumulative rejection probabilities {#rejection-summary}

Build a table showing cumulative rejections by analysis from the `first_rej` matrix.

```r
analysis_labs <- c("IA1 (18 mo)", "IA2 (24 mo)", "IA3 (30 mo)", "FA (42 mo)")

rej_tab <- data.frame(Hypothesis = c("OS", "PFS", "ORR"))
for (k in 1:4) {
  rej_tab[[analysis_labs[k]]] <- sapply(1:3, function(h) {
    sum(first_rej[, h] <= k, na.rm = TRUE)
  })
}

# Display as gt table
rej_tab |>
  gt() |>
  tab_header(
    title = paste("Cumulative Rejections by Analysis (", n_sim, "simulations)"),
    subtitle = paste0("Simulation HR: OS = ", sim_hr_os, ", PFS = ", sim_hr_pfs)
  )
```

### Interpreting results

- Each cell shows the number of simulations where the hypothesis was rejected
  by that analysis (inclusive of earlier analyses)
- Divide by `n_sim` for cumulative rejection probability
- The final column gives the overall power for each hypothesis
- Under the global null (all HRs = 1), the maximum rejection probability for
  any hypothesis should be ≤ α_total (FWER control)

---

## Correlation matrix of test statistics {#correlation}

The empirical correlation between Z-statistics across endpoint-analysis
combinations characterizes the joint distribution.

```r
col_order <- c("Z_ORR_24", "Z_PFS_18", "Z_PFS_24", "Z_PFS_30",
               "Z_OS_18", "Z_OS_24", "Z_OS_30", "Z_OS_42")
z_ordered <- z_matrix[, col_order]
cor_matrix <- cor(z_ordered, use = "pairwise.complete.obs")
```

### Expected correlation structure

- **Same endpoint, different analyses**: High correlation (0.5–0.9) reflecting
  incremental information (e.g., Z_OS_18 and Z_OS_24)
- **PFS and OS**: Moderate positive correlation (0.3–0.6) from the illness-death
  model coupling
- **ORR and PFS/OS**: Weak positive correlation (0.1–0.3) from shared response state
- **Theoretical**: For the same endpoint at analyses with $n_1$ and $n_2$ events,
  correlation $\approx \sqrt{n_1 / n_2}$

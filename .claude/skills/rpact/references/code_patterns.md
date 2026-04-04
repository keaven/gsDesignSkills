# Code Patterns for rpact

**Note**: These patterns target rpact >= 4.x (CRAN). See https://docs.rpact.org for full documentation.

## Table of Contents
1. [Group sequential designs](#group-sequential)
2. [Adaptive designs (inverse normal and Fisher)](#adaptive-designs)
3. [Sample size for survival endpoints](#ss-survival)
4. [Sample size for means and rates](#ss-means-rates)
5. [Power computation](#power)
6. [Simulation for survival](#sim-survival)
7. [Multi-arm designs and simulation](#multi-arm)
8. [Enrichment designs](#enrichment)
9. [Entering observed data with getDataset](#dataset)
10. [Stage results and analysis](#analysis)
11. [Conditional power and sample size reassessment](#conditional-power)
12. [Adjusted inference (p-values and CIs)](#inference)
13. [Piecewise survival and accrual specification](#piecewise)
14. [Comparing designs with getDesignSet](#design-set)
15. [Plotting](#plotting)
16. [Reproducible code with getObjectRCode](#reproducible)

---

## Group sequential designs {#group-sequential}

`getDesignGroupSequential()` creates classical group sequential designs
with spending functions or fixed boundary types.

### O'Brien-Fleming design (default)

```r
library(rpact)

design <- getDesignGroupSequential(
  kMax = 3,
  alpha = 0.025,
  beta = 0.2,
  sided = 1,
  typeOfDesign = "OF"
)
design
```

### Pocock design

```r
design <- getDesignGroupSequential(
  kMax = 3,
  alpha = 0.025,
  beta = 0.2,
  sided = 1,
  typeOfDesign = "P"
)
```

### Alpha-spending (Lan-DeMets O'Brien-Fleming)

```r
design <- getDesignGroupSequential(
  kMax = 3,
  alpha = 0.025,
  beta = 0.2,
  sided = 1,
  informationRates = c(0.5, 0.75, 1),
  typeOfDesign = "asOF"   # Alpha-spending OF approximation
)
```

### Alpha-spending with Hwang-Shih-DeCani (gamma)

```r
design <- getDesignGroupSequential(
  kMax = 3,
  alpha = 0.025,
  sided = 1,
  informationRates = c(0.5, 0.75, 1),
  typeOfDesign = "asHSD",
  gammaA = -4              # Gamma parameter (negative = conservative early)
)
```

### Non-binding futility with beta-spending

```r
design <- getDesignGroupSequential(
  kMax = 3,
  alpha = 0.025,
  beta = 0.2,
  sided = 1,
  informationRates = c(0.5, 0.75, 1),
  typeOfDesign = "asOF",
  typeBetaSpending = "bsOF",   # Beta-spending for futility
  bindingFutility = FALSE
)
```

### Custom futility bounds

```r
design <- getDesignGroupSequential(
  kMax = 3,
  alpha = 0.025,
  beta = 0.2,
  sided = 1,
  informationRates = c(0.5, 0.75, 1),
  typeOfDesign = "asOF",
  futilityBounds = c(0, 0.5),  # Z-scale futility bounds at IA1, IA2
  bindingFutility = FALSE
)
```

### Design characteristics

```r
# Examine properties of any design
chars <- getDesignCharacteristics(design)
chars  # Shows power, ASN, stopping probabilities
```

---

## Adaptive designs (inverse normal and Fisher) {#adaptive-designs}

### Inverse normal combination test

```r
# For adaptive designs where sample size can be recalculated mid-trial
design_in <- getDesignInverseNormal(
  kMax = 2,
  alpha = 0.025,
  beta = 0.2,
  sided = 1,
  typeOfDesign = "asOF"
)
```

### Fisher combination test

```r
design_fisher <- getDesignFisher(
  kMax = 2,
  alpha = 0.025,
  method = "equalAlpha",    # Equal alpha at each stage
  alpha0Vec = 0.5           # Futility bound (p-value scale)
)
```

### Fisher with no interaction method

```r
design_fisher <- getDesignFisher(
  kMax = 3,
  alpha = 0.025,
  method = "noInteraction",
  alpha0Vec = c(0.5, 0.5)
)
```

---

## Sample size for survival endpoints {#ss-survival}

### Basic survival design (hazard ratio)

```r
design <- getDesignGroupSequential(
  kMax = 3, alpha = 0.025, sided = 1,
  typeOfDesign = "asOF"
)

ss <- getSampleSizeSurvival(
  design = design,
  hazardRatio = 0.7,
  accrualTime = c(0, 12),
  accrualIntensity = 30,      # 30 patients/month
  followUpTime = 12,
  dropoutRate1 = 0.025,
  dropoutRate2 = 0.025
)
ss
```

### Using median survival times

```r
ss <- getSampleSizeSurvival(
  design = design,
  median1 = 18,               # Experimental arm median
  median2 = 12,               # Control arm median
  accrualTime = c(0, 6, 12),
  accrualIntensity = c(15, 30),
  followUpTime = 18,
  dropoutRate1 = 0.025,
  dropoutRate2 = 0.025
)
```

### Using event probabilities

```r
ss <- getSampleSizeSurvival(
  design = design,
  pi1 = 0.35,                 # 12-month event rate, experimental
  pi2 = 0.50,                 # 12-month event rate, control
  eventTime = 12,             # Reference time for pi1, pi2
  accrualTime = c(0, 12),
  accrualIntensity = 30,
  followUpTime = 12
)
```

### Piecewise exponential survival

```r
ss <- getSampleSizeSurvival(
  design = design,
  piecewiseSurvivalTime = c(0, 6),
  lambda2 = c(0.04, 0.03),   # Control hazard rates per period
  hazardRatio = 0.7,
  accrualTime = c(0, 12),
  accrualIntensity = 30,
  followUpTime = 18
)
```

### Fixed design (no interim)

```r
ss_fixed <- getSampleSizeSurvival(
  hazardRatio = 0.7,
  accrualTime = c(0, 12),
  accrualIntensity = 30,
  followUpTime = 12
)
```

---

## Sample size for means and rates {#ss-means-rates}

### Continuous endpoint (two-sample)

```r
design <- getDesignGroupSequential(
  kMax = 2, alpha = 0.025, sided = 1,
  typeOfDesign = "asOF"
)

ss_means <- getSampleSizeMeans(
  design = design,
  groups = 2,
  alternative = 0.5,          # Mean difference
  stDev = 1,
  normalApproximation = TRUE
)
ss_means
```

### One-sample mean

```r
ss_one <- getSampleSizeMeans(
  design = design,
  groups = 1,
  alternative = 0.3,
  stDev = 1
)
```

### Binary endpoint (two-sample)

```r
ss_rates <- getSampleSizeRates(
  design = design,
  groups = 2,
  pi1 = 0.40,                 # Experimental response rate
  pi2 = 0.20,                 # Control response rate
  normalApproximation = TRUE
)
ss_rates
```

### Binary endpoint with risk ratio

```r
ss_rr <- getSampleSizeRates(
  design = design,
  groups = 2,
  riskRatio = TRUE,
  thetaH0 = 1,                # H0: RR = 1
  pi1 = 0.40,
  pi2 = 0.20
)
```

---

## Power computation {#power}

### Power for survival endpoint

```r
design <- getDesignGroupSequential(
  kMax = 3, alpha = 0.025, sided = 1,
  typeOfDesign = "asOF"
)

pwr <- getPowerSurvival(
  design = design,
  hazardRatio = c(0.6, 0.7, 0.8),  # Multiple HR values
  maxNumberOfEvents = 300,
  maxNumberOfSubjects = 500,
  accrualTime = c(0, 12),
  accrualIntensity = 42,
  followUpTime = 12
)
pwr
```

### Power for means

```r
pwr_means <- getPowerMeans(
  design = design,
  alternative = seq(0.2, 0.8, 0.1),
  stDev = 1,
  maxNumberOfSubjects = 200
)
```

### Power for rates

```r
pwr_rates <- getPowerRates(
  design = design,
  pi1 = seq(0.3, 0.6, 0.05),
  pi2 = 0.2,
  maxNumberOfSubjects = 200
)
```

---

## Simulation for survival {#sim-survival}

### Two-arm survival simulation

```r
design <- getDesignGroupSequential(
  kMax = 3, alpha = 0.025, sided = 1,
  typeOfDesign = "asOF",
  typeBetaSpending = "bsOF",
  bindingFutility = FALSE
)

sim <- getSimulationSurvival(
  design = design,
  hazardRatio = 0.7,
  plannedEvents = c(100, 200, 300),
  maxNumberOfSubjects = 500,
  accrualTime = c(0, 12),
  accrualIntensity = 42,
  followUpTime = 18,
  maxNumberOfIterations = 10000,
  seed = 12345
)
sim
summary(sim)
```

### Simulation under multiple scenarios

```r
sim <- getSimulationSurvival(
  design = design,
  hazardRatio = c(0.6, 0.7, 0.8, 1.0),  # Multiple HR values
  plannedEvents = c(100, 200, 300),
  maxNumberOfSubjects = 500,
  accrualTime = c(0, 12),
  accrualIntensity = 42,
  followUpTime = 18,
  maxNumberOfIterations = 10000,
  seed = 12345
)
```

### Simulation with piecewise survival (delayed effect)

```r
sim_nph <- getSimulationSurvival(
  design = design,
  piecewiseSurvivalTime = c(0, 6),
  lambda2 = c(0.04, 0.03),
  lambda1 = c(0.04, 0.03 * 0.6),   # HR = 1 early, 0.6 after 6 months
  plannedEvents = c(100, 200, 300),
  maxNumberOfSubjects = 500,
  accrualTime = c(0, 12),
  accrualIntensity = 42,
  maxNumberOfIterations = 10000,
  seed = 12345
)
```

---

## Multi-arm designs and simulation {#multi-arm}

### Multi-arm survival simulation (Dunnett)

```r
design_in <- getDesignInverseNormal(
  kMax = 2,
  alpha = 0.025,
  sided = 1,
  typeOfDesign = "asOF"
)

sim_ma <- getSimulationMultiArmSurvival(
  design = design_in,
  activeArms = 3,
  plannedEvents = c(50, 100),
  intersectionTest = "Dunnett",
  typeOfSelection = "best",        # Select best arm at interim
  effectMeasure = "effectEstimate",
  successCriterion = "atLeastOne",
  hazardRatio = c(0.6, 0.7, 0.8),  # HR for each active arm
  maxNumberOfSubjects = 400,
  accrualTime = c(0, 12),
  accrualIntensity = 40,
  maxNumberOfIterations = 10000,
  seed = 12345
)
sim_ma
```

### Multi-arm with dose-response (means)

```r
sim_dose <- getSimulationMultiArmMeans(
  design = design_in,
  activeArms = 4,
  typeOfShape = "sigmoidEmax",
  muMaxVector = seq(0.2, 1, 0.2),   # Maximum effect sizes
  gED50 = 2,
  intersectionTest = "Dunnett",
  typeOfSelection = "best",
  stDev = 1,
  maxNumberOfSubjects = 300,
  maxNumberOfIterations = 10000
)
```

### Multi-arm with arm dropping

```r
sim_drop <- getSimulationMultiArmSurvival(
  design = design_in,
  activeArms = 3,
  plannedEvents = c(50, 100),
  intersectionTest = "Simes",
  typeOfSelection = "epsilon",     # Drop if within epsilon of best
  epsilonValue = 0.1,
  hazardRatio = c(0.6, 0.7, 0.8),
  maxNumberOfSubjects = 400,
  accrualTime = c(0, 12),
  accrualIntensity = 40,
  maxNumberOfIterations = 10000
)
```

---

## Enrichment designs {#enrichment}

### Enrichment survival simulation

```r
design_in <- getDesignInverseNormal(
  kMax = 2,
  alpha = 0.025,
  sided = 1,
  typeOfDesign = "asOF"
)

# Effect sizes: rows = scenarios, columns = populations (S, F\S)
effectList <- list(
  subGroups = c("S", "R"),         # Subgroup and complement
  prevalences = c(0.4, 0.6),
  effects = matrix(c(
    0.7, 1.0,     # Scenario 1: effect only in subgroup
    0.7, 0.8,     # Scenario 2: larger effect in subgroup
    0.7, 0.7      # Scenario 3: equal effect
  ), ncol = 2, byrow = TRUE)
)

sim_enrich <- getSimulationEnrichmentSurvival(
  design = design_in,
  effectList = effectList,
  intersectionTest = "Simes",
  typeOfSelection = "best",
  successCriterion = "atLeastOne",
  plannedEvents = c(80, 160),
  maxNumberOfSubjects = 400,
  maxNumberOfIterations = 10000,
  seed = 12345
)
```

---

## Entering observed data with getDataset {#dataset}

### Survival data (cumulative events and logrank statistics)

```r
# Stage-wise observed data
data_surv <- getDataset(
  cumulativeEvents = c(50, 100, 150),
  cumulativeLogRanks = c(1.2, 1.8, 2.3)
)
```

### Means data (two-arm)

```r
data_means <- getDataset(
  n1 = c(30, 30),             # Stage-wise sample sizes, arm 1
  n2 = c(30, 30),             # Stage-wise sample sizes, arm 2
  means1 = c(1.5, 1.8),      # Stage-wise means, arm 1
  means2 = c(1.0, 1.1),      # Stage-wise means, arm 2
  stDevs1 = c(2.0, 2.1),     # Stage-wise SDs, arm 1
  stDevs2 = c(1.9, 2.0)      # Stage-wise SDs, arm 2
)
```

### Rates data (two-arm)

```r
data_rates <- getDataset(
  n1 = c(50, 50),
  n2 = c(50, 50),
  events1 = c(20, 25),
  events2 = c(10, 15)
)
```

### Reading from CSV

```r
data_csv <- readDataset("trial_data.csv")
```

---

## Stage results and analysis {#analysis}

### Computing stage results

```r
design <- getDesignGroupSequential(
  kMax = 3, alpha = 0.025, sided = 1,
  typeOfDesign = "asOF"
)

data <- getDataset(
  cumulativeEvents = c(50, 100),
  cumulativeLogRanks = c(1.5, 2.1)
)

stageResults <- getStageResults(design, data)
stageResults
```

### Full analysis results

```r
results <- getAnalysisResults(
  design = design,
  dataInput = data,
  directionUpper = FALSE        # Lower is better (HR < 1)
)
results
```

### Checking test actions

```r
# What decision at each stage?
getTestActions(results)
# Returns: "continue", "reject", or "accept" (futility)
```

---

## Conditional power and sample size reassessment {#conditional-power}

### Conditional power at interim

```r
design <- getDesignInverseNormal(
  kMax = 2, alpha = 0.025, sided = 1,
  typeOfDesign = "asOF"
)

data <- getDataset(
  n1 = c(50),
  n2 = c(50),
  means1 = c(1.5),
  means2 = c(1.0),
  stDevs1 = c(2.0),
  stDevs2 = c(1.9)
)

stageResults <- getStageResults(design, data)

# Conditional power for planned n
cp <- getConditionalPower(
  stageResults,
  nPlanned = c(100)           # Planned total n for stage 2
)
cp
```

### Conditional rejection probabilities (for SSR)

```r
crp <- getConditionalRejectionProbabilities(stageResults)
crp

# Sample size reassessment: increase stage 2 n to achieve target CP
# Use conditional power to determine new n
```

---

## Adjusted inference (p-values and CIs) {#inference}

### Final p-value

```r
design <- getDesignGroupSequential(
  kMax = 3, alpha = 0.025, sided = 1,
  typeOfDesign = "asOF"
)

data <- getDataset(
  n1 = c(30, 30, 30),
  n2 = c(30, 30, 30),
  means1 = c(1.5, 1.8, 2.0),
  means2 = c(1.0, 1.1, 1.0),
  stDevs1 = c(2.0, 2.1, 1.9),
  stDevs2 = c(1.9, 2.0, 2.0)
)

stageResults <- getStageResults(design, data)

# Stage-wise adjusted p-value
pval <- getFinalPValue(stageResults)
pval
```

### Final confidence interval

```r
ci <- getFinalConfidenceInterval(
  design = design,
  dataInput = data,
  directionUpper = TRUE
)
ci
```

### Repeated p-values and confidence intervals

```r
results <- getAnalysisResults(design, data, directionUpper = TRUE)

# Repeated p-values (available at each stage)
results$repeatedPValues

# Repeated CIs
results$repeatedConfidenceIntervalLowerBounds
results$repeatedConfidenceIntervalUpperBounds
```

---

## Piecewise survival and accrual specification {#piecewise}

### Accrual time with ramp-up

```r
# Vector specification
accrual <- getAccrualTime(
  accrualTime = c(0, 3, 6, 12),
  accrualIntensity = c(10, 20, 30),   # Patients/month per period
  maxNumberOfSubjects = 500
)

# List specification (more readable)
accrual <- getAccrualTime(
  list(
    "0 - <3"  = 10,
    "3 - <6"  = 20,
    "6 - 12"  = 30
  ),
  maxNumberOfSubjects = 500
)
```

### Piecewise exponential survival

```r
# Delayed treatment effect
pw_surv <- getPiecewiseSurvivalTime(
  piecewiseSurvivalTime = c(0, 6),
  lambda2 = c(0.04, 0.03),            # Control hazard
  lambda1 = c(0.04, 0.03 * 0.6)       # Experimental hazard (HR = 1, then 0.6)
)

# Using hazard ratio (applied to all periods)
pw_surv <- getPiecewiseSurvivalTime(
  piecewiseSurvivalTime = c(0, 6),
  lambda2 = c(0.04, 0.03),
  hazardRatio = 0.7
)

# List specification
pw_surv <- getPiecewiseSurvivalTime(
  list(
    "0 - <6"  = 0.04,
    ">=6"     = 0.03
  ),
  hazardRatio = 0.7
)
```

### Using piecewise objects in sample size

```r
ss <- getSampleSizeSurvival(
  design = getDesignGroupSequential(kMax = 3, alpha = 0.025, typeOfDesign = "asOF"),
  piecewiseSurvivalTime = c(0, 6),
  lambda2 = c(0.04, 0.03),
  hazardRatio = 0.7,
  accrualTime = c(0, 3, 6, 12),
  accrualIntensity = c(10, 20, 30),
  followUpTime = 18
)
```

---

## Comparing designs with getDesignSet {#design-set}

### Varying a single parameter

```r
# Compare O'Brien-Fleming, Pocock, and HSD designs
ds <- getDesignSet(
  design = getDesignGroupSequential(kMax = 3, alpha = 0.025, sided = 1),
  typeOfDesign = c("OF", "P", "asHSD")
)
plot(ds)
```

### Comparing gamma values

```r
ds <- getDesignSet(
  design = getDesignGroupSequential(
    kMax = 3, alpha = 0.025, sided = 1,
    typeOfDesign = "asHSD"
  ),
  gammaA = c(-4, -2, 0, 2)
)
plot(ds)
```

### From pre-built designs

```r
d1 <- getDesignGroupSequential(kMax = 3, alpha = 0.025, typeOfDesign = "OF")
d2 <- getDesignGroupSequential(kMax = 3, alpha = 0.025, typeOfDesign = "P")
d3 <- getDesignGroupSequential(kMax = 3, alpha = 0.025, typeOfDesign = "asOF")
ds <- getDesignSet(designs = c(d1, d2, d3))
plot(ds)
```

---

## Plotting {#plotting}

### Design boundaries

```r
design <- getDesignGroupSequential(
  kMax = 4, alpha = 0.025, sided = 1,
  typeOfDesign = "asOF",
  typeBetaSpending = "bsOF",
  bindingFutility = FALSE
)

# Default: boundary plot
plot(design)

# Specific plot type
plot(design, type = 1)   # Boundaries (Z-scale)
plot(design, type = 3)   # Stage levels (adjusted significance levels)
plot(design, type = 4)   # Error spending
plot(design, type = 5)   # Power and early stopping
```

### Sample size results

```r
ss <- getSampleSizeSurvival(
  design = design,
  hazardRatio = seq(0.5, 0.9, 0.05),
  accrualTime = c(0, 12),
  accrualIntensity = 30,
  followUpTime = 12
)
plot(ss)
```

### Power curves

```r
pwr <- getPowerSurvival(
  design = design,
  hazardRatio = seq(0.5, 1, 0.05),
  maxNumberOfEvents = 300,
  maxNumberOfSubjects = 500,
  accrualTime = c(0, 12),
  accrualIntensity = 42,
  followUpTime = 12
)
plot(pwr)
```

### Simulation results

```r
sim <- getSimulationSurvival(
  design = design,
  hazardRatio = c(0.6, 0.7, 0.8, 1.0),
  plannedEvents = c(75, 150, 225, 300),
  maxNumberOfSubjects = 500,
  accrualTime = c(0, 12),
  accrualIntensity = 42,
  maxNumberOfIterations = 10000
)
plot(sim)
```

### Available plot types

```r
# Check what plot types are available for any object
getAvailablePlotTypes(design)
getAvailablePlotTypes(ss)
```

---

## Reproducible code with getObjectRCode {#reproducible}

```r
design <- getDesignGroupSequential(
  kMax = 3, alpha = 0.025, sided = 1,
  typeOfDesign = "asOF"
)

# Print reproducible R code
getObjectRCode(design, output = "cat")

# Include default parameters too
getObjectRCode(design, includeDefaultParameters = TRUE, output = "cat")

# Get as character vector
code <- getObjectRCode(design, output = "vector")

# Markdown output (for reports)
getObjectRCode(design, output = "markdown")
```

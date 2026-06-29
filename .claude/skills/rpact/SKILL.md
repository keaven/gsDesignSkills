---
name: rpact
description: >
  Guide users through confirmatory adaptive clinical trial design and analysis
  using the rpact R package. Use this skill when the user asks about: adaptive
  designs, sample size reassessment, conditional power, inverse normal combination
  test, Fisher combination test, multi-stage designs, or rpact design objects.
---

# Adaptive Clinical Trial Design with rpact

**Note**: This skill targets rpact >= 4.4.0 (CRAN). Dev version 4.4.0.9305 is on GitHub.

## API reference

- Full function docs: `references/llms.txt` (source: https://docs.rpact.org/reference/)
- Workflow patterns: `references/code_patterns.md`

## Key functions

### Trial design
- `getDesignGroupSequential()` - Group sequential design (O'Brien-Fleming, Pocock, alpha-spending, etc.)
- `getDesignInverseNormal()` - Inverse normal combination test for adaptive designs
- `getDesignFisher()` - Fisher combination test for adaptive designs
- `getDesignConditionalDunnett()` - Conditional Dunnett test for multi-arm
- `getDesignCharacteristics()` - Design characteristics and properties
- `getDesignSet()` - Compare multiple designs
- `getFutilityBounds()` - Convert futility bounds between scales (z-value, p-value, conditional power, predictive power, reverse conditional power, effect estimate) (new in 4.3.0)

### Sample size
- `getSampleSizeMeans()` - Continuous endpoints
- `getSampleSizeRates()` - Binary endpoints
- `getSampleSizeSurvival()` - Time-to-event endpoints
- `getSampleSizeCounts()` - Count data endpoints

### Power
- `getPowerMeans()` - Continuous endpoints
- `getPowerRates()` - Binary endpoints
- `getPowerSurvival()` - Time-to-event endpoints
- `getPowerCounts()` - Count data endpoints
- `getPowerAndAverageSampleNumber()` - Power and ASN

### Simulation
- `getSimulationMeans()` / `getSimulationRates()` / `getSimulationSurvival()` / `getSimulationCounts()` - Two-arm
- `getSimulationMultiArmMeans()` / `getSimulationMultiArmRates()` / `getSimulationMultiArmSurvival()` - Multi-arm
- `getSimulationEnrichmentMeans()` / `getSimulationEnrichmentRates()` / `getSimulationEnrichmentSurvival()` - Enrichment

### Analysis
- `getDataset()` - Enter stage-wise observed data
- `getStageResults()` - Compute stage-wise test statistics
- `getAnalysisResults()` - Final adaptive analysis results
- `getClosedCombinationTestResults()` - Closed testing for multi-arm
- `getClosedConditionalDunnettTestResults()` - Dunnett closed testing
- `getConditionalPower()` - Conditional power at interim
- `getConditionalRejectionProbabilities()` - CRP for adaptive recalculation

### Inference
- `getFinalPValue()` - Adjusted final p-value
- `getFinalConfidenceInterval()` - Adjusted confidence interval
- `getRepeatedPValues()` - Repeated p-values across stages
- `getRepeatedConfidenceIntervals()` - Repeated CIs across stages

### Utilities
- `getAccrualTime()` - Define accrual period
- `getPiecewiseSurvivalTime()` - Piecewise exponential survival
- `getEventProbabilities()` - Event probability computation
- `getNumberOfSubjects()` - Expected accrual over time
- `getObservedInformationRates()` - Observed information rates
- `getTestActions()` - Decision at each stage
- `getPerformanceScore()` - Design performance evaluation

### Plotting
- `plot()` methods for all design, power, simulation, and analysis objects
- `getAvailablePlotTypes()` - List available plot types for an object

### I/O
- `readDataset()` / `readDatasets()` - Read data from CSV
- `writeDataset()` / `writeDatasets()` - Write data to CSV
- `getObjectRCode()` - Generate reproducible R code for any rpact object

## Workflow patterns

For detailed code templates, read `references/code_patterns.md`.

Topics covered:
- Group sequential designs (OF, Pocock, alpha-spending, HSD, custom futility)
- Adaptive designs (inverse normal combination, Fisher combination)
- Sample size for survival (HR, medians, event probabilities, piecewise)
- Sample size for means and rates (one-sample, two-sample, risk ratio)
- Power computation (survival, means, rates)
- Simulation for survival (two-arm, piecewise/NPH, multiple scenarios)
- Multi-arm designs and simulation (Dunnett, arm selection, dose-response)
- Enrichment designs (subgroup-specific effects)
- Entering observed data with `getDataset()` (survival, means, rates)
- Stage results and analysis (`getStageResults()`, `getAnalysisResults()`, `getTestActions()`)
- Conditional power and sample size reassessment
- Adjusted inference (final p-values, CIs, repeated p-values/CIs)
- Piecewise survival and accrual specification
- Comparing designs with `getDesignSet()`
- Plotting (boundaries, power curves, simulation results)
- Reproducible code with `getObjectRCode()`

## Important design considerations

- **`typeOfDesign` choices**: "OF" (O'Brien-Fleming), "P" (Pocock), "asOF"/"asHSD"/"asKD" (alpha-spending variants). Alpha-spending versions ("as*") allow unequal information rate spacing
- **`getDesignInverseNormal()` vs `getDesignGroupSequential()`**: Same signature, but inverse normal is for adaptive designs where sample size can be recalculated between stages
- **`bindingFutility = FALSE`**: Non-binding futility is standard for confirmatory trials; efficacy bounds computed ignoring futility
- **Survival parameterization**: Can specify treatment effect via `hazardRatio`, `median1`/`median2`, `lambda1`/`lambda2`, or `pi1`/`pi2` (event probabilities) — use only one
- **`plannedEvents`** in simulation: Vector of cumulative event counts at each analysis (not information fractions)
- **Multi-arm `intersectionTest`**: "Dunnett" (parametric, most powerful), "Simes" (non-parametric), "Bonferroni" (conservative)
- **`typeOfSelection`**: Controls arm selection at interim — "best" (highest effect), "rBest" (r best arms), "epsilon" (within epsilon of best), "all" (no dropping)
- **`getObjectRCode()`**: Use to generate reproducible R code for any rpact object — invaluable for reporting and validation
- **`futilityBoundsScale` (v4.3.0+)**: In `getDesignInverseNormal()` and `getDesignGroupSequential()`, specify futility bounds on alternative scales (p-value, conditional power, etc.) that are internally converted
- **`alpha0Scale` (v4.3.0+)**: In `getDesignFisher()`, specify alpha0 bounds on alternative scales
- **`efficacyStops` / `futilityStops` (v4.2.1+)**: Control which stages have efficacy/futility stopping
- **Dose-response (v4.2.0+)**: Multi-arm simulations support `doseLevels` for linear or sigmoid Emax models
- **Unequal variances (v4.2.0+)**: `getSampleSizeMeans()`, `getPowerMeans()`, `getSimulationMeans()` support different variances per group
- **rpact defaults**: `alpha = 0.05`, `beta = 0.2`, `kMax = NA` (inferred), `sided = 1` — always set explicitly for clinical trial designs

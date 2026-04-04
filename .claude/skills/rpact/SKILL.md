---
name: rpact
description: >
  Guide users through confirmatory adaptive clinical trial design and analysis
  using the rpact R package. Use this skill when the user asks about: adaptive
  designs, sample size reassessment, conditional power, inverse normal combination
  test, Fisher combination test, multi-stage designs, or rpact design objects.
---

# Adaptive Clinical Trial Design with rpact

## API reference

For full function documentation (arguments, return values, examples), read `references/llms.txt`.
Source: https://docs.rpact.org/reference/

## Key functions

### Trial design
- `getDesignGroupSequential()` - Group sequential design (O'Brien-Fleming, Pocock, alpha-spending, etc.)
- `getDesignInverseNormal()` - Inverse normal combination test for adaptive designs
- `getDesignFisher()` - Fisher combination test for adaptive designs
- `getDesignConditionalDunnett()` - Conditional Dunnett test for multi-arm
- `getDesignCharacteristics()` - Design characteristics and properties
- `getDesignSet()` - Compare multiple designs

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

---
name: gsDesign
description: >
  Guide users through classical group sequential trial design using the gsDesign
  R package. Use this skill when the user asks about: group sequential boundaries,
  spending functions (sfLDOF, sfHSD, sfPoints), sample size for time-to-event or
  binomial trials, gsDesign objects, or plotting group sequential bounds.
---

# Group Sequential Design with gsDesign

## API reference

For full function documentation (arguments, return values, examples), read `references/llms.txt`.
Source: https://gsDesign.ai/packages/gsdesign/llms.txt

## Key functions

- `gsDesign()` - Core boundary computation for group sequential designs
- `gsSurvCalendar()` / `nSurv()` - Time-to-event sample size and event counts
- `nNormal()` - Sample size for normal endpoints
- `gsBinomialExact()` - Exact binomial group sequential design
- `gsBoundSummary()` - Formatted summary tables of bounds
- `gsProbability()` - Boundary crossing probabilities
- `gsCP()` / `gsBoundCP()` - Conditional power
- `ssrCP()` - Sample size re-estimation based on conditional power
- `plot.gsDesign()` - Plotting group sequential designs
- `hGraph()` - Multiplicity graph visualization
- `toInteger()` - Round sample sizes to integers
- `sequentiaPValue()` - Sequential p-value computation

## Spending functions

- `sfLDOF` - Lan-DeMets approximation to O'Brien-Fleming
- `sfHSD` - Hwang-Shih-DeCani (gamma family)
- `sfPower` - Kim-DeMets (power) spending
- `sfExponential` - Exponential spending
- `sfPoints` - Pointwise (piecewise linear) spending
- `sfLinear` - Linear spending
- `sfTDist` - t-distribution spending
- `sfXG` - Xu-Garden conditional error spending

## Output and reporting

- `as_gt()` - Convert summary tables to gt objects
- `as_rtf()` - Save summary tables as RTF files
- `as_table()` - Create summary table objects
- `xtable()` - LaTeX/HTML table output

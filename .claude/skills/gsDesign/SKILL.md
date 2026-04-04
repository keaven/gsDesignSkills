---
name: gsDesign
description: >
  Guide users through classical group sequential trial design using the gsDesign
  R package. Use this skill when the user asks about: group sequential boundaries,
  spending functions (sfLDOF, sfHSD, sfPoints), sample size for time-to-event or
  binomial trials, gsDesign objects, plotting group sequential bounds, gsSurvPower
  for power computation, or harm bounds (test.type 7/8).
---

# Group Sequential Design with gsDesign

**Note**: This skill targets gsDesign >= 3.9.0 (dev version at github.com/keaven/gsDesign).
Features marked **(dev)** may not be on CRAN yet. The CRAN `llms.txt` covers version 3.6.9.

## API reference

- Full function docs (CRAN 3.6.9): `references/llms.txt` (source: https://gsDesign.ai)
- Workflow patterns (dev 3.9.0+): `references/code_patterns.md`

## Key functions

- `gsDesign()` - Core boundary computation for group sequential designs
- `gsSurv()` / `gsSurvCalendar()` / `nSurv()` - Time-to-event design (sample size)
- `gsSurvPower()` - Power computation for survival designs **(dev)**
- `nNormal()` - Sample size for normal endpoints
- `gsBinomialExact()` - Exact binomial group sequential design
- `gsBoundSummary()` - Formatted summary tables of bounds
- `gsProbability()` - Boundary crossing probabilities
- `gsCP()` / `gsBoundCP()` - Conditional power
- `ssrCP()` - Sample size re-estimation based on conditional power
- `plot.gsDesign()` - Plotting group sequential designs
- `hGraph()` - Multiplicity graph visualization
- `toInteger()` - Round sample sizes to integers
- `sequentialPValue()` - Sequential p-value computation

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

## Workflow patterns

For detailed code templates covering common workflows, read `references/code_patterns.md`.

Topics covered:
- One-sided and two-sided designs with various spending functions
- Harm bounds with test.type 7/8 and selective bound testing **(dev)**
- Time-to-event designs with `gsSurv()` (piecewise rates, stratification)
- Calendar-based timing with `gsSurvCalendar()`
- Power computation with `gsSurvPower()` (sensitivity, what-if analyses) **(dev)**
- Normal and binomial endpoint designs
- Integer rounding, bound summaries, and reporting (gt, RTF, LaTeX)
- All 7 plot types
- Spending function selection and comparison
- Sequential p-values for graphical multiplicity
- Conditional power and sample size re-estimation (`ssrCP`)
- Updating bounds when observed timing differs from planned
- Multiplicity graphs with `hGraph()`

## Important design considerations

- **test.type = 4** (non-binding futility with beta-spending) is the most common choice for confirmatory trials
- **test.type = 8** **(dev)** adds non-binding harm bounds for monitoring experimental harm (e.g., OS in oncology)
- **gsSurvPower for what-if** **(dev)**: use `gsSurvPower(x = design, hr = ...)` for sensitivity analyses without re-solving sample size
- **Always round to integers** with `toInteger()` before reporting a design
- **gsSurv over nSurv**: prefer `gsSurv()` which combines sample size and boundary computation
- **Calendar vs event-driven in gsSurvPower**: `plannedCalendarTime` gives unconditional power (events change with HR); `targetEvents` matches the gsDesign power plot
- **Spending time vs information fraction**: use `usTime`/`lsTime` in `gsDesign()` when spending should differ from information fraction
- **sequentialPValue**: only meaningful for `test.type = 1, 4, 6` (one-sided or non-binding futility)
- **Calendar vs information spending**: `gsSurvCalendar(spending = "calendar")` spends less at early interims

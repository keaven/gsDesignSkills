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

**Note**: This skill targets gsDesign >= 3.10.0 (main branch at github.com/keaven/gsDesign).
Features marked **(dev)** are on the main branch but may not be on CRAN yet.

## API reference

- Full function docs (CRAN 3.6.9): `references/llms.txt` (source: https://gsDesign.ai)
- Workflow patterns (dev 3.9.0+): `references/code_patterns.md`

## Key functions

- `gsDesign()` - Core boundary computation for group sequential designs
- `gsSurv()` / `gsSurvCalendar()` / `nSurv()` - Time-to-event design (sample size)
- `gsSurvPower()` - Power computation for survival designs (new in 3.10.0)
- `repeatedPValueBinomialExact()` / `sequentialPValueBinomialExact()` - Exact binomial p-values (new in 3.10.0)
- `simBinomialSeasonalExact()` - Seasonal rare-event simulation (new in 3.10.0)
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

### One-parameter families
- `sfLDOF` - Lan-DeMets O'Brien-Fleming: $2\Phi(-\Phi^{-1}(1-\alpha/2)/\sqrt{t})$
- `sfHSD` - Hwang-Shih-DeCani: $\alpha(1-e^{-\gamma t})/(1-e^{-\gamma})$; $\gamma<0$ conservative, $\gamma>0$ aggressive
- `sfPower` - Kim-DeMets power: $\alpha t^\rho$; $\rho=3$ is conservative (recommended)
- `sfExponential` - Exponential: $\alpha t^{-\nu}$; $\nu=0.8$ approximates O'Brien-Fleming with simpler form; conservative early spending compared to sfHSD or sfPower at same late spending
- `sfPoints` - Pointwise (piecewise linear) spending at specified information fractions
- `sfLinear` - Linear spending: $\alpha t$
- `sfXG` - Xu-Garden conditional error spending

### Two-parameter families (Anderson & Clark, 2009)
All use the general form $\alpha(t; a, b) = \alpha F(a + b F^{-1}(t))$ where $F$ is a CDF. Given two desired spending points $(t_0, s_0)$ and $(t_1, s_1)$, solve $F^{-1}(s_i) = a + b F^{-1}(t_i)$ for the parameters $a, b$. In gsDesign, pass `param = c(t0, t1, s0, s1)` to fit automatically.

- `sfLogistic` - Logistic CDF: $\alpha c(t/(1-t))^b / (1+c(t/(1-t))^b)$ where $c=e^a$. Used in GUSTO V trial and Merck trials. General purpose.
- `sfNormal` - Standard normal CDF. Nearly identical to sfLogistic in practice.
- `sfCauchy` - Cauchy CDF. Flat between fitted points; **robust when analysis timing shifts**.
- `sfExtremeValue` - Extreme value CDF: $\alpha\exp(-e^a(-\ln t)^b)$. Conservative early spending.
- `sfExtremeValue2` - Flipped extreme value: $F(x) = 1-\exp(-\exp(x))$. Conservative early spending; useful for futility bound calibration.
- `sfTDist` - t-distribution CDF (3 parameters: $a$, $b$, and df). Cauchy at df=1, normal at df=$\infty$; df provides continuous interpolation.
- `sfBetaDist` - Incomplete beta CDF with parameters $a>0, b>0$: $\alpha F_{a,b}(t)$. Fitted via `nlminb()` (nonlinear).

**Reference**: Anderson KM, Clark JB. Fitting spending functions. *Statist. Med.* 2010; 29:321–327.

## Output and reporting

- `as_gt()` - Convert summary tables to gt objects
- `as_rtf()` - Save summary tables as RTF files
- `as_table()` - Create summary table objects
- `xtable()` - LaTeX/HTML table output

## Workflow patterns

For detailed code templates covering common workflows, read `references/code_patterns.md`.

Topics covered:
- One-sided and two-sided designs with various spending functions
- Harm bounds with test.type 7/8 and selective bound testing (new in 3.10.0)
- Time-to-event designs with `gsSurv()` (piecewise rates, stratification)
- Calendar-based timing with `gsSurvCalendar()`
- Power computation with `gsSurvPower()` (sensitivity, what-if analyses, `informationRates`, `fullSpendingAtFinal`)
- Normal and binomial endpoint designs
- Integer rounding, bound summaries, and reporting (gt, RTF, LaTeX)
- All 7 plot types
- Spending function selection and comparison
- Two-parameter spending function families (Anderson & Clark 2009): sfLogistic, sfNormal, sfCauchy, sfExtremeValue, sfExtremeValue2, sfTDist, sfBetaDist
- Fitting spending functions to desired boundary values at two information fractions
- Choosing between spending function families
- `sfExtremeValue2` for calibrating futility bounds to target HR
- `testLower`/`testUpper`/`testHarm` for selective bound testing at specific analyses
- Survival designs with `hr > hr0` (reversed HR for time-to-response or safety endpoints)
- `testBinomial()` / `nBinomial()` for binomial Z-statistics and sample size
- Sequential p-values for graphical multiplicity
- Conditional power and sample size re-estimation (`ssrCP`)
- Updating bounds when observed timing differs from planned
- Multiplicity graphs with `hGraph()`

## Important design considerations

- **test.type = 4** (non-binding futility with beta-spending) is the most common choice for confirmatory trials
- **test.type = 7/8** (new in 3.10.0) adds binding (7) or non-binding (8) harm bounds for monitoring experimental harm (e.g., OS in oncology). Controlled by `sfharm`/`sfharmparam`. Per FDA guidance (2025), all randomized oncology trials should include pre-specified OS harm assessment with interim analyses for futility/harm and justified thresholds.
- **gsSurvPower for what-if** (new in 3.10.0): use `gsSurvPower(x = design, hr = ...)` for sensitivity analyses without re-solving sample size. Supports `informationRates` to cap spending and `fullSpendingAtFinal` to force spending fraction to 1 at final analysis
- **Always round to integers** with `toInteger()` before reporting a design
- **toInteger preserves bounds**: as of 3.10.0, `toInteger()` preserves `testLower`, `testUpper`, `testHarm`, and harm-bound spending settings
- **testLower/testUpper/testHarm**: logical vectors controlling which analyses have each boundary type; `FALSE` suppresses that bound at that analysis
- **hr > hr0**: survival functions (`nEvents`, `nSurv`, `gsSurv`, `gsSurvCalendar`) support `hr > hr0` for time-to-response or safety endpoints where a larger HR is the alternative
- **Two-parameter spending families**: `param = c(t1, t2, u1, u2)` fits `sf(t1) = alpha*u1` and `sf(t2) = alpha*u2`. Choose sfCauchy for robustness to timing changes, sfExtremeValue/sfExtremeValue2 for conservative early spending, sfLogistic/sfNormal as general purpose defaults.
- **sfExtremeValue2**: 4-parameter futility spending `c(t1, t2, u1, u2)` for calibrating futility bounds to target a specific HR (e.g., HR = 0.9 to stop for futility)
- **testBinomial sign convention**: `testBinomial(x1=exp, x2=ctrl, n1=n_exp, n2=n_ctrl)` gives positive Z when experimental is better; swapping x1/x2 reverses the sign
- **gsSurv over nSurv**: prefer `gsSurv()` which combines sample size and boundary computation
- **Three ways to derive sample size and obtain desired power in gsSurv/nSurv** (Lachin & Foulkes 1986):
  1. **Vary enrollment rate** (recommended): fix `T` and `minfup`, let `gamma` scale proportionally to achieve target power
  2. **Vary enrollment duration**: fix `gamma` rates, set `T = NULL`; solves for `R` and `T` given `minfup`
  3. **Vary follow-up duration** (not recommended): fix `gamma` and `R`, set `minfup = NULL`; solves for `minfup` and `T`. Fails if enrollment alone over- or under-powers the trial with no/infinite follow-up
- **Calendar vs event-driven in gsSurvPower**: `plannedCalendarTime` gives unconditional power (events change with HR); `targetEvents` matches the gsDesign power plot
- **Spending time vs information fraction**: use `usTime`/`lsTime` in `gsDesign()` when spending should differ from information fraction
- **sequentialPValue**: only meaningful for `test.type = 1, 4, 6` (one-sided or non-binding futility). Based on Liu & Anderson (2008) Theorem 1: Type I error ≤ α for *any* stopping time τ, which is why non-binding futility and trial extension past an efficacy boundary preserve error control.
- **Sequential p-values with multiplicity**: sequential p-values can be passed directly to graphical multiplicity procedures (e.g., `graph_test_shortcut()`) to test multiple hypotheses in group sequential trials while controlling FWER. See Maurer & Bretz (2013) and the graphicalMCP-gsDesign2 skill.
- **Calendar vs information spending**: `gsSurvCalendar(spending = "calendar")` spends less at early interims

## Key references

- Anderson KM, Clark JB. Fitting spending functions. *Statist. Med.* 2010; 29:321–327.
- FDA. Approaches to assessment of overall survival in oncology clinical trials. Draft guidance for industry, August 2025.
- Liu Q, Anderson KM. On adaptive extensions of group sequential trials for clinical investigations. *JASA* 2008; 103:1621–1630.

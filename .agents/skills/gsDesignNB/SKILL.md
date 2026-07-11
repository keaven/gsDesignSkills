---
name: gsDesignNB
description: >
  Guide users through sample size calculation, group sequential design, and
  simulation for clinical trials with negative binomial (recurrent event)
  outcomes using the gsDesignNB R package. Use this skill when the user asks
  about: negative binomial sample size, recurrent event trials, overdispersed
  counts, event gaps, rate ratios, Wald test for count data, seasonal event
  rates, blinded or unblinded sample size re-estimation, group sequential
  designs for negative binomial endpoints, or the Zhu-Lakkis method.
---

# Sample Size and Simulation for Negative Binomial Outcomes with gsDesignNB

## API reference

- Full function docs: `references/llms.txt` (built from local man pages)
- Workflow patterns: `references/code_patterns.md`

## Key functions

### Sample size calculation
- `sample_size_nbinom()` - Sample size/power for negative binomial outcomes (Zhu-Lakkis Method 3)

### Group sequential design
- `gsNBCalendar()` - Group sequential design with calendar-time analysis schedule
- `compute_info_at_time()` - Statistical information at a given calendar time
- `toInteger()` - Round sample sizes to integers preserving allocation ratio
- `check_gs_bound()` - Check if group sequential bounds are crossed
- `summarize_gs_sim()` - Summarize operating characteristics from simulations

### Simulation
- `nb_sim()` - Simulate recurrent events (Gamma-Poisson mixture)
- `nb_sim_seasonal()` - Simulate recurrent events with seasonal variation
- `sim_gs_nbinom()` - Simulate multiple group sequential trials

### Data cutting and analysis timing
- `cut_data_by_date()` - Cut simulated data at a calendar date
- `get_analysis_date()` - Find date when target event count is reached
- `get_cut_date()` - Find earliest date satisfying multiple analysis criteria
- `cut_date_for_completers()` - Find date when target completers are reached
- `cut_completers()` - Cut data for completers analysis

### Statistical testing and estimation
- `mutze_test()` - Wald test for treatment rate ratio (NB or Poisson)
- `estimate_nb_mom()` - Method of moments estimation for NB parameters
- `calculate_blinded_info()` - Blinded information and dispersion estimation

### Sample size re-estimation
- `blinded_ssr()` - Blinded SSR using Friede & Schmidli method
- `unblinded_ssr()` - Unblinded SSR using observed group rates

## Workflow patterns

For detailed code templates, read `references/code_patterns.md`.

Topics covered:
- Fixed sample size calculation with piecewise accrual, dropout, event gaps
- Power calculation from a fixed design
- Non-inferiority and super-superiority designs (rr0 parameter)
- Group sequential design with calendar-time analysis schedule
- Simulation of recurrent events and group sequential trials
- Seasonal event rate simulation
- Data cutting at interim/final analyses
- Wald test (mutze_test) for treatment rate ratio
- Blinded and unblinded information estimation at interim
- Blinded and unblinded sample size re-estimation
- Completers-based interim analysis
- Verification of theoretical vs. simulated operating characteristics

## Important design considerations

- **Dispersion parameter k**: Controls overdispersion; k = 0 reduces to Poisson. Larger k means more overdispersion. Can be scalar (common) or length-2 vector (group-specific).
- **Event gaps**: After each event, patients are "off risk" for `event_gap` time units. This reduces effective exposure: `lambda_eff = lambda / (1 + lambda * gap)`. Specified in the same time units as rates.
- **Variance inflation factor Q**: When follow-up varies across patients, `Q = E[t^2] / E[t]^2` inflates the variance. `sample_size_nbinom()` handles this automatically.
- **Rate parameterization**: Rates lambda1 (control) and lambda2 (experimental) are events per unit time. The treatment effect is the rate ratio RR = lambda2/lambda1.
- **Wald test (Mütze et al.)**: `mutze_test()` fits a negative binomial GLM with offset for log exposure. Falls back to Poisson when the NB dispersion estimate is very large (> poisson_threshold).
- **Calendar-time analysis schedule**: `gsNBCalendar()` takes `analysis_times` as calendar months. Information at each analysis depends on enrollment pattern, dropout, and follow-up.
- **Spending time vs. information fraction**: For group sequential designs, `usTime`/`lsTime` control alpha spending and may differ from the information fraction. This allows calendar-based or event-based spending schedules.
- **Blinded information**: `calculate_blinded_info()` uses the blinded (pooled) rate and dispersion to estimate information. Can produce extreme values when the NB MLE is unstable — bound dispersion or use planning values as a fallback.
- **SSR**: Blinded SSR (Friede & Schmidli) maintains the blind; unblinded SSR uses observed group rates for more accurate re-estimation but requires unblinding.
- **`check_gs_bound()` info_scale**: Use "blinded" (default) or "unblinded" to select which information drives bound updates at analysis time.

# Replicate dynamic educational returns model with R

This repository uses R and replicates "Returns to Education: The Causal Effects of
Education on Earnings, Health, and Smoking" (2018) by J. Heckman, J. Humphries, and G. Veramendi especially for *returns to education*.
---

## Repository Structure
<pre>
|-- code # R scripts to perform analysis with dynamic returns model for educational returns
|-- reference # includes original paper and web appendix provided online.
</pre>

## Instruction
This section will walk through all replicatio process with R from data retrieval, variable creation, and analysis.

### Returns-to-education R replication

Run the R post-estimation script from the repository root:

```sh
Rscript code/replicate_returns_to_education.R
```

The script ports the Stata post-estimation logic for the paper's returns-to-education outcomes, focusing on `wage` and `PVwage`. If the original simulated potential-outcome file is available as `reference/2015166data/replication/ana/simulation_data.dta` or `data/simulation_data.dta`, it writes:

- `output/dynamic_returns_te_summary.csv`
- `output/dynamic_returns_option_value_summary.csv`
- `output/dynamic_returns_prte_summary.csv`
- `output/dynamic_returns_te_by_margin.png`

The current reference folder does not include `simulation_data.dta`, which the original replication creates with the C++/Stata simulation step. In that case, the script still writes `output/observed_returns_to_education_fallback.csv` as a descriptive observed-data check, but the dynamic model tables require the simulated file.

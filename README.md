# Uncertainty Analysis

Uncertainty and sensitivity analysis of crop model outputs, including local/global sensitivity, variance decomposition, and CSV-driven design points integrated with model templates.


## Repository structure:

<!--not set in stone!-->

```
├── README.md
├── 000-config.yml
├── R
│   ├── global_sensitivity.R
│   ├── local_sensitivity.R
│   └── variance_decomposition.R
├── analysis/
│   ├── global_sensitivity.qmd
│   ├── local_sensitivity.qmd
│   └── variance_decomposition.qmd
├── data_raw/   
│   ├── sa_design_points.csv
│   └── template.xml
├── scripts/
│   ├── 001_setup_design_points.R
│   ├── 011_run_local_sensitivity.R
│   ├── 012_aggregate_sensitivity.R
│   ├── 021_generate_sobol_design.R
│   ├── 022_run_global_sensitivity.R
│   ├── 023_compute_sobol_indices.R
│   ├── 031_partition_variance.R
│   └── 032_hierarchical_variance.R
├── docs/
├── tests/
└── reports/
    └── uncertainty_analysis.qmd
```

note: `data_raw` is for data of limited size (<MB) that is input to the pipeline; small outputs from these workflows can go in 'data/' but most inputs and outputs will go in one of the outdirs listed in config.yml

## Configuration

Trying something new:
- putting configuration in `000-config.yml` and reading with `config::get(file = "000-config.yml")`.
- Added PEcAn settings template (`template.xml`) is a placeholder from the workflows repository; needs sensitivity blocks added. config.yml should not duplicate content of the pecan.xml 

## Variance Decomposition Workflows

The supported Sobol workflow for this repository is the PEcAn-backed path in
`scripts/033_partition_variance_pecan_sobol.R`. In this arrangement, PEcAn is
the Sobol engine and `uncertainty` is the consumer/reporting layer:

1. PEcAn generates the Sobol design and model runs.
2. PEcAn saves `sobol.indices.*.Rdata` and `sobol.design.*.Rdata`.
3. `uncertainty` reads those saved indices, partitions variance across type-level factors such as `param`, `met`, `poolinitcond`, and `events`, then applies the existing local parameter drilldown only within the `param` contribution.

Consumer contract for the PEcAn-backed path:

- `sobol.indices.*.Rdata` is the primary source of Sobol numeric results.
- Factor/type metadata comes from `sobol_results$source_type` and `source_tag`; `sobol.design.*.Rdata` is fallback-only when those columns are absent.
- `runid`, `variable`, `start_year`, and `end_year` come from canonical PEcAn filenames, which are produced from PEcAn's analysis filename contract.
- `pecan.CONFIGS.xml` is used only for ensemble variance calculation and best-effort run-to-site mapping.

Required upstream inputs for the PEcAn-backed path:

- `pecan.CONFIGS.xml`
- `ensemble.output.*.Rdata`
- `sobol.indices.*.Rdata`
- `data/aggregated_sensitivity.csv`

The older custom Sobol path in `scripts/031_partition_variance.R` remains only
for legacy PR #2 context. It is not the supported Sobol engine for new work.

Report entry point for the PEcAn-backed path: `analysis/variance_decomposition_pecan_sobol.qmd`.

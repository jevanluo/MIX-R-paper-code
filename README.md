# MIXTree Empirical Analysis Code

This folder contains the minimal code needed to rerun the empirical Verbal Aggression MIX-R analysis.

## Contents

- `run_empirical_analysis.R`: reproducible entry point for the empirical MIX-R fit.
- `code/`: model and helper functions used by the empirical entry point.
- `data/verbal_agg.rdata`: empirical response matrix used by the analysis.

Generated posterior fit files and tables are not included. Running the analysis writes outputs to `results/`, which is ignored by git.

## Model Fitted

The wrapper fits only MIX-R, implemented as `fit_gate2_2pl()`, to all response columns in `data/verbal_agg.rdata`. It does not split the responses into `want` and `do` subsets.

## Run

From this folder:

```sh
Rscript run_empirical_analysis.R
```

For a short smoke-test run:

```sh
Rscript run_empirical_analysis.R --niter=20 --nburn=10 --thin=1 --nchains=1
```

The full empirical run is computationally heavy. The default MCMC settings reproduce the source script settings: `niter=100000`, `nburn=20000`, `thin=10`, and `nchains=4`.


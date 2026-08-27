#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) {
  gsub("~+~", " ", sub("^--file=", "", file_arg[[1]]), fixed = TRUE)
} else {
  "code/empirical/run_empirical_analysis.R"
}
script_path <- normalizePath(script_path, mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

trailing <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^", name, "="), trailing, value = TRUE)
  if (length(hit)) sub(paste0("^", name, "="), "", hit[[1]]) else default
}
as_flag <- function(x) tolower(x) %in% c("1", "true", "t", "yes", "y")

niter <- as.integer(get_arg("--niter", "24000"))
nburn <- as.integer(get_arg("--nburn", "4000"))
thin <- as.integer(get_arg("--thin", "10"))
nchains <- as.integer(get_arg("--nchains", "4"))
cores <- as.integer(get_arg("--cores", as.character(nchains)))
seed <- as.integer(get_arg("--seed", "1"))
auto_mm <- as_flag(get_arg("--auto-mm", "true"))
output_arg <- get_arg(
  "--output-dir",
  file.path(root, "results", "empirical")
)
if (!grepl("^/", output_arg)) output_arg <- file.path(root, output_arg)
results_dir <- normalizePath(output_arg, mustWork = FALSE)

if (anyNA(c(niter, nburn, thin, nchains, cores, seed))) {
  stop("MCMC arguments must be integers.", call. = FALSE)
}
if (niter <= 0 || nburn < 0 || thin <= 0 || nchains <= 0 || cores <= 0) {
  stop("Invalid MCMC settings.", call. = FALSE)
}

source(file.path(root, "code", "empirical", "prepare_verbal_aggression_data.R"))
Y <- load_publication_verbal_aggression()
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

simulation_dir <- file.path(root, "code", "simulations")
old_wd <- setwd(simulation_dir)
on.exit(setwd(old_wd), add = TRUE)

library(nimble)
library(loo)
source("Codes/functions_utilities.R")
source("Codes/functions_experiment_2pl_parallel.R")
source("Codes/functions_prevalence.R")

set.seed(seed)
message("Fitting MIX-R to all Verbal Aggression responses with fit_gate2_multicore_withIC().")
mix_R <- fit_gate2_multicore_withIC(
  Y,
  niter = niter,
  nburn = nburn,
  thin = thin,
  nchains = nchains,
  cores = cores,
  auto_mm = auto_mm
)

saveRDS(mix_R, file = file.path(results_dir, "empirical_all_mix_R_parallel_fit.rds"))
saveRDS(
  list(
    model = "MIX-R",
    data_source = "dexter::verbAggrData[, 3:26] + 1",
    item_names = colnames(Y),
    dimensions = dim(Y),
    niter = niter,
    nburn = nburn,
    thin = thin,
    nchains = nchains,
    cores = cores,
    seed = seed,
    auto_mm = auto_mm,
    package_versions = vapply(
      c("R", "nimble", "loo", "coda", "dexter", "posterior"),
      function(pkg) {
        if (pkg == "R") return(as.character(getRversion()))
        if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
        as.character(utils::packageVersion(pkg))
      },
      character(1)
    )
  ),
  file = file.path(results_dir, "empirical_all_mix_R_parallel_run_metadata.rds")
)
writeLines(
  capture.output(sessionInfo()),
  con = file.path(results_dir, "empirical_all_mix_R_parallel_session_info.txt")
)

gof <- tryCatch(report_gof(mix_R), error = function(e) {
  warning("Could not compute GOF table: ", conditionMessage(e), call. = FALSE)
  NULL
})
if (!is.null(gof)) {
  write.csv(
    data.frame(model = "MIX-R", t(gof), row.names = NULL, check.names = FALSE),
    file = file.path(results_dir, "empirical_all_mix_R_parallel_gof.csv"),
    row.names = FALSE
  )
}

prev <- tryCatch(prevalence_from_gate_draws(mix_R$samples, N = nrow(Y), I = ncol(Y)),
                 error = function(e) NULL)
if (!is.null(prev)) {
  saveRDS(prev, file = file.path(results_dir, "empirical_all_mix_R_parallel_prevalence.rds"))
}

message("Done. Results written to ", results_dir)

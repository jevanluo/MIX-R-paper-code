#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "run_empirical_analysis.R"
root <- normalizePath(dirname(script_path), mustWork = TRUE)

trailing <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^", name, "="), trailing, value = TRUE)
  if (length(hit)) sub(paste0("^", name, "="), "", hit[[1]]) else default
}

niter <- as.integer(get_arg("--niter", "100000"))
nburn <- as.integer(get_arg("--nburn", "20000"))
thin <- as.integer(get_arg("--thin", "10"))
nchains <- as.integer(get_arg("--nchains", "4"))
seed <- as.integer(get_arg("--seed", "1"))

source(file.path(root, "code", "functions_utilities.R"))
source(file.path(root, "code", "functions_experiment_2pl.R"))

load(file.path(root, "data", "verbal_agg.rdata"))
if (!exists("y")) {
  stop("Expected object `y` in data/verbal_agg.rdata.", call. = FALSE)
}

dir.create(file.path(root, "results"), showWarnings = FALSE)
set.seed(seed)

message("Fitting MIX-R to all empirical responses.")
mix_R <- fit_gate2_2pl(y, niter = niter, nburn = nburn, thin = thin, nchains = nchains)
saveRDS(mix_R, file = file.path(root, "results", "empirical_all_mix_R_fit.rds"))

gof <- tryCatch(report_gof(mix_R), error = function(e) NULL)
if (!is.null(gof)) {
  write.csv(
    data.frame(model = "MIX-R", t(gof), row.names = NULL, check.names = FALSE),
    file = file.path(root, "results", "empirical_all_mix_R_gof.csv"),
    row.names = FALSE
  )
}

message("Done.")

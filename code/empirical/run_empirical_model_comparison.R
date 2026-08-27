#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) {
  gsub("~+~", " ", sub("^--file=", "", file_arg[[1]]), fixed = TRUE)
} else {
  "code/empirical/run_empirical_model_comparison.R"
}
script_path <- normalizePath(script_path, mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

trailing <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^", name, "="), trailing, value = TRUE)
  if (length(hit)) sub(paste0("^", name, "="), "", hit[[1]]) else default
}

niter <- as.integer(get_arg("--niter", "24000"))
nburn <- as.integer(get_arg("--nburn", "4000"))
thin <- as.integer(get_arg("--thin", "10"))
nchains <- as.integer(get_arg("--nchains", "4"))
seed <- as.integer(get_arg("--seed", "1"))
models <- strsplit(tolower(get_arg("--models", "lpt,npt,mixp")), ",", fixed = TRUE)[[1]]
models <- trimws(models)
output_arg <- get_arg(
  "--output-dir",
  file.path(root, "results", "empirical", "model_comparison")
)
if (!grepl("^/", output_arg)) output_arg <- file.path(root, output_arg)
output_dir <- normalizePath(output_arg, mustWork = FALSE)

valid_models <- c("lpt", "npt", "mixp")
if (!length(models) || any(!models %in% valid_models)) {
  stop("--models must contain one or more of: ", paste(valid_models, collapse = ", "), call. = FALSE)
}
if (anyNA(c(niter, nburn, thin, nchains, seed))) {
  stop("MCMC arguments must be integers.", call. = FALSE)
}
if (niter <= 0 || nburn < 0 || thin <= 0 || nchains <= 0) {
  stop("Invalid MCMC settings.", call. = FALSE)
}

source(file.path(root, "code", "empirical", "prepare_verbal_aggression_data.R"))
Y <- load_publication_verbal_aggression()
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

simulation_dir <- file.path(root, "code", "simulations")
old_wd <- setwd(simulation_dir)
on.exit(setwd(old_wd), add = TRUE)

suppressPackageStartupMessages({
  library(nimble)
  library(loo)
})
source("Codes/functions_utilities.R")
source("Codes/functions_baseline_3cat_tree1.R")
source("Codes/functions_baseline_3cat_tree2.R")
source("Codes/functions_baseline_3cat_between_person.R")

fit_one <- function(model) {
  set.seed(seed)
  message("Fitting ", toupper(model), " with ", nchains, " chain(s).")
  fit <- switch(
    model,
    lpt = fit_BaselineA(Y, niter = niter, nburn = nburn, thin = thin, nchains = nchains),
    npt = fit_BaselineB(Y, niter = niter, nburn = nburn, thin = thin, nchains = nchains),
    mixp = fit_Gate6_personSwitch(Y, niter = niter, nburn = nburn, thin = thin, nchains = nchains)
  )
  saveRDS(fit, file.path(output_dir, paste0("empirical_", model, "_fit.rds")))
  fit
}

fits <- setNames(lapply(models, fit_one), models)

gof_rows <- lapply(names(fits), function(model) {
  gof <- report_gof(fits[[model]])
  data.frame(
    model = toupper(model),
    WAIC = unname(gof[["WAIC"]]),
    pWAIC = unname(gof[["pWAIC"]]),
    Pareto_gt_0_7 = unname(gof[["Pareto_gt_0_7"]]),
    LOOIC = unname(gof[["LOOIC"]]),
    pLOO = unname(gof[["pLOO"]])
  )
})
gof_table <- do.call(rbind, gof_rows)
write.csv(gof_table, file.path(output_dir, "empirical_model_comparison_gof.csv"), row.names = FALSE)
saveRDS(
  list(
    models = models,
    data_source = "dexter::verbAggrData[, 3:26] + 1",
    item_names = colnames(Y),
    dimensions = dim(Y),
    niter = niter,
    nburn = nburn,
    thin = thin,
    nchains = nchains,
    seed = seed,
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
  file = file.path(output_dir, "empirical_model_comparison_run_metadata.rds")
)
writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "empirical_model_comparison_session_info.txt")
)

message("Done. Results written to ", output_dir)

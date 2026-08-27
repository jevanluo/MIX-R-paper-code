# Construct the empirical response matrix exactly as in Application-III/ea_analysis.R.
#
# The historical data/verbal_agg.rdata file contains the same response profiles,
# but its rows and columns are ordered differently.  That difference leaves the
# aggregate likelihood unchanged while changing the correspondence between
# posterior indices and the published person/item labels.  The dexter data set
# is therefore the canonical source for publication reproduction.

load_publication_verbal_aggression <- function() {
  if (!requireNamespace("dexter", quietly = TRUE)) {
    stop(
      "Package `dexter` is required to load the publication data. ",
      "Install the version recorded in the run metadata and try again.",
      call. = FALSE
    )
  }

  data_env <- new.env(parent = baseenv())
  utils::data("verbAggrData", package = "dexter", envir = data_env)
  if (!exists("verbAggrData", envir = data_env, inherits = FALSE)) {
    stop("Could not load dexter::verbAggrData.", call. = FALSE)
  }

  item_names <- c(
    "S1DoCurse", "S1DoScold", "S1DoShout",
    "S1WantCurse", "S1WantScold", "S1WantShout",
    "S2DoCurse", "S2DoScold", "S2DoShout",
    "S2WantCurse", "S2WantScold", "S2WantShout",
    "S3DoCurse", "S3DoScold", "S3DoShout",
    "S3WantCurse", "S3WantScold", "S3WantShout",
    "S4DoCurse", "S4DoScold", "S4DoShout",
    "S4WantCurse", "S4WantScold", "S4WantShout"
  )
  source_data <- data_env$verbAggrData
  if (!identical(names(source_data)[3:26], item_names)) {
    stop(
      "dexter::verbAggrData does not have the item order used by the publication. ",
      "Check the installed dexter version before fitting.",
      call. = FALSE
    )
  }

  Y <- data.matrix(source_data[, 3:26, drop = FALSE]) + 1L
  storage.mode(Y) <- "integer"

  if (!identical(dim(Y), c(316L, 24L))) {
    stop("Expected a 316 x 24 publication response matrix.", call. = FALSE)
  }
  if (anyNA(Y) || !identical(sort(unique(as.vector(Y))), 1:3)) {
    stop("Publication responses must be complete and coded 1, 2, or 3.", call. = FALSE)
  }

  Y
}

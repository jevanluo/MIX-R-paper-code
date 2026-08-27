#!/usr/bin/env Rscript

# Regenerate the two empirical-study PDFs from the saved R analysis objects.
# The plotting code follows response_revise_emprical.R, ea_analysis.R, and
# ea_functions.R. The dated PDFs match the independent no-covariate MIX-R fit
# in res_all_items_26apr21.rds[[1]]. Only the color scales are changed.

suppressPackageStartupMessages({
  library(dexter)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(tibble)
})

script_path <- function() {
  command <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command, value = TRUE)
  if (length(file_arg)) {
    decoded <- gsub("~+~", " ", sub("^--file=", "", file_arg[[1]]), fixed = TRUE)
    return(normalizePath(decoded))
  }
  normalizePath(getwd())
}

arg_value <- function(name, default) {
  args <- commandArgs(trailingOnly = TRUE)
  prefix <- paste0("--", name, "=")
  match <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(match)) sub(paste0("^", prefix), "", match[[1]]) else default
}

require_file <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
  normalizePath(path, mustWork = TRUE)
}

script <- script_path()
repo_root <- normalizePath(file.path(dirname(script), "..", ".."))
fit_file <- path.expand(arg_value(
  "fit-file",
  file.path(repo_root, "results", "empirical", "empirical_all_mix_R_parallel_fit.rds")
))
output_dir <- path.expand(arg_value(
  "output-dir",
  file.path(repo_root, "empirical study")
))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data_env <- new.env(parent = baseenv())
utils::data("verbAggrProperties", package = "dexter", envir = data_env)
if (!exists("verbAggrProperties", envir = data_env, inherits = FALSE)) {
  stop("Could not load dexter::verbAggrProperties.", call. = FALSE)
}
item_info <- data_env$verbAggrProperties

##################################################
## Item gate parameter figure
##################################################

fit_object <- readRDS(require_file(fit_file))
if (is.list(fit_object) &&
    !is.null(fit_object$param_summaries) &&
    !is.null(fit_object$samples)) {
  no_covariate_fit <- fit_object
} else if (is.list(fit_object) && length(fit_object) >= 1) {
  no_covariate_fit <- fit_object[[1]]
} else {
  stop("The fit file does not contain a recognizable MIX-R fit.", call. = FALSE)
}
if (is.null(no_covariate_fit$param_summaries$beta_gate) ||
    is.null(no_covariate_fit$samples)) {
  stop(
    "The selected no-covariate fit lacks beta_gate summaries or samples.",
    call. = FALSE
  )
}
beta_gate <- no_covariate_fit$param_summaries$beta_gate
samples <- no_covariate_fit$samples
rm(fit_object, no_covariate_fit)
invisible(gc())

if (nrow(beta_gate) != nrow(item_info)) {
  stop("beta_gate and verbAggrProperties have different row counts.", call. = FALSE)
}

plot_df <- cbind(beta_gate, item_info) |>
  mutate(
    item_label = item_id,
    behavior = factor(behavior, levels = c("Curse", "Scold", "Shout")),
    mode = factor(mode, levels = c("Do", "Want")),
    blame = factor(blame, levels = c("Other", "Self"))
  )

shape_map <- c(Curse = 21, Scold = 22, Shout = 24)
gray_mode_fills <- c(Do = "grey25", Want = "grey75")

item_gate_plot <- ggplot(plot_df, aes(x = item_label, y = mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(
    aes(ymin = X2.5., ymax = X97.5.),
    width = 0.15,
    linewidth = 0.5
  ) +
  geom_point(
    aes(shape = behavior, fill = mode),
    color = "grey60",
    size = 6,
    stroke = 0.5,
    show.legend = TRUE
  ) +
  geom_point(
    aes(shape = behavior, color = blame),
    fill = NA,
    size = 6,
    stroke = 1.2,
    show.legend = TRUE
  ) +
  scale_shape_manual(values = shape_map, name = "Behavior") +
  scale_fill_manual(values = gray_mode_fills, name = "Mode") +
  scale_color_manual(
    values = c(Other = "transparent", Self = "black"),
    name = "Situation type",
    labels = c(Other = "Other-to-blame", Self = "Self-to-blame")
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21, color = "white"),
      order = 2
    ),
    shape = guide_legend(order = 1),
    color = guide_legend(
      override.aes = list(shape = 21, fill = NA),
      order = 3
    )
  ) +
  labs(x = NULL, y = expression(hat(v)[i]), title = NULL) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

item_output <- file.path(output_dir, "item_gate_26may08_grayscale.pdf")
pdf(file = item_output, height = 6, width = 10)
print(item_gate_plot)
invisible(dev.off())
message("Wrote ", normalizePath(item_output))

##################################################
## Response-level prevalence figure
##################################################

combine_draws <- function(samples) {
  if (is.list(samples) &&
      all(vapply(samples, function(x) is.matrix(x) || is.data.frame(x), logical(1)))) {
    return(as.data.frame(do.call(rbind, lapply(samples, as.matrix))))
  }
  as.data.frame(samples)
}

indexed_parameter_names <- function(names, pattern, label) {
  matched <- grep(pattern, names, value = TRUE)
  if (!length(matched)) stop("No ", label, " columns found.", call. = FALSE)
  indices <- suppressWarnings(as.integer(sub(pattern, "\\1", matched)))
  matched[order(indices, na.last = TRUE)]
}

response_prevalence <- function(samples) {
  draws <- combine_draws(samples)
  if (!"alpha0" %in% names(draws)) {
    stop("Posterior samples do not contain alpha0.", call. = FALSE)
  }

  alpha_names <- indexed_parameter_names(
    names(draws),
    "^alpha_gate\\[(\\d+)\\]$",
    "alpha_gate"
  )
  beta_names <- indexed_parameter_names(
    names(draws),
    "^beta_gate\\[(\\d+)\\]$",
    "beta_gate"
  )

  alpha0 <- as.numeric(draws$alpha0)
  alpha_gate <- t(as.matrix(draws[alpha_names]))
  beta_gate <- t(as.matrix(draws[beta_names]))
  draw_count <- length(alpha0)

  p_response_mean <- matrix(
    0,
    nrow = nrow(alpha_gate),
    ncol = nrow(beta_gate)
  )
  for (draw in seq_len(draw_count)) {
    eta <- outer(alpha_gate[, draw], beta_gate[, draw], "+") + alpha0[draw]
    p_response_mean <- p_response_mean + plogis(eta)
  }
  p_response_mean <- p_response_mean / draw_count

  eta_plugin <- outer(
    rowMeans(alpha_gate),
    rowMeans(beta_gate),
    "+"
  ) + mean(alpha0)
  person_prevalence <- rowMeans(plogis(eta_plugin))

  list(
    p_response_mean = p_response_mean,
    person_prevalence = person_prevalence
  )
}

tidy_response <- function(pmat, col_labels = NULL) {
  stopifnot(is.matrix(pmat))
  if (is.null(rownames(pmat))) rownames(pmat) <- paste0("row_", seq_len(nrow(pmat)))
  if (!is.null(col_labels)) {
    stopifnot(length(col_labels) == ncol(pmat))
    colnames(pmat) <- col_labels
  } else if (is.null(colnames(pmat))) {
    colnames(pmat) <- paste0("Item", seq_len(ncol(pmat)))
  }

  as.data.frame(pmat) |>
    rownames_to_column("person") |>
    pivot_longer(-person, names_to = "item", values_to = "score") |>
    mutate(item = factor(item, levels = colnames(pmat)))
}

gray_prevalence_scale <- function() {
  scale_fill_gradient(
    low = "grey15",
    high = "grey85",
    limits = c(0, 1),
    name = "p"
  )
}

plot_response_heatmap <- function(
    pmat,
    col_labels = NULL,
    y_show = c("all", "auto", "none"),
    y_every = NULL,
    y_as_index = FALSE,
    legend_position = "right",
    outer_lr_margin = 50,
    inner_lr_pad = 0.30) {
  y_show <- match.arg(y_show)
  df <- tidy_response(pmat, col_labels)
  df$person_lab <- as.character(df$person)
  if (y_as_index) df$person_lab <- sub(".*\\[(\\d+)\\]", "\\1", df$person_lab)
  df$person_lab <- factor(df$person_lab, levels = unique(df$person_lab))

  plot <- ggplot(df, aes(x = item, y = person_lab)) +
    geom_tile(aes(fill = score), color = "grey90", linewidth = 0.2) +
    gray_prevalence_scale() +
    scale_x_discrete(expand = expansion(mult = 0, add = inner_lr_pad))

  levels <- levels(df$person_lab)
  if (y_show == "none") {
    plot <- plot +
      scale_y_discrete(breaks = NULL, labels = NULL) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  } else {
    if (y_show == "auto" && is.null(y_every)) {
      y_every <- max(1, ceiling(length(levels) / 15))
    }
    if (!is.null(y_every) && y_every > 1) {
      plot <- plot + scale_y_discrete(
        breaks = levels[seq(1, length(levels), by = y_every)]
      )
    }
  }

  plot +
    labs(x = NULL, y = NULL, title = NULL) +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = legend_position,
      plot.margin = grid::unit(
        c(4.5, outer_lr_margin, 4.5, outer_lr_margin),
        "pt"
      )
    )
}

plot_response_points <- function(
    pmat,
    row_sel,
    col_labels = NULL,
    point_size = 8) {
  df <- tidy_response(pmat, col_labels)
  row_names <- rownames(pmat)
  if (is.null(row_names)) row_names <- paste0("row_", seq_len(nrow(pmat)))
  rows <- if (is.numeric(row_sel)) row_names[row_sel] else as.character(row_sel)

  selected <- df |>
    filter(person %in% rows) |>
    mutate(person = factor(person, levels = rows))

  indices <- suppressWarnings(as.integer(sub(".*\\[(\\d+)\\].*", "\\1", rows)))
  indices[is.na(indices)] <- match(rows[is.na(indices)], row_names)
  label_map <- setNames(paste0("person[", indices, "]"), rows)

  ggplot(selected, aes(x = item, y = person)) +
    geom_point(
      aes(fill = score),
      shape = 21,
      color = "grey20",
      size = point_size,
      stroke = 0.3
    ) +
    gray_prevalence_scale() +
    scale_y_discrete(breaks = rows, labels = label_map) +
    labs(x = NULL, y = NULL, title = NULL) +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

prevalence <- response_prevalence(samples)
rm(samples)
invisible(gc())

pmat <- prevalence$p_response_mean
person_prevalence <- prevalence$person_prevalence
item_labels <- item_info$item_id
pmat_reordered <- pmat[order(person_prevalence, decreasing = TRUE), , drop = FALSE]
extreme_people <- c(
  which(person_prevalence < 0.29),
  which(person_prevalence > 0.664)
)

response_output <- file.path(
  output_dir,
  "response.prevlance_26may08_grayscale.pdf"
)
pdf(file = response_output, height = 8, width = 12)
print(plot_response_heatmap(pmat, item_labels, y_show = "none"))
print(plot_response_heatmap(pmat_reordered, item_labels, y_show = "none"))
print(plot_response_points(pmat, seq_len(20), item_labels))
print(plot_response_points(pmat, extreme_people, item_labels, point_size = 12))
invisible(dev.off())
message("Wrote ", normalizePath(response_output))

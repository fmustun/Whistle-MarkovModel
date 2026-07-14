#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo_root <- if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."),
                mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(repo_root, "R", "statistical_utils.R"))

observed <- matrix(c(0.5, 0.2, 0, 0.3), nrow = 2)
nulls <- list(
  matrix(c(0.5, 0.1, 0.1, 0.4), nrow = 2),
  matrix(c(0.6, 0.2, 0, 0.2), nrow = 2)
)
result <- empirical_enrichment_p_values(
  observed, nulls, transition_probability_function = identity,
  progress_every = 0L
)

stopifnot(
  identical(dim(result$p_value), c(2L, 2L)),
  identical(rownames(result$p_value), c("1", "2")),
  result$p_value[1, 2] == 1,
  result$p_value[2, 1] == 2 / 3,
  result$p_value[2, 2] == 2 / 3,
  all(result$p_value > 0)
)

bh <- adjust_transition_p_values(result$p_value, "BH")
uncorrected <- adjust_transition_p_values(result$p_value, "none")
stopifnot(
  identical(dim(bh), c(2L, 2L)),
  identical(dimnames(bh), dimnames(result$p_value)),
  identical(uncorrected, result$p_value),
  all(bh >= result$p_value)
)

selection <- select_transition_matrix(observed, uncorrected, alpha = 0.9)
counts <- transition_selection_counts(selection$mask)
stopifnot(
  counts[["retained_nodes"]] == 2,
  counts[["significant_transitions_including_loops"]] == 2,
  counts[["significant_auto_loops"]] == 1,
  counts[["loop_free_edges"]] == 1
)

family_mask <- matrix(c(TRUE, FALSE, FALSE, TRUE), nrow = 2)
raw_p <- matrix(c(0.01, 1, 1, 0.02), nrow = 2)
masked_bh <- adjust_transition_p_values(raw_p, "BH", family_mask = family_mask)
stopifnot(
  isTRUE(all.equal(
    unname(masked_bh[family_mask]),
    stats::p.adjust(raw_p[family_mask], method = "BH")
  )),
  all(masked_bh[!family_mask] == 1)
)

stopifnot(
  recommended_iterations(m = 1214, alpha = 0.05, margin = 5) == 121400L,
  recommended_iterations(m = 100, alpha = 0.05, margin = 1) ==
    as.integer(ceiling(100 / 0.05))
)

message("All statistical utility tests passed")

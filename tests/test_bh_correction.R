#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
REPO_ROOT <- if (length(file_arg)) {
  normalizePath(
    file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."),
    mustWork = TRUE
  )
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
source(file.path(REPO_ROOT, "R", "statistical_utils.R"))

observed <- matrix(c(0.5, 0.5, 0, 1), nrow = 2L, byrow = TRUE)
null_probabilities <- list(
  matrix(c(0.5, 0.5, 0, 1), nrow = 2L, byrow = TRUE),
  matrix(c(0.6, 0.4, 0, 1), nrow = 2L, byrow = TRUE),
  matrix(c(0.4, 0.6, 0, 1), nrow = 2L, byrow = TRUE)
)
p_result <- empirical_enrichment_p_values(
  observed,
  null_probabilities,
  transition_probability_function = identity,
  progress_every = 0L
)
stopifnot(
  identical(unname(p_result$exceedance_count[1L, 1L]), 2L),
  identical(unname(p_result$exceedance_count[1L, 2L]), 2L),
  identical(unname(p_result$p_value[1L, 1L]), 0.75),
  identical(unname(p_result$p_value[1L, 2L]), 0.75),
  identical(unname(p_result$p_value[2L, 1L]), 1),
  identical(unname(p_result$p_value[2L, 2L]), 1)
)

raw_p <- matrix(c(0.001, 0.02, 1, 0.04), nrow = 2L)
bh_q <- adjust_transition_p_values(raw_p, correction = "BH")
stopifnot(isTRUE(all.equal(
  as.vector(bh_q), stats::p.adjust(as.vector(raw_p), method = "BH")
)))

observed_weights <- matrix(c(0.4, 0.3, 0.2, 0.1), nrow = 2L)
selection <- select_transition_matrix(observed_weights, bh_q, alpha = 0.05)
stopifnot(
  identical(selection$matrix[selection$mask], observed_weights[selection$mask]),
  all(selection$matrix[!selection$mask] == 0),
  identical(
    unname(transition_selection_counts(selection$mask)),
    c(2L, 2L, 1L, 1L)
  )
)

message("BH correction tests passed")

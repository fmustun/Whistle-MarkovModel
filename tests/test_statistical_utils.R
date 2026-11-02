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

# observed strength (diagonal stripped): node1 = 5+2=7, node2 = 4+2=6,
# node3 = 2+7=9. Both null iterations strip to strength 4 and 0 respectively
# for every node, so null_mean = 2 and null_sd = sqrt(8) for all three.
observed_counts <- matrix(c(1, 2, 3, 0, 1, 4, 2, 0, 1), nrow = 3, byrow = TRUE)
null_counts <- list(
  matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3, byrow = TRUE),
  matrix(0, nrow = 3, ncol = 3)
)
strength_null <- strength_null_distribution(
  observed_counts, null_counts, progress_every = 0L
)
stopifnot(
  isTRUE(all.equal(strength_null$observed_strength, c(7, 6, 9))),
  isTRUE(all.equal(strength_null$null_mean, c(2, 2, 2))),
  isTRUE(all.equal(strength_null$null_sd, rep(sqrt(8), 3))),
  isTRUE(all.equal(
    strength_null$z_score, c(5, 4, 7) / sqrt(8)
  )),
  # Neither null iteration (strength 4 or 0) ever reaches the observed
  # strength (7, 6, 9), so exceedance_count = 0 and p = (1 + 0) / (2 + 1).
  isTRUE(all.equal(strength_null$exceedance_count, c(0L, 0L, 0L))),
  isTRUE(all.equal(strength_null$p_value, rep(1 / 3, 3)))
)

# exclude_self = FALSE keeps the diagonal in: observed becomes
# node1 = 6+3=9, node2 = 5+3=8, node3 = 3+8=11 (both null matrices already
# have a zero diagonal, so their null_mean/null_sd are unchanged at 2 and
# sqrt(8)).
strength_null_with_self <- strength_null_distribution(
  observed_counts, null_counts, progress_every = 0L, exclude_self = FALSE
)
stopifnot(
  isTRUE(all.equal(strength_null_with_self$observed_strength, c(9, 8, 11))),
  isTRUE(all.equal(strength_null_with_self$null_mean, c(2, 2, 2))),
  isTRUE(all.equal(strength_null_with_self$null_sd, rep(sqrt(8), 3))),
  isTRUE(all.equal(
    strength_null_with_self$z_score, c(7, 6, 9) / sqrt(8)
  )),
  isTRUE(all.equal(strength_null_with_self$exceedance_count, c(0L, 0L, 0L))),
  isTRUE(all.equal(strength_null_with_self$p_value, rep(1 / 3, 3)))
)

message("All statistical utility tests passed")

#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo_root <- if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."),
                mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(ccber))
source(file.path(repo_root, "R", "config.R"))
source(file.path(repo_root, "R", "graph_utils.R"))
source(file.path(repo_root, "R", "multiloop_utils.R"))

# Two recordings covering: a short gap (<250ms), a gap right at the
# threshold (>=250ms, new chain), a negative/overlapping gap (still
# chained), a long gap (new chain), and a trailing singleton.
whistles_list <- data.frame(
  recording = c(rep("recA", 4), rep("recB", 3)),
  whistle_type_chr = c(
    "SW_Neo_Category_1", "SW_Luna_Category_1",
    "SW_Yosefa_Category_1", "SW_Neo_Category_1",
    "SW_Neo_Category_1", "SW_Neo_Category_1", "SW_Luna_Category_1"
  ),
  start_time = c(0.0, 0.6, 1.5, 2.05, 10.0, 10.4, 15.0),
  end_time = c(0.5, 1.0, 2.0, 2.3, 10.5, 10.9, 15.2),
  whistle_type = c(1, 2, 3, 1, 1, 1, 2),
  stringsAsFactors = FALSE
)

result <- compute_multiloop_chains(whistles_list, threshold = 0.25)

stopifnot(
  nrow(result) == 7,
  identical(
    result$chain_id,
    c("recA_1", "recA_1", "recA_2", "recA_2", "recB_1", "recB_1", "recB_2")
  ),
  identical(result$chain_length, c(2L, 2L, 2L, 2L, 2L, 2L, 1L)),
  identical(result$position_in_chain, c(1L, 2L, 1L, 2L, 1L, 2L, 1L)),
  identical(
    result$category,
    c(
      "SW_Neo", "SW_Luna", "SW_Yosefa", "SW_Neo",
      "SW_Neo", "SW_Neo", "SW_Luna"
    )
  )
)

chains <- summarize_multiloop_chains(result)

stopifnot(
  nrow(chains) == 4,
  identical(chains$chain_id, c("recA_1", "recA_2", "recB_1", "recB_2")),
  identical(chains$chain_length, c(2L, 2L, 2L, 1L)),
  identical(
    chains$leading_category,
    c("SW_Neo", "SW_Yosefa", "SW_Neo", "SW_Luna")
  ),
  identical(chains$start_time, c(0.0, 1.5, 10.0, 15.0)),
  identical(chains$end_time, c(1.0, 2.3, 10.9, 15.2))
)

overlap <- compute_markov_transition_multiloop_overlap(
  whistles_list, time_window = c(0, 6), threshold = 0.25
)

# time_window = c(0, 6) is wide enough that every whistle sees all later
# whistles in its recording as candidates, so the only reason a transition
# is skipped or redirected is the overlap handling itself:
#  - recA: w1->w2 (gap .1, same chain), w2->w3 (gap .5, different chain),
#    w3->w4 (gap .05, same chain) - straightforward, no overlaps.
#  - recB: w5/w6 overlap (gap -0.1); since w7 is also a following candidate,
#    MarkovWhistle() skips the overlapping w6 and instead records w5->w7,
#    then separately records w6->w7. Note w5/w6 themselves (same multi-loop
#    chain) never appear as a recorded transition pair at all.
stopifnot(
  nrow(overlap) == 5,
  identical(overlap$recording, c(rep("recA", 3), rep("recB", 2))),
  identical(overlap$from_whistle_type, c(1, 2, 3, 1, 1)),
  identical(overlap$to_whistle_type, c(2, 3, 1, 2, 2)),
  identical(overlap$within_multiloop_chain, c(TRUE, FALSE, TRUE, FALSE, FALSE)),
  isTRUE(all.equal(overlap$transition_interval, c(0.1, 0.5, 0.05, 4.5, 4.1))),
  # from_row_id/to_row_id are 1-based positions in whistles_list as passed
  # in (recA occupies rows 1-4, recB rows 5-7 in the fixture above); recB's
  # skip-the-overlap case (w5->w7, w6->w7) is the case that would break a
  # naive "to_row_id = from_row_id + 1" assumption.
  identical(overlap$from_row_id, c(1L, 2L, 3L, 5L, 6L)),
  identical(overlap$to_row_id, c(2L, 3L, 4L, 7L, 7L))
)

self_repeats <- compute_self_repeat_counts(result)

# Only recB's w5->w6 (both whistle_type 1) is an adjacent same-type pair;
# recA_1 (type1 -> type2), recA_2 (type3 -> type1), and singleton recB_2
# contribute none.
stopifnot(
  nrow(self_repeats) == 1,
  self_repeats$whistle_type == 1L,
  self_repeats$n_self_repeats == 1L
)


# compute_multiloop_type_significance(): the sparse per-observed-edge
# exceedance counting only exists for performance (see the function's
# comment in R/multiloop_utils.R) - it must be statistically identical to
# the straightforward dense n_nodes x n_nodes computation of the same
# empirical enrichment test. Same seed before each run makes the
# iteration-by-iteration sample() draws (and therefore the null
# permutations themselves) bit-for-bit identical between the two
# implementations, so their p-values must match exactly, not just
# approximately.
set.seed(123)
chain_node_ids_sig <- sample(1:4, size = 40, replace = TRUE)
within_window_sig <- rep(TRUE, 39)
n_nodes_sig <- 4
from_idx_sig <- which(within_window_sig)
to_idx_sig <- from_idx_sig + 1L
counts_sig <- matrix(0L, n_nodes_sig, n_nodes_sig)
for (k in seq_along(from_idx_sig)) {
  fi <- chain_node_ids_sig[from_idx_sig[k]]
  ti <- chain_node_ids_sig[to_idx_sig[k]]
  counts_sig[fi, ti] <- counts_sig[fi, ti] + 1L
}

set.seed(456)
result_sparse <- compute_multiloop_type_significance(
  chain_node_ids_sig, within_window_sig, n_nodes_sig, counts_sig,
  iterations = 200, progress_every = 0
)

set.seed(456)
observed_probability_sig <- CalcTransitionMatrix(counts_sig)
exceedance_dense <- matrix(0L, n_nodes_sig, n_nodes_sig)
for (it in seq_len(200)) {
  shuffled <- sample(chain_node_ids_sig)
  null_counts <- matrix(0L, n_nodes_sig, n_nodes_sig)
  for (k in seq_along(from_idx_sig)) {
    fi <- shuffled[from_idx_sig[k]]
    ti <- shuffled[to_idx_sig[k]]
    null_counts[fi, ti] <- null_counts[fi, ti] + 1L
  }
  null_probability_sig <- CalcTransitionMatrix(null_counts)
  exceedance_dense <- exceedance_dense +
    (null_probability_sig >= observed_probability_sig)
}
p_value_dense <- (1 + exceedance_dense) / 201
p_value_dense[counts_sig == 0] <- 1

stopifnot(
  isTRUE(all.equal(result_sparse$p_value_matrix, p_value_dense)),
  all(result_sparse$p_value_matrix[counts_sig == 0] == 1),
  all(result_sparse$p_value_matrix > 0 & result_sparse$p_value_matrix <= 1)
)

message("All multi-loop utility tests passed")

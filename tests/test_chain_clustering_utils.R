#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo_root <- if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."),
                mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(repo_root, "R", "multiloop_utils.R"))
source(file.path(repo_root, "R", "chain_clustering_utils.R"))

# -- build_chain_substitution_cost -------------------------------------------

sub_levels <- c("A1", "A2", "B1")
category_of_sub <- c(A1 = "X", A2 = "X", B1 = "Y")
cost <- build_chain_substitution_cost(sub_levels, category_of_sub, epsilon = 0.3)

stopifnot(
  identical(dim(cost), c(3L, 3L)),
  identical(dimnames(cost), list(sub_levels, sub_levels)),
  all(diag(cost) == 0),
  cost["A1", "A2"] == 0.3,
  cost["A2", "A1"] == 0.3,
  cost["A1", "B1"] == 1,
  cost["B1", "A1"] == 1,
  cost["A2", "B1"] == 1
)

cost_half <- build_chain_substitution_cost(sub_levels, category_of_sub, epsilon = 0.5)
stopifnot(cost_half["A1", "A2"] == 0.5)

tryCatch({
  build_chain_substitution_cost(sub_levels, category_of_sub, epsilon = 1.5)
  stop("expected an error for epsilon out of [0, cross_category_cost]")
}, error = function(e) invisible(NULL))

# cross_category_cost overrides the default-1 "different main category" cost
# (e.g. the 0/0.3/0.6 scheme, where a cross-category swap is cheaper than
# the default gap_penalty = 1) - identical/same-category costs are untouched.
cost_06 <- build_chain_substitution_cost(
  sub_levels, category_of_sub, epsilon = 0.3, cross_category_cost = 0.6
)
stopifnot(
  all(diag(cost_06) == 0),
  cost_06["A1", "A2"] == 0.3,
  cost_06["A1", "B1"] == 0.6,
  cost_06["B1", "A1"] == 0.6,
  cost_06["A2", "B1"] == 0.6
)

tryCatch({
  build_chain_substitution_cost(
    sub_levels, category_of_sub, epsilon = 0.5, cross_category_cost = 0.3
  )
  stop("expected an error when epsilon exceeds cross_category_cost")
}, error = function(e) invisible(NULL))

tryCatch({
  build_chain_substitution_cost(
    sub_levels, category_of_sub, epsilon = 0.3, cross_category_cost = -1
  )
  stop("expected an error for a negative cross_category_cost")
}, error = function(e) invisible(NULL))

# -- nw_alignment_cost / nw_alignment_distance -------------------------------

# seq_a = A1 -> B1 (idx 1, 3); seq_b = A1 -> A2 -> B1 (idx 1, 2, 3).
# Optimal alignment: A1~A1 (0), insert A2 (gap=1), B1~B1 (0) -> cost 1.
seq_a_idx <- c(1L, 3L)
seq_b_idx <- c(1L, 2L, 3L)
raw_cost <- nw_alignment_cost(seq_a_idx, seq_b_idx, cost, gap_penalty = 1)
stopifnot(isTRUE(all.equal(raw_cost, 1)))

norm_dist <- nw_alignment_distance(seq_a_idx, seq_b_idx, cost, gap_penalty = 1)
stopifnot(isTRUE(all.equal(norm_dist, 1 / 3)))

# Identical sequences align for free.
stopifnot(
  nw_alignment_cost(seq_a_idx, seq_a_idx, cost, gap_penalty = 1) == 0,
  nw_alignment_distance(seq_a_idx, seq_a_idx, cost, gap_penalty = 1) == 0
)

# Empty vs. non-empty: pure gap cost, normalized distance == gap_penalty.
stopifnot(
  nw_alignment_cost(integer(0), seq_b_idx, cost, gap_penalty = 1) == 3,
  nw_alignment_distance(integer(0), seq_b_idx, cost, gap_penalty = 1) == 1
)

# A same-category substitution should be cheaper than a cross-category one
# over otherwise-identical sequences of equal length.
same_cat_cost <- nw_alignment_cost(c(1L, 1L), c(1L, 2L), cost, gap_penalty = 1)
cross_cat_cost <- nw_alignment_cost(c(1L, 1L), c(1L, 3L), cost, gap_penalty = 1)
stopifnot(same_cat_cost < cross_cat_cost)

# -- nw_alignment_ops ----------------------------------------------------------

# seq_a (A1,B1) is "the sequence being encoded"; seq_b (A1,A2,B1) is "the
# reference/medoid" - A1~A1 (match), then seq_b's extra "A2" is skipped (a
# DELETION: the medoid has content seq_a doesn't, no target label needed),
# then B1~B1 (match).
ops1 <- nw_alignment_ops(seq_a_idx, seq_b_idx, cost, gap_penalty = 1)
stopifnot(
  identical(ops1$counts, c(match = 2L, substitution = 0L, insertion = 0L, deletion = 1L)),
  identical(ops1$substitution_targets, integer(0)),
  identical(ops1$insertion_targets, integer(0))
)

# Same length, one substitution: A1~A1 (match), A1~B1 (substitution, target =
# seq_a's own label, index 1 = "A1") - cost 1, cheaper than the 2 indels
# (cost 2) a fully-unaligned reading would need.
ops2 <- nw_alignment_ops(c(1L, 1L), c(1L, 3L), cost, gap_penalty = 1)
stopifnot(
  identical(ops2$counts, c(match = 1L, substitution = 1L, insertion = 0L, deletion = 0L)),
  identical(ops2$substitution_targets, 1L),
  identical(ops2$insertion_targets, integer(0))
)

# Identical sequences align as all matches.
ops3 <- nw_alignment_ops(seq_a_idx, seq_a_idx, cost, gap_penalty = 1)
stopifnot(identical(
  ops3$counts, c(match = 2L, substitution = 0L, insertion = 0L, deletion = 0L)
))

# Empty seq_a (nothing being encoded) vs. non-empty seq_b: every seq_b
# position is skipped - all deletions, no insertions (seq_a has no content
# of its own to report) and so no target labels either way.
ops4 <- nw_alignment_ops(integer(0), seq_b_idx, cost, gap_penalty = 1)
stopifnot(
  identical(ops4$counts, c(match = 0L, substitution = 0L, insertion = 0L, deletion = 3L)),
  identical(ops4$substitution_targets, integer(0)),
  identical(ops4$insertion_targets, integer(0))
)

# -- chain_alignment_distance_matrix -----------------------------------------

chain_sequences <- list(
  p1 = c("A1", "B1"),
  p2 = c("A1", "A2", "B1"),
  p3 = c("B1", "B1")
)
dm <- chain_alignment_distance_matrix(chain_sequences, cost, gap_penalty = 1, progress_every = 0L)
stopifnot(
  inherits(dm, "dist"),
  attr(dm, "Size") == 3L,
  identical(attr(dm, "Labels"), c("p1", "p2", "p3"))
)
full <- as.matrix(dm)
stopifnot(
  isTRUE(all.equal(full["p1", "p2"], 1 / 3)),
  isTRUE(all.equal(full["p1", "p3"], nw_alignment_distance(c(1L, 3L), c(3L, 3L), cost, 1))),
  all(full == t(full)),
  all(diag(full) == 0)
)

# normalize = FALSE returns raw (unnormalized) cost - p1 (length 2) vs. p2
# (length 3) is 1/3 normalized but 1 raw (matches nw_alignment_cost()
# directly, not nw_alignment_distance()'s /longer version). This is what
# mdl_k_sweep() relies on for medoid selection: normalized- and raw-cost
# argmin can disagree once cluster members span more than one chain length,
# so medoid selection must use the raw matrix, not clustering$dist_obj.
dm_raw <- chain_alignment_distance_matrix(
  chain_sequences, cost, gap_penalty = 1, normalize = FALSE, progress_every = 0L
)
full_raw <- as.matrix(dm_raw)
stopifnot(
  isTRUE(all.equal(full_raw["p1", "p2"], 1)),
  isTRUE(all.equal(full_raw["p1", "p3"], nw_alignment_cost(c(1L, 3L), c(3L, 3L), cost, 1))),
  all(full_raw == t(full_raw)),
  all(diag(full_raw) == 0)
)

tryCatch({
  chain_alignment_distance_matrix(list(p1 = c("A1")), cost)
  stop("expected an error for fewer than two chain sequences")
}, error = function(e) invisible(NULL))

# -- silhouette_k_sweep -------------------------------------------------------

# Two well-separated triplets: {p1,p2,p3} and {p4,p5,p6}, within-group
# distance 0.1, cross-group distance 0.9 - k = 2 should clearly maximize
# mean silhouette width over k = 2..4.
labels <- paste0("p", 1:6)
full_dist <- matrix(0.9, 6, 6, dimnames = list(labels, labels))
full_dist[1:3, 1:3] <- 0.1
full_dist[4:6, 4:6] <- 0.1
diag(full_dist) <- 0
synthetic_dist <- as.dist(full_dist)

hc <- hclust(synthetic_dist, method = "average")
sweep <- silhouette_k_sweep(hc, synthetic_dist, k_range = 2:4)
stopifnot(
  identical(sweep$k, 2:4),
  sweep$k[which.max(sweep$mean_silhouette_width)] == 2
)

# -- cluster_multiloop_chains -------------------------------------------------

# Two composition-cluster fixtures: c1/c2/c3/c4 all built from A1/A2 (same
# main category "A"; c1 and c4 are exact duplicates, to check dedup carries
# counts), c5/c6/c7 all built from B1/B2 ("B"); c8 is a singleton
# (chain_length 1), excluded from clustering entirely. Cross-category
# distance is always 1 (maximal) while within-category distance is <= 0.3,
# so a k = 2 cut should cleanly separate A from B and score far better than
# k = 3 (which would needlessly split one of the two cohesive groups).
chain_fixture <- data.frame(
  chain_id = c(
    rep("c1", 2), rep("c2", 2), rep("c3", 2), rep("c4", 2),
    rep("c5", 2), rep("c6", 2), rep("c7", 2), "c8"
  ),
  chain_length = c(rep(2L, 14), 1L),
  position_in_chain = c(rep(1:2, 7), 1L),
  whistle_type_chr = c(
    "A1", "A1",
    "A1", "A2",
    "A2", "A1",
    "A1", "A1",
    "B1", "B1",
    "B1", "B2",
    "B2", "B1",
    "A2"
  ),
  stringsAsFactors = FALSE
)
fixture_sub_levels <- c("A1", "A2", "B1", "B2")
fixture_category_of_sub <- c(A1 = "A", A2 = "A", B1 = "B", B2 = "B")
fixture_cost <- build_chain_substitution_cost(
  fixture_sub_levels, fixture_category_of_sub, epsilon = 0.3
)

clustering <- cluster_multiloop_chains(
  chain_fixture, fixture_cost, gap_penalty = 1, linkage_method = "average",
  k_range = 2:3, progress_every = 0L
)

stopifnot(
  identical(sort(clustering$chain_level$chain_id), paste0("c", 1:7)),
  clustering$best_k == 2L,
  nrow(clustering$pattern_table) == 6L,
  # c1 and c4 are identical chains ("A1 -> A1") and must dedupe into one
  # pattern with n_chains == 2; every other chain is its own distinct
  # pattern with n_chains == 1.
  clustering$pattern_table$n_chains[
    clustering$pattern_table$ordered_pattern == "A1 -> A1"
  ] == 2L,
  sum(clustering$pattern_table$n_chains) == 7L
)

by_chain <- setNames(clustering$chain_level$cluster, clustering$chain_level$chain_id)
stopifnot(
  length(unique(by_chain[c("c1", "c2", "c3", "c4")])) == 1L,
  length(unique(by_chain[c("c5", "c6", "c7")])) == 1L,
  by_chain[["c1"]] != by_chain[["c5"]]
)

tryCatch({
  cluster_multiloop_chains(chain_fixture[chain_fixture$chain_id == "c8", ], fixture_cost)
  stop("expected an error when no multi-loop chains (chain_length >= 2) are present")
}, error = function(e) invisible(NULL))

# force_k cuts the same tree at a caller-chosen k instead of silhouette's
# argmax - forcing k = 6 (every one of the 6 distinct patterns its own
# cluster, the max valid force_k for n_patterns = 6) should give every chain
# a cluster matching its own pattern_id, unlike the silhouette default (k=2).
clustering_k6 <- cluster_multiloop_chains(
  chain_fixture, fixture_cost, gap_penalty = 1, linkage_method = "average",
  k_range = 2:3, progress_every = 0L, force_k = 5L
)
stopifnot(
  clustering_k6$best_k == 5L,
  # c1 and c4 share a pattern and so must still land in the same cluster
  # even forced this fine (5 clusters over 6 patterns - only one pair can
  # ever be merged).
  clustering_k6$chain_level$cluster[clustering_k6$chain_level$chain_id == "c1"] ==
    clustering_k6$chain_level$cluster[clustering_k6$chain_level$chain_id == "c4"],
  length(unique(clustering_k6$chain_level$cluster)) == 5L
)

tryCatch({
  cluster_multiloop_chains(chain_fixture, fixture_cost, force_k = 6L, progress_every = 0L)
  stop("expected an error for force_k >= n_patterns")
}, error = function(e) invisible(NULL))

tryCatch({
  cluster_multiloop_chains(chain_fixture, fixture_cost, force_k = 1L, progress_every = 0L)
  stop("expected an error for force_k < 2")
}, error = function(e) invisible(NULL))

# -- mdl_k_sweep -----------------------------------------------------------

# Same 6-pattern A/B fixture, swept from k = 1 (which silhouette/gap can't
# evaluate) through k = 5 (cluster::pam() itself rejects k == n_patterns,
# the trivial every-point-its-own-cluster case - see the k_range-clipping
# check below). The partition is re-optimized by PAM at every k (not
# inherited from one fixed tree), and each cluster's BIC penalty scales with
# how many distinct substitution-/insertion-target labels it actually
# involves, not a flat per-cluster fee - see mdl_k_sweep()'s own docstring
# for why an earlier flat-2-parameter version could never find an interior
# optimum. At k = 1 here, the one cluster mixes A- and B-family labels as
# substitution targets, so its adaptive label code (and its BIC dof) is
# expensive; splitting cleanly separates that diversity, so best_k = 2,
# matching silhouette_k_sweep()'s answer on the same well-separated A/B
# fixture.
mdl_result <- mdl_k_sweep(clustering, fixture_cost, k_range = 1:5, progress_every = 0L)

stopifnot(
  identical(mdl_result$sweep$k, 1:5),
  nrow(mdl_result$sweep) == 5L,
  mdl_result$best_k == 2L,
  all(mdl_result$sweep$model_bits > 0),
  all(mdl_result$sweep$data_bits >= 0),
  isTRUE(all.equal(
    mdl_result$sweep$total_bits,
    mdl_result$sweep$model_bits + mdl_result$sweep$data_bits
  )),
  # model_bits grows strictly with k (one more medoid + at least one more
  # BIC penalty term to pay for each time).
  all(diff(mdl_result$sweep$model_bits) > 0),
  # This fixture's data_bits is non-increasing too under PAM re-optimization.
  all(diff(mdl_result$sweep$data_bits) <= 1e-9)
)

# Deterministic - cluster::pam()'s default BUILD+SWAP is not randomized, and
# there's no other RNG involved.
mdl_result_repeat <- mdl_k_sweep(clustering, fixture_cost, k_range = 1:5, progress_every = 0L)
stopifnot(identical(mdl_result$sweep, mdl_result_repeat$sweep))

# k_range is silently clipped to [1, n_patterns) rather than erroring (note
# the open upper end - k == n_patterns is invalid for cluster::pam()), same
# convention as cluster_multiloop_chains()'s own k_range handling.
mdl_result_clipped <- mdl_k_sweep(clustering, fixture_cost, k_range = 0:100, progress_every = 0L)
stopifnot(identical(mdl_result_clipped$sweep$k, 1:5))

message("All chain clustering utility tests passed")

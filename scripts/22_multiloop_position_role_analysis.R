#!/usr/bin/env Rscript

# Response-letter analysis: replaces the whistle-level Spearman correlation
# between P(position | in loop) and network measures (scripts/
# 10_multiloop_position_subcategory.R's positions_vs_network_measure plots)
# with a group-comparison permutation test. Motivation (from the reviewer):
# a correlation computed across ALL nodes can be flat while a small set of
# hub nodes is specifically enriched at the terminal position - the
# reviewer's question is about those specific hub nodes, not about a
# population-wide linear trend. The per-node P(first)/P(last) computation
# itself (script 10's n_first / n_multiloop, n_last / n_multiloop) is
# unchanged; what changes is how it's tested against network role.
#
# NODE SETS
#   Set A: every sub-category with >= 10 occurrences inside multi-loop
#     chains (chain_length >= 2) that also has a node in the BH-significant
#     whistle-level network (scripts/01/02) - i.e. exactly script 10's
#     min_n_for_plot filter on n_multiloop, reproducing its n = 51.
#   Set B: Set A restricted to nodes above the 90th percentile of EITHER
#     betweenness or standardized strength in Figure 5c
#     (outputs/plots/Figure_network_analysis.png/svg, panel c =
#     outputs/plots/markov_betweenness_vs_strength_zscore.png, produced by
#     plot_betweenness_strength_scatter() in R/graph_utils.R). That function
#     normalizes betweenness by dividing by its own max (equivalent to
#     min-max since betweenness has a hard floor at 0) and min-max
#     normalizes the z-score strength (it can be negative), then thresholds
#     each at quantile(., 0.9) - reproduced here from outputs/node_measures.csv
#     rather than re-deriving betweenness/strength from scratch, since that
#     CSV is the exact data plot_betweenness_strength_scatter() was called
#     on (scripts/02_node_network_measures.R, same script run: both files'
#     mtimes are seconds apart).
#
#   The instructions accompanying this analysis included node numbers "read
#   off the figure" as a manual cross-check:
#     high strength:    21, 26, 46, 16, 23, 27, 8, 1
#     high betweenness: 26, 44, 34, 29, 45, 42, 41, 38, 55
#   Recomputing from node_measures.csv (below) does NOT exactly reproduce
#   these - see the printed discrepancy report. Cropping the actual
#   manuscript panel (outputs/plots/Figure_network_analysis.png, panel c)
#   and reading its point labels directly confirms the RECOMPUTED sets
#   (high strength adds node 83; high betweenness is 20/58/64 in place of
#   29/44/55), so the recomputed sets - not the typed reference lists - are
#   used as Set B's membership below. This is reported, not silently
#   substituted.
#
#   Turn-taking's "high" criterion (Figure 5b,
#   outputs/plots/network_measures/turn_taking.png) is not a percentile: the
#   figure's own dashed reference line is a fixed threshold at
#   prob_diff_cat_outbound > 0.5, reused as-is.
#
# TEST
#   Within each node set (A, B), nodes are split into "high" vs "rest" under
#   each of 3 role definitions (strength, betweenness, turn-taking) - for
#   Set B, "rest" means "in Set B (i.e. a hub by SOME measure) but not high
#   on this particular role", not "not a hub at all". Null: permute the
#   assignment of sub-category labels to positions WITHIN each multi-loop
#   chain, independently per chain, 10,000 times - identical in mechanism to
#   the within-chain shuffle null already used in script 10 (see its
#   "permutation null: shuffle whistle_type_chr within chain_id" section),
#   which preserves chain length, each chain's own sub-category multiset,
#   and hence every node's total in-chain occurrence count (n_multiloop)
#   exactly; only which position each occurrence lands in is randomized.
#   Statistic: delta = mean P(last) of the "high" group - mean P(last) of
#   "rest"; two-sided empirical p = (1 + #{|delta_null| >= |delta_obs|}) /
#   (iterations + 1) (the same "+1" smoothing convention script 10 already
#   uses for its permutation p-values). Repeated identically for P(first)
#   [introductory] and P(middle) [central] - each of the three positions gets
#   its own 6-test (role x node set) family, BH-corrected within itself.
#
# OUTPUT (outputs/multiloop_analysis/position_analysis/)
#   position_per_node.csv               - one row per Set A node
#   position_group_tests.csv            - 18 rows: {first, middle, last} x
#                                          3 roles x 2 node sets, each
#                                          position's 6 BH-corrected within
#                                          itself
#   position_role_definitions.csv       - Fig 5b/5c high-node sets + the
#                                          reference-list discrepancy
#   position_diagnostics_excl_L2.csv    - all three positions' tests re-run
#                                          on chain_length >= 3 chains only
#   position_diagnostics_heterogeneous.csv - all three positions re-run on
#                                          heterogeneous-only chains, IF
#                                          most chains are homogeneous
#   position_diagnostics_summary.txt    - diagnostics 1/2/4, dropped-node
#                                          notes, and the 5 quoted nodes

suppressPackageStartupMessages({
  library(stats)
})

options(width = 200) # keeps the group-test tables from wrapping in printed/captured output

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  REPO_ROOT <- normalizePath(
    file.path(dirname(sub("^--file=", "", file_arg[1])), ".."),
    mustWork = TRUE
  )
} else {
  REPO_ROOT <- normalizePath(getwd(), mustWork = TRUE)
}
setwd(REPO_ROOT)

source(file.path(REPO_ROOT, "R", "config.R"))
source(file.path(REPO_ROOT, "R", "load_whistles.R"))
source(file.path(REPO_ROOT, "R", "graph_utils.R"))
source(file.path(REPO_ROOT, "R", "multiloop_utils.R"))

position_dir <- file.path(REPO_ROOT, OUTPUT_DIR, "multiloop_analysis", "position_analysis")
dir.create(position_dir, recursive = TRUE, showWarnings = FALSE)

ITERATIONS <- 10000L
MIN_N_MULTILOOP <- 10L # unchanged from script 10's min_n_for_plot, applied to n_multiloop
QUOTE_NODES <- c(21L, 1L, 36L, 26L, 46L)

# ============================================================================
# 1. Load whistles / multi-loop chains; per-occurrence position (UNCHANGED
#    from scripts/10_multiloop_position_subcategory.R)
# ============================================================================

message("Loading whistles from ", DATA_PATH)
whistles_list <- load_whistles(DATA_PATH)

message(
  "Grouping whistles into multi-loop chains (IWI < ",
  MULTILOOP_IWI_THRESHOLD * 1000, "ms)..."
)
whistle_chains <- compute_multiloop_chains(
  whistles_list, threshold = MULTILOOP_IWI_THRESHOLD, list_names = LIST_NAMES
)
ml <- whistle_chains[whistle_chains$chain_length >= 2, ]
ml$position_category <- ifelse(
  ml$position_in_chain == 1L, "first",
  ifelse(ml$position_in_chain == ml$chain_length, "last", "middle")
)

node_measures_csv <- file.path(REPO_ROOT, OUTPUT_DIR, "node_measures.csv")
if (!file.exists(node_measures_csv)) {
  stop("Run scripts/02_node_network_measures.R first. Missing: ", node_measures_csv, call. = FALSE)
}
node_measures <- read.csv(node_measures_csv, stringsAsFactors = FALSE)
node_measures$node_id <- as.integer(node_measures$nodes)

# node_id (integer, Markov-model numbering) <-> whistle_type_chr (sub-category
# string) lookup, both already columns of whistle_chains.
node_lookup <- unique(whistle_chains[, c("whistle_type", "whistle_type_chr", "category")])
chr_by_node_id <- setNames(node_lookup$whistle_type_chr, as.character(node_lookup$whistle_type))
node_id_by_chr <- setNames(node_lookup$whistle_type, node_lookup$whistle_type_chr)

all_sub_levels <- sort(unique(whistle_chains$whistle_type_chr))
n_sub_all <- length(all_sub_levels)
node_id_all <- as.integer(node_id_by_chr[all_sub_levels])
category_all <- node_lookup$category[match(all_sub_levels, node_lookup$whistle_type_chr)]

# ============================================================================
# 2. Per-node baseline descriptives from the FULL chain_length >= 2 data:
#    n_multiloop (occurrence count - unchanged from script 10),
#    n_chains (DISTINCT chain count - new), mean_chain_length (new).
# ============================================================================

sub_code_full <- match(ml$whistle_type_chr, all_sub_levels)
n_multiloop_full <- tabulate(sub_code_full, nbins = n_sub_all)

# One row per (chain_id, node) actually present in a chain, deduplicating a
# node that repeats within one chain - R/multiloop_utils.R's own
# compute_multiloop_type_membership(), already built for exactly this.
membership <- compute_multiloop_type_membership(whistle_chains)
chain_length_by_id <- whistle_chains$chain_length[
  match(membership$chain_id, whistle_chains$chain_id)
]
membership$chain_length <- chain_length_by_id
membership_sub_code <- match(as.integer(membership$whistle_type), node_id_all)

n_chains_full <- tabulate(membership_sub_code, nbins = n_sub_all)
mean_chain_length_full <- as.vector(tapply(
  membership$chain_length,
  factor(membership_sub_code, levels = seq_len(n_sub_all)),
  mean
))

# ============================================================================
# 3. Node role definitions - Figure 5b (turn-taking) and Figure 5c
#    (betweenness / standardized strength), reproduced from node_measures.csv
#    exactly as plot_betweenness_strength_scatter() computes them.
# ============================================================================

betweenness_norm <- node_measures$betweenness / max(node_measures$betweenness)
strength_range <- range(node_measures$z_score_incl_self)
strength_norm <- (node_measures$z_score_incl_self - strength_range[1]) / diff(strength_range)
betweenness_threshold <- quantile(betweenness_norm, 0.9, names = FALSE)
strength_threshold <- quantile(strength_norm, 0.9, names = FALSE)

high_strength_ids <- sort(node_measures$node_id[strength_norm > strength_threshold])
high_betweenness_ids <- sort(node_measures$node_id[betweenness_norm > betweenness_threshold])
high_turntaking_ids <- sort(node_measures$node_id[node_measures$prob_diff_cat_outbound > 0.5])

ref_strength_ids <- c(21, 26, 46, 16, 23, 27, 8, 1)
ref_betweenness_ids <- c(26, 44, 34, 29, 45, 42, 41, 38, 55)

strength_extra_recomputed <- setdiff(high_strength_ids, ref_strength_ids)
strength_extra_reference <- setdiff(ref_strength_ids, high_strength_ids)
betweenness_extra_recomputed <- setdiff(high_betweenness_ids, ref_betweenness_ids)
betweenness_extra_reference <- setdiff(ref_betweenness_ids, high_betweenness_ids)

message("\n================ Figure 5b/5c role-definition check ================")
message(
  "Recomputed high-strength (z-score, min-max norm, >90th pct = ",
  round(strength_threshold, 4), "), n=", length(high_strength_ids), ": ",
  paste(high_strength_ids, collapse = ", ")
)
message(
  "Recomputed high-betweenness (normalized, >90th pct = ",
  round(betweenness_threshold, 4), "), n=", length(high_betweenness_ids), ": ",
  paste(high_betweenness_ids, collapse = ", ")
)
if (length(strength_extra_recomputed) == 0 && length(strength_extra_reference) == 0) {
  message("high-strength: recomputed set matches the reference list exactly.")
} else {
  message(
    "DISCREPANCY high-strength: recomputed-only = {",
    paste(strength_extra_recomputed, collapse = ", "), "}; reference-only = {",
    paste(strength_extra_reference, collapse = ", "), "}. ",
    "Cropping outputs/plots/Figure_network_analysis.png panel c and reading ",
    "its point labels directly confirms the RECOMPUTED set (node 83 is ",
    "labeled above the strength threshold line in the actual figure) - the ",
    "reference list appears to have omitted it. Using the recomputed set."
  )
}
if (length(betweenness_extra_recomputed) == 0 && length(betweenness_extra_reference) == 0) {
  message("high-betweenness: recomputed set matches the reference list exactly.")
} else {
  message(
    "DISCREPANCY high-betweenness: recomputed-only = {",
    paste(betweenness_extra_recomputed, collapse = ", "), "}; reference-only = {",
    paste(betweenness_extra_reference, collapse = ", "), "}. ",
    "Reference-only nodes' actual normalized betweenness (recomputed): ",
    paste0(
      betweenness_extra_reference, "=",
      round(betweenness_norm[match(betweenness_extra_reference, node_measures$node_id)], 4),
      collapse = ", "
    ),
    " - all far below the ", round(betweenness_threshold, 4), " threshold (node 29 has raw ",
    "betweenness = 0), so this isn't a marginal rounding difference. Cropping ",
    "outputs/plots/Figure_network_analysis.png panel c and reading its point ",
    "labels directly confirms the RECOMPUTED set (20, 58, 64 are the labeled ",
    "points beyond the betweenness threshold line, not 29/44/55). Using the ",
    "recomputed set; the reference list appears to be a transcription error."
  )
}
message(
  "Turn-taking high (Figure 5b fixed criterion, prob_diff_cat_outbound > 0.5), n=",
  length(high_turntaking_ids), ": ", paste(high_turntaking_ids, collapse = ", ")
)
message("======================================================================\n")

role_definition_rows <- rbind(
  data.frame(role = "strength", node_id = high_strength_ids, source = "recomputed"),
  data.frame(role = "betweenness", node_id = high_betweenness_ids, source = "recomputed"),
  data.frame(role = "turn_taking", node_id = high_turntaking_ids, source = "fixed_0.5_threshold"),
  data.frame(role = "strength", node_id = ref_strength_ids, source = "reference_list_as_given"),
  data.frame(role = "betweenness", node_id = ref_betweenness_ids, source = "reference_list_as_given")
)
write.csv(
  role_definition_rows, file.path(position_dir, "position_role_definitions.csv"),
  row.names = FALSE
)
message("Saved ", file.path(position_dir, "position_role_definitions.csv"))

# ============================================================================
# 4. Node sets A and B
# ============================================================================

in_significant_network <- node_id_all %in% node_measures$node_id
setA_mask <- (n_multiloop_full >= MIN_N_MULTILOOP) & in_significant_network
setA_ids <- sort(node_id_all[setA_mask])
message("Set A (n_multiloop >= ", MIN_N_MULTILOOP, ", in BH-significant network): n = ", length(setA_ids))
if (length(setA_ids) != 51L) {
  message(
    "NOTE: Set A size (", length(setA_ids), ") differs from the n = 51 quoted in the ",
    "instructions - underlying data may have changed since that number was reported."
  )
}

setB_raw_ids <- sort(union(high_strength_ids, high_betweenness_ids))
setB_ids <- intersect(setB_raw_ids, setA_ids)
setB_dropped_ids <- sort(setdiff(setB_raw_ids, setA_ids))
setB_dropped_frac <- length(setB_dropped_ids) / length(setB_raw_ids)

message(
  "Set B raw pool (Figure 5c, high strength OR high betweenness), n = ",
  length(setB_raw_ids), ": ", paste(setB_raw_ids, collapse = ", ")
)
if (length(setB_dropped_ids) > 0) {
  dropped_n_multiloop <- n_multiloop_full[match(chr_by_node_id[as.character(setB_dropped_ids)], all_sub_levels)]
  message(
    "Set B nodes dropped for n_multiloop < ", MIN_N_MULTILOOP, " (", length(setB_dropped_ids),
    " / ", length(setB_raw_ids), " = ", round(100 * setB_dropped_frac, 1), "%): ",
    paste0(setB_dropped_ids, " (n_multiloop=", dropped_n_multiloop, ")", collapse = ", ")
  )
}
if (setB_dropped_frac > 0.5) {
  message(
    "Set B: more than half the Figure-5c high-centrality pool has < ", MIN_N_MULTILOOP,
    " multi-loop occurrences and cannot be assigned a stable P(first)/P(last) - Set B is ",
    "underpowered by construction. Proceeding is still requested for transparency, but this ",
    "should be read as exploratory, not confirmatory."
  )
} else {
  message(
    "Set B usable pool after the n_multiloop filter: n = ", length(setB_ids),
    " / ", length(setB_raw_ids), " (", round(100 * (1 - setB_dropped_frac), 1),
    "% retained - not \"most of the group\", so the test proceeds)."
  )
}
message("Set B (usable): ", paste(setB_ids, collapse = ", "), "\n")

node_sets <- list(SetA = setA_ids, SetB = setB_ids)
role_high_sets <- list(
  strength = high_strength_ids, betweenness = high_betweenness_ids, turn_taking = high_turntaking_ids
)

# Report the composition of each (node set x role) split up front - a group
# of 1 isn't just underpowered, it's degenerate (its "mean" is one node's
# own value), so bh_correct_family() below drops that cell from its BH
# family entirely rather than correcting across it.
message("Group sizes for each (node set x role) split:")
for (set_name in names(node_sets)) {
  universe <- node_sets[[set_name]]
  for (role_name in names(role_high_sets)) {
    hi <- intersect(universe, role_high_sets[[role_name]])
    rest <- setdiff(universe, role_high_sets[[role_name]])
    flag <- if (length(hi) <= 1 || length(rest) <= 1) " <- DEGENERATE (excluded from BH family)" else ""
    message(
      "  ", set_name, " x ", role_name, ": high n=", length(hi),
      ", rest n=", length(rest), flag
    )
  }
}
message("")

# ============================================================================
# 5. Permutation machinery
#    Shuffles WHICH sub-category label occupies each chain slot, uniformly
#    at random within each chain (order(chain_code, runif(n)) is the same
#    idiom script 10 already uses for its within-chain shuffle null) - this
#    leaves chain length, each chain's own label multiset, and hence every
#    node's n_multiloop exactly fixed; only per-node n_first/n_last (and so
#    P(first)/P(last)) vary across draws.
# ============================================================================

compute_observed_positions <- function(ml_subset, sub_levels) {
  n_sub <- length(sub_levels)
  sub_code <- match(ml_subset$whistle_type_chr, sub_levels)
  is_first <- ml_subset$position_category == "first"
  is_last <- ml_subset$position_category == "last"
  n_multiloop <- tabulate(sub_code, nbins = n_sub)
  n_first <- tabulate(sub_code[is_first], nbins = n_sub)
  n_last <- tabulate(sub_code[is_last], nbins = n_sub)
  n_middle <- n_multiloop - n_first - n_last
  list(
    sub_levels = sub_levels, n_multiloop = n_multiloop,
    n_first = n_first, n_last = n_last, n_middle = n_middle,
    P_first = n_first / n_multiloop, P_last = n_last / n_multiloop,
    P_middle = n_middle / n_multiloop,
    # Analytic chance level for P(first) = P(last) under this same null (see
    # script 10's derivation: E[1/chain_length], averaged over each node's
    # own in-chain occurrences) - used below as the per-node test's null
    # center, and as a sanity check against the empirical permutation mean.
    p_chance = as.vector(tapply(
      1 / ml_subset$chain_length, factor(sub_code, levels = seq_len(n_sub)), mean
    ))
  )
}

run_permutation_null <- function(ml_subset, sub_levels, iterations, seed) {
  n_sub <- length(sub_levels)
  chain_code <- as.integer(factor(ml_subset$chain_id))
  sort_idx <- order(chain_code)
  chain_code_s <- chain_code[sort_idx]
  sub_code_s <- match(ml_subset$whistle_type_chr, sub_levels)[sort_idx]
  is_first_s <- ml_subset$position_category[sort_idx] == "first"
  is_last_s <- ml_subset$position_category[sort_idx] == "last"
  n_rows <- length(sub_code_s)
  n_multiloop <- tabulate(sub_code_s, nbins = n_sub)

  null_P_first <- matrix(NA_real_, nrow = iterations, ncol = n_sub)
  null_P_last <- matrix(NA_real_, nrow = iterations, ncol = n_sub)

  set.seed(seed)
  progress_every <- 2000L
  for (b in seq_len(iterations)) {
    ord <- order(chain_code_s, runif(n_rows))
    shuffled <- sub_code_s[ord]
    null_P_first[b, ] <- tabulate(shuffled[is_first_s], nbins = n_sub) / n_multiloop
    null_P_last[b, ] <- tabulate(shuffled[is_last_s], nbins = n_sub) / n_multiloop
    if (b %% progress_every == 0L || b == iterations) {
      message("  permutation ", b, " / ", iterations)
    }
  }
  # P(middle) isn't drawn separately - first/middle/last partition every
  # occurrence in every permutation (the shuffle only reassigns labels to
  # fixed slots), so P_middle = 1 - P_first - P_last exactly, every draw.
  null_P_middle <- 1 - null_P_first - null_P_last
  list(
    sub_levels = sub_levels, n_multiloop = n_multiloop,
    null_P_first = null_P_first, null_P_middle = null_P_middle, null_P_last = null_P_last
  )
}

# node_ids (integer) -> column index into a sub_levels-aligned vector/matrix
node_ids_to_idx <- function(node_ids, sub_levels) {
  match(chr_by_node_id[as.character(node_ids)], sub_levels)
}

emp_p <- function(null_vals, observed) (1 + sum(abs(null_vals) >= abs(observed))) / (length(null_vals) + 1)

# One row per (node_set x role), for one position ("first", "middle", or
# "last"): observed delta, and its permutation p from the pre-computed null
# matrix (P_obs/null_mat must already be the pair for that position).
# Also reports each group's observed mean AND its own null mean (mean of the
# group-mean's permutation distribution) side by side with the delta, so a
# reader can see whether an effect comes from the "high" group sitting above
# ITS null level, the "rest" group sitting below, or both - and n_retained_*
# (nodes with a finite P_obs, i.e. actually present in this ml_subset) in
# case a restricted subset (e.g. L>=3 only) leaves some group members with
# zero occurrences and hence an undefined P(position).
group_tests_for_position <- function(P_obs, null_mat, sub_levels, node_sets, role_high_sets, position_label) {
  rows <- list()
  for (set_name in names(node_sets)) {
    universe <- node_sets[[set_name]]
    for (role_name in names(role_high_sets)) {
      high_ids <- intersect(universe, role_high_sets[[role_name]])
      rest_ids <- setdiff(universe, role_high_sets[[role_name]])
      high_idx <- node_ids_to_idx(high_ids, sub_levels)
      rest_idx <- node_ids_to_idx(rest_ids, sub_levels)

      high_vals <- P_obs[high_idx]
      rest_vals <- P_obs[rest_idx]
      n_retained_high <- sum(is.finite(high_vals))
      n_retained_rest <- sum(is.finite(rest_vals))
      mean_high <- mean(high_vals, na.rm = TRUE)
      mean_rest <- mean(rest_vals, na.rm = TRUE)
      delta_obs <- mean_high - mean_rest

      null_high <- if (length(high_idx) > 0) {
        rowMeans(null_mat[, high_idx, drop = FALSE], na.rm = TRUE)
      } else {
        rep(NA_real_, nrow(null_mat))
      }
      null_rest <- if (length(rest_idx) > 0) {
        rowMeans(null_mat[, rest_idx, drop = FALSE], na.rm = TRUE)
      } else {
        rep(NA_real_, nrow(null_mat))
      }
      null_delta <- null_high - null_rest
      p_val <- if (n_retained_high > 0 && n_retained_rest > 0) emp_p(null_delta, delta_obs) else NA_real_

      rows[[length(rows) + 1]] <- data.frame(
        node_set = set_name, role = role_name, position = position_label,
        n_high = length(high_idx), n_rest = length(rest_idx),
        n_retained_high = n_retained_high, n_retained_rest = n_retained_rest,
        mean_high = round(mean_high, 4), mean_rest = round(mean_rest, 4),
        null_mean_high = round(mean(null_high, na.rm = TRUE), 4),
        null_mean_rest = round(mean(null_rest, na.rm = TRUE), 4),
        delta_observed = round(delta_obs, 4), p_raw = signif(p_val, 4),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# A group of 1 isn't underpowered so much as degenerate - its "mean" is a
# single node's own value, not a group statistic - so it's dropped from the
# multiple-testing family entirely (BH corrected over the remaining 5)
# rather than diluting the correction with a test that can't really answer
# the group-comparison question. The row itself is kept (p_raw, delta, group
# sizes all still visible) with p_bh left NA and excluded_from_bh = TRUE, so
# nothing silently disappears - it's just not part of what gets corrected.
# Currently this is exactly SetB x turn_taking (n_high = 1); written as a
# general n < 2 rule rather than a one-off case so it generalizes correctly
# if a future re-run (different node sets, different data) produces another
# degenerate cell.
bh_correct_family <- function(df) {
  eligible <- df$n_high >= 2 & df$n_rest >= 2
  df$excluded_from_bh <- !eligible
  df$p_bh <- NA_real_
  if (any(eligible)) {
    df$p_bh[eligible] <- signif(stats::p.adjust(df$p_raw[eligible], method = "BH"), 4)
  }
  df
}

run_full_test_suite <- function(ml_subset, sub_levels, iterations, seed, label) {
  message("Running ", iterations, " within-chain permutations (", label, ")...")
  observed <- compute_observed_positions(ml_subset, sub_levels)
  null <- run_permutation_null(ml_subset, sub_levels, iterations, seed)
  # last/introductory-first/middle-central each get their own 6-candidate
  # (role x node set) family, BH-corrected within itself over whichever of
  # those 6 are eligible (see bh_correct_family) - "primary"/"secondary"/
  # "tertiary" naming kept only for the pre-existing P(last)/P(first) fields
  # other sections already refer to; all three are run and reported
  # identically now.
  primary <- bh_correct_family(
    group_tests_for_position(observed$P_last, null$null_P_last, sub_levels, node_sets, role_high_sets, "last")
  )
  secondary <- bh_correct_family(
    group_tests_for_position(observed$P_first, null$null_P_first, sub_levels, node_sets, role_high_sets, "first")
  )
  tertiary <- bh_correct_family(
    group_tests_for_position(observed$P_middle, null$null_P_middle, sub_levels, node_sets, role_high_sets, "middle")
  )
  list(observed = observed, null = null, primary = primary, secondary = secondary, tertiary = tertiary)
}

# ============================================================================
# 6. DIAGNOSTICS (reported before results - they determine interpretability)
# ============================================================================

diag_lines <- c("Diagnostics", "===========", "")

# --- 1. Chain length distribution -------------------------------------------
chain_level <- whistle_chains[
  !duplicated(whistle_chains$chain_id) & whistle_chains$chain_length >= 2,
  c("chain_id", "chain_length")
]
n_chains_total <- nrow(chain_level)
len_bucket <- cut(
  chain_level$chain_length, breaks = c(1, 2, 3, 4, Inf),
  labels = c("L=2", "L=3", "L=4", "L=5+"), right = TRUE
)
len_table <- table(len_bucket)
diag_lines <- c(
  diag_lines,
  "1. Chain length distribution (chains with length >= 2):",
  paste0("   n chains total: ", n_chains_total),
  paste0(
    "   ", names(len_table), ": ", as.vector(len_table),
    " (", round(100 * as.vector(len_table) / n_chains_total, 1), "%)"
  ),
  ""
)
message(paste(diag_lines, collapse = "\n"))

# --- 2. Fraction of chains with >= 2 distinct sub-categories ---------------
chain_n_distinct <- tapply(ml$whistle_type_chr, ml$chain_id, function(x) length(unique(x)))
heterogeneous_chain_ids <- names(chain_n_distinct)[chain_n_distinct >= 2]
frac_heterogeneous <- length(heterogeneous_chain_ids) / length(chain_n_distinct)
most_homogeneous <- frac_heterogeneous < 0.5
diag2_lines <- c(
  paste0(
    "2. Chains with >= 2 distinct sub-categories: ", length(heterogeneous_chain_ids),
    " / ", length(chain_n_distinct), " (", round(100 * frac_heterogeneous, 1), "%)"
  ),
  paste0(
    "   Homogeneous chains (single sub-category, permutation is degenerate for these - ",
    "shuffling changes nothing observable): ", length(chain_n_distinct) - length(heterogeneous_chain_ids),
    " (", round(100 * (1 - frac_heterogeneous), 1), "%)"
  ),
  if (most_homogeneous) {
    "   Most chains are homogeneous -> re-running the full test suite on heterogeneous-only chains below."
  } else {
    "   Most chains are heterogeneous -> the primary (all-chains) run is already representative; heterogeneous-only re-run skipped."
  },
  ""
)
diag_lines <- c(diag_lines, diag2_lines)
message(paste(diag2_lines, collapse = "\n"))

# ============================================================================
# 7. PRIMARY run (all chains, chain_length >= 2)
# ============================================================================

primary_run <- run_full_test_suite(ml, all_sub_levels, ITERATIONS, SEED, "primary, all chains")

message("\n---- P(last) [terminal] group tests: 3 roles x 2 node sets (BH across the 6) ----")
print(primary_run$primary, row.names = FALSE)
message("\n---- P(first) [introductory] group tests: 3 roles x 2 node sets (BH across the 6) ----")
print(primary_run$secondary, row.names = FALSE)
message("\n---- P(middle) [central] group tests: 3 roles x 2 node sets (BH across the 6) ----")
print(primary_run$tertiary, row.names = FALSE)

# --- 3. Re-run all three positions excluding L=2 chains ---------------------
# At L=2, first and last are the only two slots (no middle at all - P(middle)
# computed on L=2-only occurrences would be identically 0), so for a node's
# occurrences restricted to L=2 chains, P(first) + P(last) = 1 exactly and
# P(middle) = 0 exactly - all three are mechanically constrained there. This
# re-run checks whether the full-data result depends on that.
ml_L3 <- ml[ml$chain_length >= 3, ]
message(
  "\n3. Re-running first/middle/last excluding L=2 chains: ",
  length(unique(ml_L3$chain_id)), " chains remain (of ", n_chains_total, ")."
)
excl_L2_run <- run_full_test_suite(ml_L3, all_sub_levels, ITERATIONS, SEED, "diagnostic, chain_length >= 3 only")
message("\n---- Diagnostic (chain_length >= 3 only): P(last) ----")
print(excl_L2_run$primary, row.names = FALSE)
message("\n---- Diagnostic (chain_length >= 3 only): P(first) ----")
print(excl_L2_run$secondary, row.names = FALSE)
message("\n---- Diagnostic (chain_length >= 3 only): P(middle) ----")
print(excl_L2_run$tertiary, row.names = FALSE)
write.csv(
  rbind(excl_L2_run$primary, excl_L2_run$secondary, excl_L2_run$tertiary),
  file.path(position_dir, "position_diagnostics_excl_L2.csv"), row.names = FALSE
)
message("Saved ", file.path(position_dir, "position_diagnostics_excl_L2.csv"))

# --- heterogeneous-only re-run (conditional on diagnostic 2) ---------------
heterogeneous_run <- NULL
if (most_homogeneous) {
  ml_heterogeneous <- ml[ml$chain_id %in% heterogeneous_chain_ids, ]
  message(
    "\nRe-running full suite (P(last) + P(first) + P(middle)) on heterogeneous-only chains: ",
    length(heterogeneous_chain_ids), " chains."
  )
  heterogeneous_run <- run_full_test_suite(
    ml_heterogeneous, all_sub_levels, ITERATIONS, SEED, "diagnostic, heterogeneous chains only"
  )
  message("\n---- Diagnostic (heterogeneous chains only): P(last) ----")
  print(heterogeneous_run$primary, row.names = FALSE)
  message("\n---- Diagnostic (heterogeneous chains only): P(first) ----")
  print(heterogeneous_run$secondary, row.names = FALSE)
  message("\n---- Diagnostic (heterogeneous chains only): P(middle) ----")
  print(heterogeneous_run$tertiary, row.names = FALSE)
  write.csv(
    rbind(heterogeneous_run$primary, heterogeneous_run$secondary, heterogeneous_run$tertiary),
    file.path(position_dir, "position_diagnostics_heterogeneous.csv"), row.names = FALSE
  )
  message("Saved ", file.path(position_dir, "position_diagnostics_heterogeneous.csv"))
}

# --- 4. Correlation between n_multiloop (occurrence count) and each of
#        P(first)/P(middle)/P(last) -----------------------------------------
setA_idx <- node_ids_to_idx(setA_ids, all_sub_levels)
n_for_corr <- primary_run$observed$n_multiloop[setA_idx]
diag4_lines <- c(
  "4. Spearman correlation, n_multiloop (occurrence count) vs. position rate, across Set A:",
  paste0(
    "   (tests whether a position rate regresses toward the mean for nodes with more multi-loop ",
    "occurrences - a large |rho| here would mean frequent nodes systematically pull toward the ",
    "population average, which would matter for interpreting any node-level outlier.)"
  )
)
for (pos in c("first", "middle", "last")) {
  P_vec <- switch(pos,
    first = primary_run$observed$P_first, middle = primary_run$observed$P_middle,
    last = primary_run$observed$P_last
  )[setA_idx]
  freq_cor <- suppressWarnings(cor.test(n_for_corr, P_vec, method = "spearman"))
  diag4_lines <- c(
    diag4_lines,
    paste0(
      "   P(", pos, "): rho = ", round(unname(freq_cor$estimate), 3),
      ", p = ", signif(unname(freq_cor$p.value), 3), " (n=", length(setA_idx), ")"
    )
  )
}
diag4_lines <- c(diag4_lines, "")
diag_lines <- c(diag_lines, diag4_lines)
message(paste(diag4_lines, collapse = "\n"))

# ============================================================================
# 8. PER-NODE supplementary table (Set A), with a per-node permutation p for
#    each of P(first)/P(middle)/P(last): two-sided, comparing how far the
#    OBSERVED rate sits from that node's own analytic chance level
#    (p_chance = E[1/chain_length] for first/last, 1 - 2*p_chance for
#    middle, the same null already derived in script 10) against how far the
#    null draws sit from that same chance level.
# ============================================================================

obs <- primary_run$observed
per_node_perm_p <- function(P_obs_vec, null_mat, chance_vec, idx) {
  vapply(idx, function(j) {
    dev_obs <- P_obs_vec[j] - chance_vec[j]
    dev_null <- null_mat[, j] - chance_vec[j]
    emp_p(dev_null, dev_obs)
  }, numeric(1))
}
p_chance_middle <- 1 - 2 * obs$p_chance # middle's chance level under the same null (first/last symmetric)

per_node <- data.frame(
  node_id = setA_ids,
  category = category_all[setA_idx],
  n_chains = n_chains_full[setA_idx],
  mean_chain_length = round(mean_chain_length_full[setA_idx], 2),
  n_multiloop = obs$n_multiloop[setA_idx],
  P_first = round(obs$P_first[setA_idx], 4),
  P_middle = round(obs$P_middle[setA_idx], 4),
  P_last = round(obs$P_last[setA_idx], 4),
  p_chance_first_last = round(obs$p_chance[setA_idx], 4),
  p_chance_middle = round(p_chance_middle[setA_idx], 4),
  permutation_p_first = signif(
    per_node_perm_p(obs$P_first, primary_run$null$null_P_first, obs$p_chance, setA_idx), 4
  ),
  permutation_p_middle = signif(
    per_node_perm_p(obs$P_middle, primary_run$null$null_P_middle, p_chance_middle, setA_idx), 4
  ),
  permutation_p_last = signif(
    per_node_perm_p(obs$P_last, primary_run$null$null_P_last, obs$p_chance, setA_idx), 4
  ),
  stringsAsFactors = FALSE
)
per_node$BH_q_first <- signif(stats::p.adjust(per_node$permutation_p_first, method = "BH"), 4)
per_node$BH_q_middle <- signif(stats::p.adjust(per_node$permutation_p_middle, method = "BH"), 4)
per_node$BH_q_last <- signif(stats::p.adjust(per_node$permutation_p_last, method = "BH"), 4)
per_node$in_set_B <- per_node$node_id %in% setB_ids
per_node$high_strength <- per_node$node_id %in% high_strength_ids
per_node$high_betweenness <- per_node$node_id %in% high_betweenness_ids
per_node$high_turn_taking <- per_node$node_id %in% high_turntaking_ids
per_node <- per_node[order(-per_node$P_last), ]

write.csv(per_node, file.path(position_dir, "position_per_node.csv"), row.names = FALSE)
message("Saved ", file.path(position_dir, "position_per_node.csv"))

# ============================================================================
# 9. Full group-test table: first/middle/last x 3 roles x 2 node sets (18 rows)
# ============================================================================

group_tests <- rbind(primary_run$primary, primary_run$secondary, primary_run$tertiary)
write.csv(group_tests, file.path(position_dir, "position_group_tests.csv"), row.names = FALSE)
message("Saved ", file.path(position_dir, "position_group_tests.csv"))

# ============================================================================
# 10. Plain-text block: nodes 21, 1, 36, 26, 46 (for the response letter)
# ============================================================================

quote_lines <- c("", "Per-node values for the response letter:", "")
for (nd in QUOTE_NODES) {
  row <- per_node[per_node$node_id == nd, ]
  if (nrow(row) == 0) {
    quote_lines <- c(quote_lines, paste0("  node ", nd, ": not in Set A (n_multiloop < ", MIN_N_MULTILOOP, " or not in the significant network)."))
    next
  }
  quote_lines <- c(
    quote_lines,
    sprintf(
      paste0(
        "  node %d: P(first) = %.3f (p=%s), P(middle) = %.3f (p=%s), P(last) = %.3f (p=%s), n chains = %d"
      ),
      row$node_id,
      row$P_first, formatC(row$permutation_p_first, format = "f", digits = 4),
      row$P_middle, formatC(row$permutation_p_middle, format = "f", digits = 4),
      row$P_last, formatC(row$permutation_p_last, format = "f", digits = 4),
      row$n_chains
    )
  )
}
quote_lines <- c(quote_lines, "")
message(paste(quote_lines, collapse = "\n"))
diag_lines <- c(diag_lines, quote_lines)

# ============================================================================
# 11. Write diagnostics summary text file
# ============================================================================

summary_lines <- c(
  diag_lines,
  "Set A / Set B definitions",
  "==========================",
  paste0("Set A (n_multiloop >= ", MIN_N_MULTILOOP, ", in BH-significant network): n = ", length(setA_ids)),
  paste0("Set B raw pool (Figure 5c, high strength OR high betweenness): n = ", length(setB_raw_ids)),
  paste0(
    "Set B usable (after n_multiloop filter): n = ", length(setB_ids), " (",
    length(setB_dropped_ids), " dropped: ", paste(setB_dropped_ids, collapse = ", "), ")"
  ),
  "",
  "P(last) [terminal] group tests (3 roles x 2 node sets, BH within these 6):",
  utils::capture.output(print(primary_run$primary, row.names = FALSE)),
  "",
  "P(first) [introductory] group tests (3 roles x 2 node sets, BH within these 6):",
  utils::capture.output(print(primary_run$secondary, row.names = FALSE)),
  "",
  "P(middle) [central] group tests (3 roles x 2 node sets, BH within these 6):",
  utils::capture.output(print(primary_run$tertiary, row.names = FALSE)),
  "",
  "Diagnostic: excluding L=2 chains, P(last):",
  utils::capture.output(print(excl_L2_run$primary, row.names = FALSE)),
  "",
  "Diagnostic: excluding L=2 chains, P(first):",
  utils::capture.output(print(excl_L2_run$secondary, row.names = FALSE)),
  "",
  "Diagnostic: excluding L=2 chains, P(middle):",
  utils::capture.output(print(excl_L2_run$tertiary, row.names = FALSE))
)
if (!is.null(heterogeneous_run)) {
  summary_lines <- c(
    summary_lines, "",
    "Diagnostic: heterogeneous chains only (P(last)):",
    utils::capture.output(print(heterogeneous_run$primary, row.names = FALSE)),
    "",
    "Diagnostic: heterogeneous chains only (P(first)):",
    utils::capture.output(print(heterogeneous_run$secondary, row.names = FALSE)),
    "",
    "Diagnostic: heterogeneous chains only (P(middle)):",
    utils::capture.output(print(heterogeneous_run$tertiary, row.names = FALSE))
  )
}
# ============================================================================
# 12. P(middle), L>=3 chains only, promoted to a PRIMARY test (own printed
#     block + dedicated CSV) - same 6-test family (3 roles x 2 node sets),
#     same within-family BH correction excl_L2_run already applied to
#     $tertiary; this just re-surfaces it with primary-level prominence
#     rather than folded into the general L>=3 diagnostic dump.
# ============================================================================

message("\n================ PRIMARY: P(middle), L>=3 chains only ================")
message(
  length(unique(ml_L3$chain_id)), " chains (of ", n_chains_total,
  " total chain_length >= 2 chains); null_mean_high/rest = mean of that ",
  "group's own permutation-null distribution, for comparison against mean_high/rest."
)
print(excl_L2_run$tertiary, row.names = FALSE)
write.csv(
  excl_L2_run$tertiary, file.path(position_dir, "position_group_tests_L3_middle_primary.csv"),
  row.names = FALSE
)
message("Saved ", file.path(position_dir, "position_group_tests_L3_middle_primary.csv"))
message("=========================================================================\n")

# ============================================================================
# 13. Full membership of the high-strength / high-betweenness sets, with each
#     node's standardized strength (raw z-score), normalized betweenness,
#     and each measure's percentile rank across the full BH-significant
#     network (node_measures.csv, n = nrow(node_measures)).
# ============================================================================

percentile_rank <- function(x) 100 * rank(x, ties.method = "average") / length(x)
strength_percentile_all <- percentile_rank(strength_norm)
betweenness_percentile_all <- percentile_rank(betweenness_norm)

build_role_membership <- function(ids, order_by) {
  idx <- match(ids, node_measures$node_id)
  tbl <- data.frame(
    node_id = ids,
    category = node_measures$category[idx],
    z_score_strength = round(node_measures$z_score_incl_self[idx], 3),
    strength_norm = round(strength_norm[idx], 4),
    strength_percentile = round(strength_percentile_all[idx], 1),
    betweenness_raw = node_measures$betweenness[idx],
    betweenness_norm = round(betweenness_norm[idx], 4),
    betweenness_percentile = round(betweenness_percentile_all[idx], 1),
    stringsAsFactors = FALSE
  )
  tbl[order(-tbl[[order_by]]), ]
}

high_strength_membership <- build_role_membership(high_strength_ids, "strength_norm")
high_betweenness_membership <- build_role_membership(high_betweenness_ids, "betweenness_norm")

message("\n---- High-strength set membership (n=", length(high_strength_ids), "), sorted by strength_norm ----")
print(high_strength_membership, row.names = FALSE)
message("\n---- High-betweenness set membership (n=", length(high_betweenness_ids), "), sorted by betweenness_norm ----")
print(high_betweenness_membership, row.names = FALSE)

write.csv(
  rbind(
    cbind(role = "strength", high_strength_membership),
    cbind(role = "betweenness", high_betweenness_membership)
  ),
  file.path(position_dir, "position_role_membership_detail.csv"), row.names = FALSE
)
message("Saved ", file.path(position_dir, "position_role_membership_detail.csv"))

# ============================================================================
# 14. Top 5 nodes by normalized betweenness across the full BH-significant
#     network.
# ============================================================================

betweenness_rank_all <- data.frame(
  node_id = node_measures$node_id,
  category = node_measures$category,
  betweenness_raw = node_measures$betweenness,
  betweenness_norm = round(betweenness_norm, 4),
  betweenness_percentile = round(betweenness_percentile_all, 1),
  stringsAsFactors = FALSE
)
betweenness_rank_all <- betweenness_rank_all[order(-betweenness_rank_all$betweenness_norm), ]
top5_betweenness <- head(betweenness_rank_all, 5)
message("\n---- Top 5 nodes by normalized betweenness (full BH-significant network, n=", nrow(node_measures), ") ----")
print(top5_betweenness, row.names = FALSE)
write.csv(top5_betweenness, file.path(position_dir, "position_betweenness_top5.csv"), row.names = FALSE)
message("Saved ", file.path(position_dir, "position_betweenness_top5.csv"))

# ============================================================================
# 15. Per-node P(first)/P(middle)/P(last) + null/chance expectations for all
#     9 high-strength nodes (all already in Set A - confirmed above), so the
#     node-level contributions behind the SetA/SetB x strength group effects
#     are visible individually.
# ============================================================================

high_strength_breakdown <- per_node[per_node$high_strength, ]
high_strength_breakdown <- high_strength_breakdown[order(-high_strength_breakdown$P_last), ]
message("\n---- Per-node breakdown, high-strength nodes (n=", nrow(high_strength_breakdown), ") ----")
print(
  high_strength_breakdown[, c(
    "node_id", "category", "n_chains", "n_multiloop",
    "P_first", "p_chance_first_last", "permutation_p_first",
    "P_middle", "p_chance_middle", "permutation_p_middle",
    "P_last", "permutation_p_last"
  )],
  row.names = FALSE
)
write.csv(
  high_strength_breakdown, file.path(position_dir, "position_high_strength_breakdown.csv"),
  row.names = FALSE
)
message("Saved ", file.path(position_dir, "position_high_strength_breakdown.csv"))

# ============================================================================
# 16. L>=3-only: chains and nodes retained per group (from the enriched
#     n_retained_high/n_retained_rest columns already in excl_L2_run's
#     tables - a group member with zero occurrences among L>=3 chains has an
#     undefined P(position) and is silently excluded from that group's mean
#     unless flagged here).
# ============================================================================

l3_retained <- rbind(excl_L2_run$primary, excl_L2_run$secondary, excl_L2_run$tertiary)[, c(
  "node_set", "role", "position", "n_high", "n_rest", "n_retained_high", "n_retained_rest"
)]
message(
  "\n---- L>=3 only: ", length(unique(ml_L3$chain_id)), " chains (of ", n_chains_total,
  " total); nodes retained per group (n_high/n_rest = group size in the full-data node set; ",
  "n_retained_* = how many of those actually have >=1 occurrence among L>=3 chains) ----"
)
print(l3_retained, row.names = FALSE)

# ============================================================================
# Append all of the above to the diagnostics summary text file
# ============================================================================

extra_lines <- c(
  "", "PRIMARY: P(middle), L>=3 chains only", "======================================",
  paste0(length(unique(ml_L3$chain_id)), " chains (of ", n_chains_total, " total chain_length >= 2 chains)."),
  utils::capture.output(print(excl_L2_run$tertiary, row.names = FALSE)),
  "",
  paste0("High-strength set membership (n=", length(high_strength_ids), "):"),
  utils::capture.output(print(high_strength_membership, row.names = FALSE)),
  "",
  paste0("High-betweenness set membership (n=", length(high_betweenness_ids), "):"),
  utils::capture.output(print(high_betweenness_membership, row.names = FALSE)),
  "",
  paste0("Top 5 nodes by normalized betweenness (full network, n=", nrow(node_measures), "):"),
  utils::capture.output(print(top5_betweenness, row.names = FALSE)),
  "",
  "Per-node breakdown, high-strength nodes:",
  utils::capture.output(print(
    high_strength_breakdown[, c(
      "node_id", "category", "n_chains", "n_multiloop",
      "P_first", "p_chance_first_last", "permutation_p_first",
      "P_middle", "p_chance_middle", "permutation_p_middle",
      "P_last", "permutation_p_last"
    )],
    row.names = FALSE
  )),
  "",
  "L>=3 only: chains and nodes retained per group:",
  utils::capture.output(print(l3_retained, row.names = FALSE))
)
summary_lines <- c(summary_lines, extra_lines)

writeLines(summary_lines, file.path(position_dir, "position_diagnostics_summary.txt"))
message("Saved ", file.path(position_dir, "position_diagnostics_summary.txt"))
message("\nDone. All outputs in ", position_dir)

#!/usr/bin/env Rscript

# Sub-clusters five of the k = 9 top-level composition clusters (cut from
# the same substitution-cost/alignment tree cluster_multiloop_chains()
# builds from MULTILOOP_CLUSTER_EPSILON/CROSS_CATEGORY_COST/GAP_PENALTY/
# LINKAGE in R/config.R, just forced to k = 9 instead of that function's own
# silhouette-selected k) - cluster 6 (SW_Luna, 492 chains), 3 (SW_Neo, 266),
# 5 (SW_Nikita, 119), 2 (NSW_9, 106), 7 (SW_Nana, 95), the five biggest/
# cleanest ones - using the "hierarchical" substitution cost scheme
# (epsilon = 0.3, cross_category_cost = 1, same as the top level): within an
# already mostly-single-category cluster, this still treats a substitution
# between two sub-categories of the SAME main category as cheaper than one
# that crosses to a different main category - e.g. cluster 7 is 75% SW_Nana
# but 16% SW_Yosefa, so that distinction still has real minority structure
# to act on. (An alternative "uniform" scheme - no same-/different-main-
# category distinction at all - was tried and compared against this one;
# hierarchical was kept as more informative and is the only one here now.)
#
# Each cluster's own sub-k is chosen by MDL (mdl_k_sweep(),
# R/chain_clustering_utils.R - the PAM/k-medoids-re-optimized, BIC-penalized
# criterion used for the top-level k too), not by silhouette: k = 1 ("don't
# split further") is a native, well-defined point on the MDL curve, unlike
# silhouette (undefined at k = 1). Since cluster_multiloop_chains() itself
# always drives its own cluster assignment from silhouette, each
# sub-clustering here is done in two steps: build the tree or "scaffold"
# (whatever k_range satisfies that function's own internal requirement, its
# own best_k is unused), run mdl_k_sweep() against that scaffold, then cut
# the SAME tree again at the MDL-chosen k (force_k) - or leave it unsplit if
# MDL's own choice is k = 1.
#
# Stops at sub-cluster assignments + diagnostics - scripts/21 (the combined
# whistle/cluster/sub-cluster transition network) reads this script's own
# multiloop_subcluster_mdl_chain_assignments.csv output, so re-running this
# script and then scripts/21 keeps both in sync.

suppressPackageStartupMessages({
  library(cluster)
})

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
source(file.path(REPO_ROOT, "R", "chain_clustering_utils.R"))
source(file.path(REPO_ROOT, "R", "plot_io.R"))

ensure_plots_dir(REPO_ROOT, "multiloops")
multiloop_dir <- file.path(REPO_ROOT, OUTPUT_DIR, "multiloop_analysis")
dir.create(multiloop_dir, recursive = TRUE, showWarnings = FALSE)

TOP_K <- 9L
SUBCLUSTER_TARGET_CATEGORIES <- c("SW_Luna", "SW_Neo", "SW_Nikita", "NSW_9", "SW_Nana")
SUBCLUSTER_EPSILON <- 0.3
SUBCLUSTER_CROSS_CATEGORY_COST <- 1
# Both self-clip to each cluster's own valid range (cluster_multiloop_
# chains()/mdl_k_sweep() internals), so one shared upper bound - matching
# MULTILOOP_CLUSTER_K_RANGE's own top-level range - is enough for clusters
# this size (95-492 chains).
SUBCLUSTER_SIL_K_RANGE <- 2:40
SUBCLUSTER_MDL_K_RANGE <- 1:40

message("Loading whistles from ", DATA_PATH)
whistles_list <- load_whistles(DATA_PATH)

message(
  "Grouping whistles into multi-loop chains (IWI < ",
  MULTILOOP_IWI_THRESHOLD * 1000, "ms)..."
)
whistle_chains <- compute_multiloop_chains(
  whistles_list, threshold = MULTILOOP_IWI_THRESHOLD, list_names = LIST_NAMES
)

# -- top-level k = 9 clustering ----------------------------------------------

type_lookup <- unique(whistle_chains[, c("whistle_type_chr", "category")])
if (anyDuplicated(type_lookup$whistle_type_chr)) {
  stop("A whistle_type_chr maps to more than one main category", call. = FALSE)
}
sub_levels <- sort(type_lookup$whistle_type_chr)
category_of_sub <- setNames(type_lookup$category, type_lookup$whistle_type_chr)
top_cost_matrix <- build_chain_substitution_cost(
  sub_levels, category_of_sub, epsilon = MULTILOOP_CLUSTER_EPSILON,
  cross_category_cost = MULTILOOP_CLUSTER_CROSS_CATEGORY_COST
)
top_clustering <- cluster_multiloop_chains(
  whistle_chains, top_cost_matrix, gap_penalty = MULTILOOP_CLUSTER_GAP_PENALTY,
  linkage_method = MULTILOOP_CLUSTER_LINKAGE, k_range = MULTILOOP_CLUSTER_K_RANGE,
  force_k = TOP_K, progress_every = 0L
)
chain_level <- top_clustering$chain_level
message("Top-level clustering: k = ", TOP_K)

whistle_level_multiloop <- whistle_chains[whistle_chains$chain_length >= 2, ]
whistle_cluster <- chain_level$cluster[
  match(whistle_level_multiloop$chain_id, chain_level$chain_id)
]
clusters_present <- sort(unique(chain_level$cluster))
composition_counts <- table(
  category = factor(whistle_level_multiloop$category, levels = LIST_NAMES),
  cluster = factor(whistle_cluster, levels = clusters_present)
)
composition_pct <- prop.table(composition_counts, margin = 2) * 100
dominant_category <- setNames(
  LIST_NAMES[apply(composition_pct, 2, which.max)], as.character(clusters_present)
)

target_clusters <- clusters_present[
  dominant_category[as.character(clusters_present)] %in% SUBCLUSTER_TARGET_CATEGORIES
]
message(
  length(target_clusters), " target clusters: cluster ",
  paste0(
    target_clusters, " (", dominant_category[as.character(target_clusters)], ")",
    collapse = ", "
  )
)

# -- sub-cluster each target cluster (hierarchical cost, MDL k-selection) ---

subcluster_one <- function(cluster_id) {
  member_chain_ids <- chain_level$chain_id[chain_level$cluster == cluster_id]
  member_whistle_chains <- whistle_chains[
    whistle_chains$chain_length >= 2 & whistle_chains$chain_id %in% member_chain_ids,
  ]
  member_sub_levels <- sort(unique(member_whistle_chains$whistle_type_chr))
  scheme_cost_matrix <- build_chain_substitution_cost(
    member_sub_levels, category_of_sub[member_sub_levels],
    epsilon = SUBCLUSTER_EPSILON, cross_category_cost = SUBCLUSTER_CROSS_CATEGORY_COST
  )

  scaffold <- cluster_multiloop_chains(
    member_whistle_chains, scheme_cost_matrix,
    gap_penalty = MULTILOOP_CLUSTER_GAP_PENALTY, linkage_method = MULTILOOP_CLUSTER_LINKAGE,
    k_range = SUBCLUSTER_SIL_K_RANGE, progress_every = 0L
  )
  mdl_result <- mdl_k_sweep(
    scaffold, scheme_cost_matrix, k_range = SUBCLUSTER_MDL_K_RANGE, progress_every = 0L
  )
  mdl_best_k <- mdl_result$best_k

  if (mdl_best_k == 1L) {
    assignment <- setNames(rep(1L, length(member_chain_ids)), member_chain_ids)
    pattern_table <- scaffold$pattern_table
    pattern_table$cluster <- 1L
  } else {
    sub <- cluster_multiloop_chains(
      member_whistle_chains, scheme_cost_matrix,
      gap_penalty = MULTILOOP_CLUSTER_GAP_PENALTY, linkage_method = MULTILOOP_CLUSTER_LINKAGE,
      k_range = SUBCLUSTER_SIL_K_RANGE, progress_every = 0L, force_k = mdl_best_k
    )
    assignment <- setNames(sub$chain_level$cluster, sub$chain_level$chain_id)
    pattern_table <- sub$pattern_table
  }

  list(
    cluster_id = cluster_id,
    sub_k = mdl_best_k,
    total_bits = mdl_result$sweep$total_bits[mdl_result$sweep$k == mdl_best_k],
    assignment = assignment,
    pattern_table = pattern_table
  )
}

message("Sub-clustering each target cluster (hierarchical cost, MDL k-selection)...")
sub_results <- list()
for (cluster_id in target_clusters) {
  res <- subcluster_one(cluster_id)
  sub_results[[as.character(cluster_id)]] <- res
  message(
    "  cluster ", cluster_id, " (", dominant_category[as.character(cluster_id)],
    "): sub-k = ", res$sub_k, " (", round(res$total_bits, 1), " bits)"
  )
}

# -- combine into chain-level + comparison + per-sub-group tables -----------

target_chain_level <- chain_level[chain_level$cluster %in% target_clusters, ]
target_chain_level$dominant_category <- dominant_category[as.character(target_chain_level$cluster)]
target_chain_level$sub_cluster_hierarchical <- NA_integer_
for (cluster_id in target_clusters) {
  res <- sub_results[[as.character(cluster_id)]]
  member_idx <- target_chain_level$cluster == cluster_id
  target_chain_level$sub_cluster_hierarchical[member_idx] <-
    as.integer(res$assignment[target_chain_level$chain_id[member_idx]])
}
target_chain_level <- target_chain_level[
  order(target_chain_level$cluster, target_chain_level$sub_cluster_hierarchical), c(
    "chain_id", "cluster", "dominant_category", "chain_length", "ordered_pattern",
    "sub_cluster_hierarchical"
  )
]

cluster_summary <- do.call(rbind, lapply(target_clusters, function(cluster_id) {
  res <- sub_results[[as.character(cluster_id)]]
  data.frame(
    cluster = cluster_id,
    dominant_category = dominant_category[as.character(cluster_id)],
    n_chains = sum(chain_level$cluster == cluster_id),
    sub_k = res$sub_k,
    total_bits = round(res$total_bits, 1),
    stringsAsFactors = FALSE
  )
}))

# One row per (cluster, sub-cluster): size + dominant sub-category, derived
# entirely from that run's own pattern_table (weighted by n_chains) rather
# than re-touching whistle_chains.
build_sub_group_summary <- function(res, cluster_id) {
  pt <- res$pattern_table
  sub_ids <- sort(unique(pt$cluster))
  do.call(rbind, lapply(sub_ids, function(sc) {
    sub_pt <- pt[pt$cluster == sc, ]
    token_list <- strsplit(sub_pt$ordered_pattern, " -> ", fixed = TRUE)
    tokens <- unlist(token_list)
    weights <- rep(sub_pt$n_chains, lengths(token_list))
    tok_counts <- tapply(weights, tokens, sum)
    dominant_tok <- names(which.max(tok_counts))
    dominant_pct <- round(100 * max(tok_counts) / sum(tok_counts), 1)
    top_pattern_row <- sub_pt[which.max(sub_pt$n_chains), ]
    data.frame(
      cluster = cluster_id, sub_cluster = sc,
      n_chains = sum(sub_pt$n_chains), n_unique_patterns = nrow(sub_pt),
      mean_chain_length = round(
        sum(sub_pt$n_chains * sub_pt$pattern_length) / sum(sub_pt$n_chains), 2
      ),
      dominant_sub_category = dominant_tok, dominant_sub_category_pct = dominant_pct,
      top_pattern = top_pattern_row$ordered_pattern,
      top_pattern_n_chains = top_pattern_row$n_chains,
      stringsAsFactors = FALSE
    )
  }))
}
sub_group_summary <- do.call(rbind, lapply(target_clusters, function(cluster_id) {
  build_sub_group_summary(sub_results[[as.character(cluster_id)]], cluster_id)
}))
sub_group_summary <- sub_group_summary[
  order(sub_group_summary$cluster, -sub_group_summary$n_chains),
]

write.csv(
  target_chain_level,
  file.path(multiloop_dir, "multiloop_subcluster_mdl_chain_assignments.csv"),
  row.names = FALSE
)
write.csv(
  cluster_summary,
  file.path(multiloop_dir, "multiloop_subcluster_mdl_cluster_summary.csv"),
  row.names = FALSE
)
write.csv(
  sub_group_summary,
  file.path(multiloop_dir, "multiloop_subcluster_mdl_groups.csv"),
  row.names = FALSE
)
message("Saved sub-cluster assignment/summary/group CSVs under ", multiloop_dir)

# -- text summary --------------------------------------------------------------

summary_lines <- c(
  "Sub-clustering of five k = 9 top-level composition clusters (hierarchical",
  "substitution cost: epsilon = 0.3, cross_category_cost = 1; sub-k chosen by MDL):",
  "",
  paste0(
    "  cluster ", cluster_summary$cluster, " (", cluster_summary$dominant_category,
    ", ", cluster_summary$n_chains, " chains): sub-k = ", cluster_summary$sub_k,
    " (", cluster_summary$total_bits, " bits)"
  ),
  "",
  "per sub-group breakdown, largest first within each cluster:",
  paste0(
    "  cluster ", sub_group_summary$cluster, " sub-", sub_group_summary$sub_cluster, ": ",
    sub_group_summary$n_chains, " chains (", sub_group_summary$n_unique_patterns,
    " distinct compositions), mean length ", sub_group_summary$mean_chain_length,
    ", dominant sub-category ", sub_group_summary$dominant_sub_category,
    " (", sub_group_summary$dominant_sub_category_pct, "%)"
  ),
  "",
  "full breakdown: multiloop_subcluster_mdl_chain_assignments.csv (one row per chain),",
  "multiloop_subcluster_mdl_cluster_summary.csv (one row per cluster),",
  "multiloop_subcluster_mdl_groups.csv (one row per cluster x sub-cluster)"
)
writeLines(summary_lines, file.path(multiloop_dir, "multiloop_subcluster_mdl_summary.txt"))
message("Saved ", file.path(multiloop_dir, "multiloop_subcluster_mdl_summary.txt"))

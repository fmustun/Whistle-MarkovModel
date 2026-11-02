#!/usr/bin/env Rscript

# Manuscript Fig. S6 (Methods 4.11): stand-alone whistle / multi-loop
# composition network (a) and its betweenness vs. standardized strength (b).
#
# Combined whistle/multi-loop-cluster transition network at the finest
# granularity kept for this analysis: the five target clusters scripts/09
# sub-clustered (2 = NSW_9, 3 = SW_Neo, 5 = SW_Nikita, 6 = SW_Luna,
# 7 = SW_Nana) are each split into their own HIERARCHICAL sub-clusters
# (epsilon = 0.3/cross_category_cost = 1, MDL-selected k) instead of staying
# as one node apiece; the other four k = 9 clusters (1 = NSW_3, 4 = SW_Dana,
# 8 = SW_Shy, 9 = SW_Yosefa) - not sub-clustered - keep their single node.
# Reads scripts/09's own chain-level assignment CSV
# (multiloop_subcluster_mdl_chain_assignments.csv) rather than re-running
# the sub-clustering, since that file already has exactly the
# sub_cluster_hierarchical column needed here.
#
# Node universe: 85 standalone whistle sub-categories (script 01's
# granularity, unchanged) + 4 un-split top-level clusters + 25 hierarchical
# sub-clusters (7 + 3 + 3 + 10 + 2, from scripts/09) = 114 nodes total.

suppressPackageStartupMessages({
  library(igraph)
  library(ccber)
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
source(file.path(REPO_ROOT, "R", "statistical_utils.R"))
source(file.path(REPO_ROOT, "R", "plot_io.R"))

ensure_plots_dir(REPO_ROOT, "multiloops")
multiloop_dir <- file.path(REPO_ROOT, OUTPUT_DIR, "multiloop_analysis")
dir.create(multiloop_dir, recursive = TRUE, showWarnings = FALSE)

register_fcircle_shape()

TOP_K <- 9L
subcluster_csv <- file.path(multiloop_dir, "multiloop_subcluster_mdl_chain_assignments.csv")
if (!file.exists(subcluster_csv)) {
  stop("Run scripts/09_multiloop_composition_clustering.R first. Missing: ", subcluster_csv, call. = FALSE)
}

message("Loading whistles from ", DATA_PATH)
whistles_list <- load_whistles(DATA_PATH)

message(
  "Grouping whistles into multi-loop chains (IWI < ",
  MULTILOOP_IWI_THRESHOLD * 1000, "ms)..."
)
whistle_chains <- compute_multiloop_chains(
  whistles_list, threshold = MULTILOOP_IWI_THRESHOLD, list_names = LIST_NAMES
)

# -- top-level k = 9 clustering (same tree/cost as scripts/09) --------------

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
message("Clustering multi-loop chain compositions, forced to k = ", TOP_K, "...")
top_clustering <- cluster_multiloop_chains(
  whistle_chains, top_cost_matrix, gap_penalty = MULTILOOP_CLUSTER_GAP_PENALTY,
  linkage_method = MULTILOOP_CLUSTER_LINKAGE, k_range = MULTILOOP_CLUSTER_K_RANGE,
  force_k = TOP_K
)
chain_level <- top_clustering$chain_level

# -- graft scripts/09's hierarchical sub-clusters onto the target clusters --

subclusters <- read.csv(subcluster_csv, stringsAsFactors = FALSE)
sub_lookup <- setNames(subclusters$sub_cluster_hierarchical, subclusters$chain_id)
target_clusters <- sort(unique(subclusters$cluster))

is_target <- chain_level$cluster %in% target_clusters
chain_level$final_key <- ifelse(
  is_target,
  paste0("cluster_", chain_level$cluster, "_", sub_lookup[chain_level$chain_id]),
  paste0("cluster_", chain_level$cluster)
)
n_final_clusters <- length(unique(chain_level$final_key))
message(
  nrow(chain_level), " multi-loop chains -> ", length(unique(chain_level$cluster)),
  " top-level clusters, ", length(target_clusters), " of them sub-clustered -> ",
  n_final_clusters, " final multi-loop-cluster units"
)

# -- one row per chain of *any* length, singletons included -----------------

chains <- summarize_multiloop_chains(whistle_chains)
chains <- merge(
  chains, chain_level[, c("chain_id", "final_key")],
  by = "chain_id", all.x = TRUE, sort = FALSE
)
chains$is_multiloop <- chains$chain_length >= 2
chains$key <- ifelse(chains$is_multiloop, chains$final_key, chains$leading_whistle_type_chr)
chains$node_type <- ifelse(chains$is_multiloop, "multiloop_cluster", "whistle")

# -- one node per final multi-loop-cluster unit ------------------------------

ml_whistles <- whistle_chains[whistle_chains$chain_length >= 2, ]
ml_whistles <- merge(
  ml_whistles, chain_level[, c("chain_id", "final_key")], by = "chain_id"
)
final_keys_present <- sort(unique(chain_level$final_key))
cluster_category_counts <- table(
  final_key = factor(ml_whistles$final_key, levels = final_keys_present),
  category = factor(ml_whistles$category, levels = LIST_NAMES)
)
cluster_category_pct <- prop.table(cluster_category_counts, margin = 1) * 100
cluster_n_occurrences <- table(factor(chain_level$final_key, levels = final_keys_present))

multiloop_node_table <- data.frame(
  key = final_keys_present,
  n_categories = as.integer(rowSums(cluster_category_counts > 0)),
  main_category = LIST_NAMES[apply(cluster_category_pct, 1, which.max)],
  dominant_category_pct = round(apply(cluster_category_pct, 1, max), 1),
  n_occurrences = as.integer(cluster_n_occurrences),
  node_type = "multiloop_cluster",
  stringsAsFactors = FALSE
)

# -- one node per whistle sub-category that ever occurs on its own ----------

category_by_type_chr <- category_of_sub

whistle_freq <- as.data.frame(table(chains$leading_whistle_type_chr[!chains$is_multiloop]))
names(whistle_freq) <- c("key", "n_occurrences")
whistle_freq$key <- as.character(whistle_freq$key)
whistle_node_table <- data.frame(
  key = whistle_freq$key,
  n_categories = 1L,
  main_category = unname(category_by_type_chr[whistle_freq$key]),
  dominant_category_pct = 100,
  n_occurrences = whistle_freq$n_occurrences,
  node_type = "whistle",
  stringsAsFactors = FALSE
)

# -- combined node universe --------------------------------------------------

category_color_lookup <- setNames(LIST_COLORS, LIST_NAMES)

node_table <- rbind(
  whistle_node_table[, c(
    "key", "node_type", "n_categories", "main_category",
    "dominant_category_pct", "n_occurrences"
  )],
  multiloop_node_table[, c(
    "key", "node_type", "n_categories", "main_category",
    "dominant_category_pct", "n_occurrences"
  )]
)
node_table <- node_table[
  order(node_table$node_type, -node_table$n_occurrences, node_table$key),
]
node_table$node_id <- seq_len(nrow(node_table))
node_table$color <- category_color_lookup[node_table$main_category]

# Node number actually drawn in the plot: a whistle node keeps its
# whistle_type id from the original per-whistle-type network (script 01),
# so the same node reads as the same number across both plots; a multi-loop
# cluster/sub-cluster node (no equivalent in that original network) gets a
# number continuing right after the whistle_type universe, in the same
# largest-first order node_id above already uses. This is purely a display
# label - node_id (not display_id) remains the dense 1..n_nodes index every
# matrix/lookup below is built against.
type_id_map <- unique(whistles_list[, c("whistle_type_chr", "whistle_type")])
type_id_lookup <- setNames(type_id_map$whistle_type, type_id_map$whistle_type_chr)
max_type_id <- max(whistles_list$whistle_type)
is_whistle_node <- node_table$node_type == "whistle"
node_table$display_id <- NA_integer_
node_table$display_id[is_whistle_node] <-
  as.integer(type_id_lookup[node_table$key[is_whistle_node]])
node_table$display_id[!is_whistle_node] <-
  max_type_id + seq_len(sum(!is_whistle_node))

message(
  nrow(node_table), " nodes total: ",
  sum(node_table$node_type == "whistle"), " standalone whistle sub-categories, ",
  sum(node_table$node_type == "multiloop_cluster"), " multi-loop composition cluster/sub-cluster units"
)

node_csv <- file.path(multiloop_dir, "multiloop_cluster_network_subclusters_nodes.csv")
write.csv(node_table, node_csv, row.names = FALSE)
message("Saved ", node_csv)

# -- transitions between consecutive units -----------------------------------

ordered_chains <- chains[order(chains$recording, chains$start_time), ]
n <- nrow(ordered_chains)
same_recording_next <- ordered_chains$recording[-n] == ordered_chains$recording[-1]
gap <- ordered_chains$start_time[-1] - ordered_chains$end_time[-n]
within_window <- same_recording_next &
  gap >= TIME_WINDOW[1] & gap < TIME_WINDOW[2]

message(
  sum(within_window), " / ", sum(same_recording_next),
  " consecutive same-recording unit pairs fall within ",
  "TIME_WINDOW (", TIME_WINDOW[1], "-", TIME_WINDOW[2], "s)"
)

ordered_chains$node_key <- paste(ordered_chains$node_type, ordered_chains$key)
node_table$node_key <- paste(node_table$node_type, node_table$key)

from_node_key <- ordered_chains$node_key[-n][within_window]
to_node_key <- ordered_chains$node_key[-1][within_window]

node_key_to_id <- setNames(node_table$node_id, node_table$node_key)
n_nodes <- nrow(node_table)
from_node <- factor(node_key_to_id[from_node_key], levels = seq_len(n_nodes))
to_node <- factor(node_key_to_id[to_node_key], levels = seq_len(n_nodes))
transition_counts_matrix <- unclass(table(from_node, to_node))
dimnames(transition_counts_matrix) <- NULL

transition_probabilities_matrix <- CalcTransitionMatrix(transition_counts_matrix)

# -- significance test (same idea as script 01, see R/multiloop_utils.R) ----

unit_node_ids <- as.integer(node_key_to_id[ordered_chains$node_key])

m <- sum(transition_counts_matrix > 0)
iterations_cluster <- recommended_iterations(m, GRAPH_ALPHA, margin = 5)
message(
  "Compute combined whistle/multi-loop-cluster null model (", m, " observed edges, ",
  "GRAPH_ALPHA = ", GRAPH_ALPHA, " -> ", iterations_cluster,
  " shuffle-null iterations)..."
)
set.seed(SEED)
significance <- compute_multiloop_type_significance(
  unit_node_ids, within_window, n_nodes, transition_counts_matrix,
  iterations = iterations_cluster
)

family_mask <- transition_probabilities_matrix > 0
adjusted_p_value <- adjust_transition_p_values(
  significance$p_value_matrix,
  correction = MULTIPLE_TESTING_CORRECTION, family_mask = family_mask
)
selection <- select_transition_matrix(
  transition_probabilities_matrix, adjusted_p_value, alpha = GRAPH_ALPHA
)
selection_counts <- transition_selection_counts(selection$mask)
message(paste(
  names(selection_counts), selection_counts, sep = "=", collapse = "; "
))

edge_table <- which(transition_counts_matrix > 0, arr.ind = TRUE)
edge_table <- data.frame(
  from_node_id = edge_table[, 1],
  to_node_id = edge_table[, 2],
  from_node_type = node_table$node_type[edge_table[, 1]],
  to_node_type = node_table$node_type[edge_table[, 2]],
  from_key = node_table$key[edge_table[, 1]],
  to_key = node_table$key[edge_table[, 2]],
  n = transition_counts_matrix[edge_table],
  probability = transition_probabilities_matrix[edge_table],
  is_self_transition = edge_table[, 1] == edge_table[, 2],
  empirical_p = significance$p_value_matrix[edge_table],
  adjusted_p = adjusted_p_value[edge_table],
  significant = selection$mask[edge_table],
  stringsAsFactors = FALSE
)
if (MULTIPLE_TESTING_CORRECTION == "BH") {
  names(edge_table)[names(edge_table) == "adjusted_p"] <- "bh_q"
}
edge_table <- edge_table[order(-edge_table$n), ]
edge_csv <- file.path(multiloop_dir, "multiloop_cluster_network_subclusters_edges.csv")
write.csv(edge_table, edge_csv, row.names = FALSE)
message("Saved ", edge_csv, " (", nrow(edge_table), " edges, ",
        sum(selection$mask), " significant)")

# -- graph + plot ------------------------------------------------------------

gra_cluster <- graph_from_adjacency_matrix(
  selection$matrix, mode = "directed", weighted = TRUE
)
V(gra_cluster)$name <- as.character(node_table$node_id)
V(gra_cluster)$color <- node_table$color
V(gra_cluster)$node_type <- node_table$node_type
V(gra_cluster)$sub_category <- node_table$key
V(gra_cluster)$display_id <- as.character(node_table$display_id)

n_before <- vcount(gra_cluster)
gra_cluster_plot <- delete_vertices(gra_cluster, V(gra_cluster)[degree(gra_cluster) == 0])

node_id_plot <- as.integer(V(gra_cluster_plot)$name)
has_self_loop <- diag(selection$matrix)[node_id_plot] > 0
V(gra_cluster_plot)$vertex.frame.col <- ifelse(has_self_loop, "black", "gray")
gra_cluster_plot <- delete_edges(gra_cluster_plot, which(which_loop(gra_cluster_plot)))

# A node whose only significant transition was a self-loop is degree-0 once
# that loop is stripped above (self-loops are shown via a black border
# instead, never as a drawn edge) - drop those too, so a plotted node always
# has at least one edge to another node (matching what the caption's own
# "nodes with no significant transition are omitted" line claims).
n_after_degree_filter <- vcount(gra_cluster_plot)
gra_cluster_plot <- delete_vertices(
  gra_cluster_plot, V(gra_cluster_plot)[degree(gra_cluster_plot) == 0]
)
message(
  vcount(gra_cluster_plot), " / ", n_before, " nodes plotted (",
  n_before - n_after_degree_filter, " with no significant transition at all, ",
  n_after_degree_filter - vcount(gra_cluster_plot),
  " with only a significant self-transition, both omitted)"
)

is_cluster_node <- V(gra_cluster_plot)$node_type == "multiloop_cluster"
V(gra_cluster_plot)$label_color <- ifelse(is_cluster_node, "black", "darkblue")
V(gra_cluster_plot)$label_font <- ifelse(is_cluster_node, 2, 1)

CLUSTER_NODE_SIZE <- 8
WHISTLE_NODE_SIZE <- 5
CLUSTER_NODE_FRAME_WIDTH <- 4.5
WHISTLE_NODE_FRAME_WIDTH <- 2.5
V(gra_cluster_plot)$vertex_size <- ifelse(
  is_cluster_node, CLUSTER_NODE_SIZE, WHISTLE_NODE_SIZE
)
V(gra_cluster_plot)$vertex_frame_width <- ifelse(
  is_cluster_node, CLUSTER_NODE_FRAME_WIDTH, WHISTLE_NODE_FRAME_WIDTH
)

message("Plotting combined whistle / multi-loop-cluster network (with hierarchical sub-clusters)...")
with_pdf_plot(
  plot_path("multiloop_cluster_network_subclusters", REPO_ROOT, subdir = "multiloops"),
  width = 12,
  height = 12,
  plot_markov_graph(
    gra_cluster_plot,
    seed = PLOT_SEED,
    show_legend = FALSE,
    vertex_label = V(gra_cluster_plot)$display_id,
    vertex_label_color = V(gra_cluster_plot)$label_color,
    vertex_label_font = V(gra_cluster_plot)$label_font,
    vertex_size = V(gra_cluster_plot)$vertex_size,
    vertex_frame_width = V(gra_cluster_plot)$vertex_frame_width,
    caption = c(
      paste0(
        "each node is either one whistle sub-category produced on its own (as in script 01), one un-split k = 9 multi-loop cluster (NSW_3/SW_Dana/SW_Shy/SW_Yosefa), ",
        "or one hierarchical (epsilon = ", MULTILOOP_CLUSTER_EPSILON, ", cross-category = ",
        MULTILOOP_CLUSTER_CROSS_CATEGORY_COST, ") sub-cluster of the other five k = 9 clusters (Luna/Neo/Nikita/NSW_9/Nana - MDL-selected sub-k, see scripts/09)"
      ),
      "edge width = transition probability between nodes; node color = dominant main category (exact for a whistle node, plurality for a cluster/sub-cluster node - see multiloop_cluster_network_subclusters_nodes.csv for dominant_category_pct)",
      "multi-loop cluster/sub-cluster nodes are drawn larger with a thicker border and a bold black node number; standalone whistle sub-categories are smaller/thin-bordered with a navy node number",
      paste0(
        "node numbers: standalone whistle sub-categories keep their whistle_type id from the original network (script 01, 1-", max_type_id,
        "); multi-loop cluster/sub-cluster nodes are numbered ", max_type_id + 1, "-", max_type_id + sum(!is_whistle_node), " (see multiloop_cluster_network_subclusters_nodes.csv's display_id column)"
      ),
      "black border = node has a significant self-transition (BH-adjusted, not drawn as an edge)",
      paste0(
        "edges require the next unit within TIME_WINDOW (",
        TIME_WINDOW[1], "-", TIME_WINDOW[2], "s) and ",
        MULTIPLE_TESTING_CORRECTION, "-adjusted p < ", GRAPH_ALPHA,
        " vs. a shuffle-null model (", iterations_cluster, " iterations); ",
        "nodes with no significant transition are omitted"
      )
    )
  )
)

# Fig. S6b: betweenness vs. standardized strength (z-score against the same
# shuffle-null model already run above for edge significance - significance$strength_z_score, self-transitions included,
# see the comment on compute_permutation_transition_significance() in
# R/multiloop_utils.R) instead of plot_betweenness_strength_scatter()'s own
# strength(gra, mode = "all"). Mirrors scripts/04's markov_betweenness_vs_
# strength_zscore.png for the whistle-level network, computed here for the
# combined whistle/multi-loop-cluster network instead.
message("Plotting betweenness vs. standardized strength (z-score)...")
plot_betweenness_strength_scatter(
  gra_cluster_plot,
  plot_path("multiloop_cluster_betweenness_vs_strength_zscore_subclusters", REPO_ROOT, subdir = "multiloops"),
  width = 12,
  vertex_label = V(gra_cluster_plot)$display_id,
  vertex_label_color = V(gra_cluster_plot)$label_color,
  vertex_label_font = V(gra_cluster_plot)$label_font,
  vertex_cex = ifelse(is_cluster_node, 2.6, 1.6),
  vertex_lwd = ifelse(is_cluster_node, 2.5, 1),
  strength_values = significance$strength_z_score[as.integer(V(gra_cluster_plot)$name)],
  strength_label = "standardized strength (z-score vs. shuffle null)",
  # Points at the observed betweenness/z-score min or max would otherwise
  # sit exactly on the panel border (xaxs/yaxs = "i", the default axis_
  # padding = 0) and read as half-clipped - matches scripts/04's own
  # markov_betweenness_vs_strength_zscore.png, which uses the same padding
  # for the same reason.
  axis_padding = 0.05,
  # No caption, bigger text - matches scripts/04's own
  # markov_betweenness_vs_strength_zscore.png styling.
  show_caption = FALSE,
  text_cex = 1.5
)

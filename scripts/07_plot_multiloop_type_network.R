#!/usr/bin/env Rscript

# Builds a Markov-style transition network the same way script 01 does
# (transition counts -> row-normalized probabilities -> directed graph,
# tested against a "shuffle" null model, BH-corrected, only the
# significant edges plotted), but at a different unit: each node is a
# "multi-loop type" - a distinct multi-loop composition (the same
# sub-categories, in any order/position; see compute_multiloop_patterns()'s
# unordered_pattern) - rather than an individual whistle sub-category. The
# significance test itself is compute_multiloop_type_significance()
# (R/multiloop_utils.R): the same empirical-p/BH-correction idea as script
# 01's MarkovWhistle(), but computed only at the observed-edge cells that
# can ever end up significant, since this network has far more nodes (630)
# relative to observed edges (~300) than the whistle-level one does.

suppressPackageStartupMessages({
  library(igraph)
  library(ccber)
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
source(file.path(REPO_ROOT, "R", "statistical_utils.R"))
source(file.path(REPO_ROOT, "R", "plot_io.R"))

ensure_plots_dir(REPO_ROOT, "multiloops")
multiloop_dir <- file.path(REPO_ROOT, OUTPUT_DIR, "multiloop_analysis")
dir.create(multiloop_dir, recursive = TRUE, showWarnings = FALSE)

register_fcircle_shape()

message("Loading whistles from ", DATA_PATH)
whistles_list <- load_whistles(DATA_PATH)

message(
  "Grouping whistles into multi-loop chains (IWI < ",
  MULTILOOP_IWI_THRESHOLD * 1000, "ms)..."
)
whistle_chains <- compute_multiloop_chains(
  whistles_list, threshold = MULTILOOP_IWI_THRESHOLD, list_names = LIST_NAMES
)
chains <- summarize_multiloop_chains(whistle_chains)
multiloop_chains <- chains[chains$chain_length >= 2, ]
patterns <- compute_multiloop_patterns(whistle_chains)
multiloop_chains <- merge(
  multiloop_chains, patterns[, c("chain_id", "unordered_pattern")],
  by = "chain_id"
)

# -- one node per distinct multi-loop composition --------------------------

ml_whistles <- whistle_chains[whistle_chains$chain_length >= 2, ]
ml_whistles <- merge(
  ml_whistles, patterns[, c("chain_id", "unordered_pattern")], by = "chain_id"
)
n_categories <- aggregate(
  category ~ unordered_pattern, data = ml_whistles,
  FUN = function(x) length(unique(x))
)
names(n_categories)[2] <- "n_categories"
main_category <- aggregate(
  category ~ unordered_pattern, data = ml_whistles,
  FUN = function(x) if (length(unique(x)) == 1) unique(x) else NA_character_
)
names(main_category)[2] <- "main_category"
pattern_freq <- as.data.frame(table(multiloop_chains$unordered_pattern))
names(pattern_freq) <- c("unordered_pattern", "n_chains")

node_table <- merge(n_categories, main_category, by = "unordered_pattern")
node_table <- merge(node_table, pattern_freq, by = "unordered_pattern")
node_table <- node_table[order(-node_table$n_chains, node_table$unordered_pattern), ]
node_table$node_id <- seq_len(nrow(node_table))

# Single main category -> that category's usual color; a composition that
# mixes multiple main categories (e.g. an SW_Neo whistle immediately
# followed by an SW_Luna one within the same multi-loop chain) gets a
# distinct fixed color instead, since it doesn't belong to just one. Plain
# white was tried first but read as too close to NSW_9's very pale
# lavender (#e6e6fa) at a glance; a mid gray has no near neighbor in
# LIST_COLORS (all saturated hues or pale pastels) and reads unambiguously
# as "not a real category".
MIXED_CATEGORY_COLOR <- "gray50"
category_color_lookup <- setNames(LIST_COLORS, LIST_NAMES)
node_table$color <- ifelse(
  node_table$n_categories == 1,
  category_color_lookup[node_table$main_category],
  MIXED_CATEGORY_COLOR
)

message(
  nrow(node_table), " distinct multi-loop types (",
  sum(node_table$n_categories == 1), " single-category, ",
  sum(node_table$n_categories > 1), " mixed-category)"
)

node_csv <- file.path(multiloop_dir, "multiloop_type_nodes.csv")
write.csv(node_table, node_csv, row.names = FALSE)
message("Saved ", node_csv)

# -- transitions between consecutive multi-loop chains ----------------------

# Mirrors script 01's MarkovWhistle(): a transition is only counted when the
# next multi-loop chain in the same recording starts within TIME_WINDOW of
# the previous chain's end (same gap definition as MarkovWhistle's
# time_interval). Unlike individual whistles, multi-loop chains never
# overlap in time (their own gap-based definition requires a >=
# MULTILOOP_IWI_THRESHOLD silence before a new chain starts), so there's no
# need for MarkovWhistle's overlap/skip-to-next-candidate handling here -
# just a straightforward window check.
ordered_chains <- multiloop_chains[
  order(multiloop_chains$recording, multiloop_chains$start_time),
]
n <- nrow(ordered_chains)
same_recording_next <- ordered_chains$recording[-n] == ordered_chains$recording[-1]
gap <- ordered_chains$start_time[-1] - ordered_chains$end_time[-n]
within_window <- same_recording_next &
  gap >= TIME_WINDOW[1] & gap < TIME_WINDOW[2]

message(
  sum(within_window), " / ", sum(same_recording_next),
  " consecutive same-recording multi-loop chain pairs fall within ",
  "TIME_WINDOW (", TIME_WINDOW[1], "-", TIME_WINDOW[2], "s)"
)

from_pattern <- ordered_chains$unordered_pattern[-n][within_window]
to_pattern <- ordered_chains$unordered_pattern[-1][within_window]

pattern_to_node <- setNames(node_table$node_id, node_table$unordered_pattern)
n_nodes <- nrow(node_table)
from_node <- factor(pattern_to_node[from_pattern], levels = seq_len(n_nodes))
to_node <- factor(pattern_to_node[to_pattern], levels = seq_len(n_nodes))
transition_counts_matrix <- unclass(table(from_node, to_node))
dimnames(transition_counts_matrix) <- NULL

transition_probabilities_matrix <- CalcTransitionMatrix(transition_counts_matrix)

# -- significance test (same idea as script 01, see R/multiloop_utils.R) ----

# Chain-level counterpart of MarkovWhistle()'s NULL_MODEL = "shuffle": a
# global permutation of which composition label sits on which multi-loop
# chain, keeping every chain's timing - and therefore `within_window` -
# fixed. Built once here since it's shared by every null iteration inside
# compute_multiloop_type_significance().
chain_node_ids <- as.integer(pattern_to_node[ordered_chains$unordered_pattern])

m <- sum(transition_counts_matrix > 0)
iterations_type <- recommended_iterations(m, GRAPH_ALPHA, margin = 5)
message(
  "Compute multi-loop type null model (", m, " observed edges, ",
  "GRAPH_ALPHA = ", GRAPH_ALPHA, " -> ", iterations_type,
  " shuffle-null iterations)..."
)
set.seed(SEED)
significance <- compute_multiloop_type_significance(
  chain_node_ids, within_window, n_nodes, transition_counts_matrix,
  iterations = iterations_type
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
  from_pattern = node_table$unordered_pattern[edge_table[, 1]],
  to_pattern = node_table$unordered_pattern[edge_table[, 2]],
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
edge_csv <- file.path(multiloop_dir, "multiloop_type_edges.csv")
write.csv(edge_table, edge_csv, row.names = FALSE)
message("Saved ", edge_csv, " (", nrow(edge_table), " edges, ",
        sum(selection$mask), " significant)")

# -- graph + plot ------------------------------------------------------------

# Built from the significant matrix only (non-significant/never-observed
# cells are already zero there) - same as script 01's gra1, which is built
# from selection$matrix rather than the full transition probability matrix.
gra_type <- graph_from_adjacency_matrix(
  selection$matrix, mode = "directed", weighted = TRUE
)
V(gra_type)$name <- as.character(node_table$node_id)
V(gra_type)$color <- node_table$color

# A multi-loop type with no significant transition in either direction -
# including a significant self-transition, since degree() counts loop
# edges - has nothing to show on this plot; drop it here (only from the
# plot; multiloop_type_nodes.csv above still lists every composition).
# Mirrors compute_graph()'s isolate filter (R/graph_utils.R), which also
# runs before self-loop edges are stripped, so a type whose only
# significant transition is a self-loop stays in the plot as an isolated,
# black-bordered dot rather than being dropped.
n_before <- vcount(gra_type)
gra_type_plot <- delete_vertices(gra_type, V(gra_type)[degree(gra_type) == 0])
message(
  vcount(gra_type_plot), " / ", n_before,
  " multi-loop types have at least one significant transition (plotted)"
)

# Black border = this multi-loop type has a *significant* self-transition
# (a nonzero diagonal entry in selection$matrix) - the same "black node
# border, no drawn loop edge" convention script 01 uses for significant
# Markov self-transitions (see compute_graph_from_selection(),
# R/graph_utils.R). Assigned after the isolate filter (matching diag() up
# by name, since node ordering is unaffected by delete_vertices()) and the
# self-loop edges themselves are dropped last so they never render as arcs.
node_id_plot <- as.integer(V(gra_type_plot)$name)
has_self_loop <- diag(selection$matrix)[node_id_plot] > 0
V(gra_type_plot)$vertex.frame.col <- ifelse(has_self_loop, "black", "gray")
gra_type_plot <- delete_edges(gra_type_plot, which(which_loop(gra_type_plot)))

message("Plotting multi-loop type network...")
with_pdf_plot(
  plot_path("multiloop_type_network", REPO_ROOT, subdir = "multiloops"),
  width = 9,
  height = 9,
  plot_markov_graph(
    gra_type_plot,
    seed = PLOT_SEED,
    show_legend = FALSE,
    caption = c(
      "each node = one multi-loop composition (same sub-categories, any order/position); edge width = transition probability between multi-loop types",
      "node color = single main category (as in the whistle-type network); gray = composition mixes more than one main category",
      "black border = type has a significant self-transition (BH-adjusted, not drawn as an edge)",
      paste0(
        "edges require the next multi-loop chain within TIME_WINDOW (",
        TIME_WINDOW[1], "-", TIME_WINDOW[2], "s) and ",
        MULTIPLE_TESTING_CORRECTION, "-adjusted p < ", GRAPH_ALPHA,
        " vs. a shuffle-null model (", iterations_type, " iterations); ",
        "types with no significant transition are omitted"
      )
    )
  )
)

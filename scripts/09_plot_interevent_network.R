#!/usr/bin/env Rscript

# The complement of script 08: same node universe (every whistle
# sub-category MarkovWhistle() would consider) and the exact same
# BH-corrected empirical-shuffle significance test, but the transition
# counts feeding that test are restricted to Markov transition instances
# where the two whistles do NOT belong to the same multi-loop chain
# (compute_markov_transition_multiloop_overlap()'s within_multiloop_chain
# flag, negated) - "inter-event" transitions, i.e. every Markov transition
# except the multi-loop-internal ones script 08 isolates. Together, script
# 08's edge list and this one's exactly partition every Markov transition
# instance MarkovWhistle() would ever count.

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
# Same node universe MarkovWhistle() dimensions its matrices to - every
# whistle sub-category in the data, not just the subset script 01's own
# significance test happens to retain (gra1 is already a pruned subset).
whistle_number <- max(whistles_list$whistle_type)

# -- transition instances excluding multi-loop-internal pairs ---------------

message(
  "Grouping whistles into multi-loop chains (IWI < ",
  MULTILOOP_IWI_THRESHOLD * 1000, "ms) and finding Markov transition ",
  "instances (TIME_WINDOW ", TIME_WINDOW[1], "-", TIME_WINDOW[2], "s)..."
)
overlap <- compute_markov_transition_multiloop_overlap(
  whistles_list, time_window = TIME_WINDOW,
  threshold = MULTILOOP_IWI_THRESHOLD, list_names = LIST_NAMES
)
filtered <- overlap[!overlap$within_multiloop_chain, ]
message(
  nrow(filtered), " / ", nrow(overlap),
  " Markov transition instances are inter-event (not within a multi-loop chain) (",
  round(100 * nrow(filtered) / nrow(overlap), 1), "%)"
)

from_type <- factor(filtered$from_whistle_type, levels = seq_len(whistle_number))
to_type <- factor(filtered$to_whistle_type, levels = seq_len(whistle_number))
transition_counts_matrix <- unclass(table(from_type, to_type))
dimnames(transition_counts_matrix) <- NULL
transition_probabilities_matrix <- CalcTransitionMatrix(transition_counts_matrix)

# -- significance test (same idea as script 01, see R/multiloop_utils.R) ----

# NULL_MODEL = "shuffle" permutes whistle_type globally across every
# whistle in whistles_list, keeping timing fixed - and since multi-loop
# chain membership is purely a function of gaps between start/end times,
# within_multiloop_chain (and therefore which transition instances are
# even candidates here) is unaffected by that permutation. So, exactly as
# in scripts 07/08, one null iteration only needs to permute the label
# vector and re-look-up the fixed set of transition-instance row pairs -
# `filtered$from_row_id`/`to_row_id` - rather than rebuilding chain
# structure per iteration the way MarkovWhistle() rebuilds transitions.
labels <- whistles_list$whistle_type
m <- sum(transition_counts_matrix > 0)
iterations_ie <- recommended_iterations(m, GRAPH_ALPHA, margin = 5)
message(
  "Compute inter-event null model (", m, " observed edges, ",
  "GRAPH_ALPHA = ", GRAPH_ALPHA, " -> ", iterations_ie,
  " shuffle-null iterations)..."
)
set.seed(SEED)
significance <- compute_permutation_transition_significance(
  labels, filtered$from_row_id, filtered$to_row_id, whistle_number,
  transition_counts_matrix, iterations = iterations_ie
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

whistle_chr_map <- tapply(
  whistles_list$whistle_type_chr, whistles_list$whistle_type, function(x) x[1]
)
edge_table <- which(transition_counts_matrix > 0, arr.ind = TRUE)
edge_table <- data.frame(
  from_whistle_type = edge_table[, 1],
  to_whistle_type = edge_table[, 2],
  from_sub_category = unname(whistle_chr_map[as.character(edge_table[, 1])]),
  to_sub_category = unname(whistle_chr_map[as.character(edge_table[, 2])]),
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
edge_csv <- file.path(multiloop_dir, "interevent_network_edges.csv")
write.csv(edge_table, edge_csv, row.names = FALSE)
message("Saved ", edge_csv, " (", nrow(edge_table), " edges, ",
        sum(selection$mask), " significant)")

# -- graph + plot -------------------------------------------------------------

# Same node/edge construction path as script 01's gra1 (and script 08's) -
# compute_graph() assigns category colors and drops zero-degree nodes,
# compute_graph_from_selection() marks significant self-transitions with a
# black border and drops the self-loop edges themselves - just fed
# transition_matrix = selection$matrix built from the inter-event-filtered
# counts above instead of MarkovWhistle()'s full transition matrix.
whistle_occurrence <- table(whistles_list$whistle_type)
selected_graphs <- compute_graph_from_selection(
  transition_matrix = selection$matrix,
  selection_values = adjusted_p_value,
  Markov_Model = list(whistle_occurrence = whistle_occurrence),
  Whistles_List = whistles_list,
  list_names = LIST_NAMES,
  list_colors = LIST_COLORS,
  edge_attribute_name = if (MULTIPLE_TESTING_CORRECTION == "BH") {
    "bh_q"
  } else {
    "empirical_p"
  },
  raw_p_values = significance$p_value_matrix
)
gra_interevent <- selected_graphs$no_loops

message(
  vcount(gra_interevent), " / ", whistle_number,
  " whistle sub-categories have at least one significant",
  " inter-event transition (plotted)"
)

message("Plotting inter-event significant network...")
node_table <- plot_network_by_component(
  gra_interevent,
  plot_path("interevent_significant_network", REPO_ROOT, subdir = "multiloops"),
  caption = c(
    "same node universe and significance test as the full Markov network (script 01)",
    paste0(
      "transition counts restricted to instances where the two whistles are NOT in the same multi-loop chain (IWI < ",
      MULTILOOP_IWI_THRESHOLD * 1000, "ms) - ", nrow(filtered), " / ", nrow(overlap), " (",
      round(100 * nrow(filtered) / nrow(overlap), 1), "%) of all Markov transition instances qualify"
    ),
    paste0(
      MULTIPLE_TESTING_CORRECTION, "-adjusted p < ", GRAPH_ALPHA,
      " vs. a shuffle-null model (", iterations_ie, " iterations); ",
      "black border = significant self-transition (not drawn as an edge)"
    )
  )
)
node_csv <- file.path(multiloop_dir, "interevent_network_nodes.csv")
write.csv(node_table, node_csv, row.names = FALSE)
message("Saved ", node_csv)

message("Plotting betweenness vs. strength...")
plot_betweenness_strength_scatter(
  gra_interevent,
  plot_path("interevent_betweenness_vs_strength", REPO_ROOT, subdir = "multiloops")
)

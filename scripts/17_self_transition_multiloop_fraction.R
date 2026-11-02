#!/usr/bin/env Rscript

# For each node's self-transitions (the excess connectivity that
# scripts/15_node_strength_frequency_null.R found for the most frequent
# whistle types, e.g. nodes 1 and 21), what fraction actually fall inside a
# tight-IWI multi-loop chain (MULTILOOP_IWI_THRESHOLD, R/config.R) rather
# than being independent repeats separated by a longer gap? A self-loop
# that's really just two elements of the same multi-loop bout is a
# segmentation artifact of chain-vs-Markov-transition counting, not
# evidence the animal chose to repeat the whistle type in a fresh instance.

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
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
source(file.path(REPO_ROOT, "R", "graph_utils.R"))
source(file.path(REPO_ROOT, "R", "multiloop_utils.R"))
source(file.path(REPO_ROOT, "R", "plot_io.R"))

ensure_plots_dir(REPO_ROOT, "network_measures")

rds_path <- file.path(REPO_ROOT, OUTPUT_DIR, "markov_model.rds")
if (!file.exists(rds_path)) {
  stop("Run scripts/01_run_markov_model.R first. Missing: ", rds_path, call. = FALSE)
}

saved <- readRDS(rds_path)
MarkovModel <- saved$MarkovModel
whistles_list <- saved$whistles_list
list_names <- saved$config$list_names

message(
  "Grouping whistles into multi-loop chains (IWI < ",
  MULTILOOP_IWI_THRESHOLD * 1000, "ms) and replicating Markov transition ",
  "pairing to flag each transition instance..."
)
whistle_chains <- compute_multiloop_chains(
  whistles_list, threshold = MULTILOOP_IWI_THRESHOLD, list_names = list_names
)
transition_overlap <- compute_markov_transition_multiloop_overlap(
  whistles_list,
  time_window = MarkovModel$time_window,
  threshold = MULTILOOP_IWI_THRESHOLD,
  list_names = list_names
)

self_transitions <- transition_overlap[
  transition_overlap$from_whistle_type == transition_overlap$to_whistle_type,
]

# compute_markov_transition_multiloop_overlap() replicates MarkovWhistle()'s
# own pairing logic exactly, so this instance-level self-transition count
# must equal transition_counts_matrix's diagonal - a direct check that this
# script is slicing the same self-loops the strength/z-score figures in
# scripts/15 are built from, not a different or approximate count.
n_nodes <- nrow(MarkovModel$transition_counts_matrix)
node_ids <- as.character(seq_len(n_nodes))
observed_self_counts <- as.integer(
  table(factor(self_transitions$from_whistle_type, levels = seq_len(n_nodes)))
)
stopifnot(identical(observed_self_counts, as.integer(diag(MarkovModel$transition_counts_matrix))))

n_self_multiloop <- as.integer(
  table(factor(
    self_transitions$from_whistle_type[self_transitions$within_multiloop_chain],
    levels = seq_len(n_nodes)
  ))
)

# Whistle-level (not transition-level) fraction: of every occurrence of a
# node's whistle type, how many sit in a singleton chain (chain_length == 1,
# i.e. no neighbor within MULTILOOP_IWI_THRESHOLD on either side) rather
# than any multi-loop chain at all, self-repeat or not. This is a broader
# question than fraction_self_multiloop above: it counts every whistle, not
# just the ones that happen to be a same-type Markov self-transition.
n_whistles_total <- as.integer(
  table(factor(whistle_chains$whistle_type, levels = seq_len(n_nodes)))
)
stopifnot(identical(
  n_whistles_total, as.integer(MarkovModel$whistle_occurrence[node_ids])
))
n_whistles_not_multiloop <- as.integer(
  table(factor(
    whistle_chains$whistle_type[whistle_chains$chain_length == 1],
    levels = seq_len(n_nodes)
  ))
)
n_whistles_multiloop <- n_whistles_total - n_whistles_not_multiloop

# Strength "including self" is the same total-degree quantity plotted in
# strength_vs_frequency_null_incl_self.png (scripts/15), recomputed directly
# here (a single rowSums+colSums, no null loop) so this script only depends
# on markov_model.rds, not on scripts/15 having already been run.
strength_incl_self <- rowSums(MarkovModel$transition_counts_matrix) +
  colSums(MarkovModel$transition_counts_matrix)

node_graph <- make_empty_graph(n = n_nodes, directed = FALSE)
V(node_graph)$name <- node_ids
seq_category <- build_seq_category(whistles_list, list_names)
node_graph <- assign_vertex_category_color(
  node_graph, seq_category, list_names, LIST_COLORS
)
whistle_chr_map <- tapply(
  whistles_list$whistle_type_chr, whistles_list$whistle_type,
  function(x) x[1]
)

df <- data.frame(
  node_id = node_ids,
  sub_category = unname(whistle_chr_map[node_ids]),
  category = V(node_graph)$category,
  color = V(node_graph)$color,
  occurrences = as.numeric(MarkovModel$whistle_occurrence[node_ids]),
  strength_incl_self = strength_incl_self,
  n_self_transitions = observed_self_counts,
  n_self_transitions_multiloop = n_self_multiloop,
  fraction_self_multiloop = n_self_multiloop / observed_self_counts,
  n_whistles_multiloop = n_whistles_multiloop,
  fraction_whistles_multiloop = n_whistles_multiloop / n_whistles_total
)

csv_path <- file.path(REPO_ROOT, OUTPUT_DIR, "self_transition_multiloop_fraction.csv")
write.csv(df, csv_path, row.names = FALSE, quote = FALSE)
message("Saved ", csv_path)

# Both bar charts share the same styling and differ only in node order and
# which column/labels they draw, so one helper draws each rather than
# duplicating the ggplot call. `node_order` is the caller-supplied character
# vector of node_id, left to right.
plot_node_fraction_bar <- function(df, y_col, node_order, x_label, title, y_label, caption, name) {
  df$node_id <- factor(df$node_id, levels = node_order)
  save_ggplot(
    ggplot(df) +
      aes(x = node_id, y = .data[[y_col]], fill = category) +
      geom_col(width = 0.75, color = "black") +
      scale_fill_manual(breaks = unique(df$category), values = unique(df$color)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45)) +
      labs(title = title, x = x_label, y = y_label, caption = caption),
    name,
    REPO_ROOT,
    width = 12,
    subdir = "network_measures"
  )
}

# Nodes with zero self-transitions have no fraction to plot at all (0/0);
# excluded here rather than shown as a misleading 0. The whistle-level
# fraction below has no such gap - every node has occurrences > 0 - so it
# plots all 86 nodes.
df_self <- df[df$n_self_transitions > 0, ]
n_excluded <- nrow(df) - nrow(df_self)
message(n_excluded, " / ", nrow(df), " nodes have zero self-transitions and are excluded from that plot")

plot_node_fraction_bar(
  df_self, "fraction_self_multiloop",
  node_order = df_self$node_id[order(-df_self$strength_incl_self)],
  x_label = "node (ordered by strength, self-transitions included; left = highest)",
  title = "Fraction of each node's self-transitions inside a multi-loop chain",
  y_label = "fraction of self-transitions with IWI < 0.25s",
  caption = paste0(n_excluded, " node(s) with zero self-transitions omitted"),
  name = "self_transition_multiloop_fraction"
)

# Category order requested for this figure specifically - not LIST_NAMES'
# own order (which has SW_Yosefa before SW_Nikita/SW_Nana) - reusing
# order_nodes_by_category() (R/graph_utils.R) with that order instead of the
# global one; within each category, nodes still sort by strength descending.
category_plot_order <- c(
  "SW_Neo", "SW_Luna", "SW_Nikita", "SW_Nana", "SW_Yosefa", "SW_Dana",
  "SW_Shy", "NSW_3", "NSW_6", "NSW_9"
)
category_ordered <- order_nodes_by_category(
  data.frame(nodes = df$node_id, category = df$category, strength_incl_self = df$strength_incl_self),
  category_plot_order, "strength_incl_self"
)

plot_node_fraction_bar(
  df, "fraction_whistles_multiloop",
  node_order = levels(category_ordered$nodes),
  x_label = "node (grouped by category: Neo, Luna, Nikita, Nana, Yosefa, Dana, Shy, NSW_3, NSW_6, NSW_9)",
  title = "Fraction of each node's whistles involved in a multi-loop chain",
  y_label = "fraction of occurrences with a neighbor within IWI < 0.25s",
  caption = "all 86 nodes; a whistle counts as 'multi-loop' if its chain length is >= 2",
  name = "whistles_in_multiloop_fraction"
)

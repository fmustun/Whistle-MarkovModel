#!/usr/bin/env Rscript

# Builds a Markov-style transition network the same way script 01 does
# (transition counts -> row-normalized probabilities -> directed graph,
# tested against a "shuffle" null model, BH-corrected, only the
# significant edges plotted), but at a mixed unit: a whistle produced on
# its own (not part of a multi-loop chain) is its own node - the same
# sub-category node script 01 would give it - while a whistle produced
# inside a multi-loop chain (IWI < MULTILOOP_IWI_THRESHOLD) is collapsed
# into a single node for its chain's *composition*: the same distinct
# sub-categories present, regardless of order, position, or how many times
# each repeats - e.g. "A + A + B" and "A + B" collapse to the same node,
# since both chains ever contain exactly {A, B} (see
# compute_multiloop_patterns()'s composition_pattern). Two whistles of the
# same sub-category therefore end up in different nodes whenever one
# occurs standalone and the other occurs inside a loop - that split is the
# point: this network asks "does this whistle behave differently when
# produced in a loop?" rather than pooling every occurrence into one node
# as script 01 does.
#
# Both node kinds sit in one combined node universe and one combined set of
# transition counts, so a single significance test (the same empirical-p/
# BH-correction idea as script 01's MarkovWhistle(), computed only at the
# observed-edge cells that can ever end up significant - see
# compute_multiloop_type_significance(), R/multiloop_utils.R) covers
# transitions between any pair of nodes: whistle -> whistle, whistle ->
# multi-loop type, multi-loop type -> whistle, or multi-loop type ->
# multi-loop type.

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
# One row per chain of *any* length - unlike the original version of this
# script, singleton chains (chain_length == 1, i.e. a whistle with no
# neighbor close enough to loop with) are kept: each becomes its own
# "unit" below, exactly like a standalone whistle in script 01's network.
chains <- summarize_multiloop_chains(whistle_chains)
patterns <- compute_multiloop_patterns(whistle_chains)
chains <- merge(
  chains, patterns[, c("chain_id", "composition_pattern")],
  by = "chain_id", all.x = TRUE
)
# composition_pattern is NA for singleton chains (compute_multiloop_
# patterns() only computes it for chain_length >= 2); `key` picks the right
# identity for either kind of unit, and `node_type` records which kind it
# is so a whistle-type key can never collide with a composition_pattern
# string that happens to spell the same thing.
chains$is_multiloop <- chains$chain_length >= 2
chains$key <- ifelse(
  chains$is_multiloop, chains$composition_pattern, chains$leading_whistle_type_chr
)
chains$node_type <- ifelse(chains$is_multiloop, "multiloop", "whistle")

# -- one node per distinct multi-loop composition ---------------------------

# Grouped by composition_pattern (the *set* of distinct sub-categories in
# the chain), not unordered_pattern (which also keeps repeat counts) - so
# e.g. a 2xA+1xB chain and a 1xA+1xB chain land in the same node here, per
# the composition_pattern definition above.
ml_whistles <- whistle_chains[whistle_chains$chain_length >= 2, ]
ml_whistles <- merge(
  ml_whistles, patterns[, c("chain_id", "composition_pattern")], by = "chain_id"
)
n_categories <- aggregate(
  category ~ composition_pattern, data = ml_whistles,
  FUN = function(x) length(unique(x))
)
names(n_categories)[2] <- "n_categories"
main_category <- aggregate(
  category ~ composition_pattern, data = ml_whistles,
  FUN = function(x) if (length(unique(x)) == 1) unique(x) else NA_character_
)
names(main_category)[2] <- "main_category"
pattern_freq <- as.data.frame(table(chains$composition_pattern[chains$is_multiloop]))
names(pattern_freq) <- c("composition_pattern", "n_occurrences")

multiloop_node_table <- merge(n_categories, main_category, by = "composition_pattern")
multiloop_node_table <- merge(multiloop_node_table, pattern_freq, by = "composition_pattern")
multiloop_node_table$key <- multiloop_node_table$composition_pattern
multiloop_node_table$node_type <- "multiloop"

# -- one node per whistle sub-category that ever occurs on its own ----------

# Same sub-category granularity script 01's network uses (whistle_type_chr),
# restricted to occurrences that are *not* part of a multi-loop chain -
# singleton chains are exactly those occurrences (see compute_multiloop_
# chains(): a chain of length 1 is a whistle with no close-enough
# neighbor). category comes along for free via compute_multiloop_chains()'s
# own type-to-category lookup, so no re-derivation is needed here.
type_chr_category <- unique(whistle_chains[, c("whistle_type_chr", "category")])
category_by_type_chr <- setNames(
  type_chr_category$category, type_chr_category$whistle_type_chr
)

whistle_freq <- as.data.frame(table(
  chains$leading_whistle_type_chr[!chains$is_multiloop]
))
names(whistle_freq) <- c("key", "n_occurrences")
whistle_freq$key <- as.character(whistle_freq$key)
whistle_node_table <- data.frame(
  key = whistle_freq$key,
  n_categories = 1L,
  main_category = unname(category_by_type_chr[whistle_freq$key]),
  n_occurrences = whistle_freq$n_occurrences,
  node_type = "whistle",
  stringsAsFactors = FALSE
)

# -- combined node universe --------------------------------------------------

# Single main category -> that category's usual color; a composition that
# mixes multiple main categories (e.g. an SW_Neo whistle immediately
# followed by an SW_Luna one within the same multi-loop chain) gets a
# distinct fixed color instead, since it doesn't belong to just one. Plain
# white was tried first but read as too close to NSW_9's very pale
# lavender (#e6e6fa) at a glance; a mid gray has no near neighbor in
# LIST_COLORS (all saturated hues or pale pastels) and reads unambiguously
# as "not a real category". A whistle node's n_categories is always 1 (it
# is a single sub-category by construction), so this same rule colors it
# with its own category, exactly like script 01.
MIXED_CATEGORY_COLOR <- "gray50"
category_color_lookup <- setNames(LIST_COLORS, LIST_NAMES)

node_table <- rbind(
  whistle_node_table[, c("key", "node_type", "n_categories", "main_category", "n_occurrences")],
  multiloop_node_table[, c("key", "node_type", "n_categories", "main_category", "n_occurrences")]
)
node_table <- node_table[
  order(node_table$node_type, -node_table$n_occurrences, node_table$key),
]
node_table$node_id <- seq_len(nrow(node_table))
node_table$color <- ifelse(
  node_table$n_categories == 1,
  category_color_lookup[node_table$main_category],
  MIXED_CATEGORY_COLOR
)

message(
  nrow(node_table), " nodes total: ",
  sum(node_table$node_type == "whistle"), " standalone whistle sub-categories, ",
  sum(node_table$node_type == "multiloop"), " distinct multi-loop compositions (",
  sum(node_table$node_type == "multiloop" & node_table$n_categories == 1),
  " single-category, ",
  sum(node_table$node_type == "multiloop" & node_table$n_categories > 1),
  " mixed-category)"
)

node_csv <- file.path(multiloop_dir, "multiloop_type_nodes.csv")
write.csv(node_table, node_csv, row.names = FALSE)
message("Saved ", node_csv)

# -- transitions between consecutive units -----------------------------------

# Mirrors script 01's MarkovWhistle(): a transition is only counted when the
# next unit in the same recording starts within TIME_WINDOW of the
# previous unit's end (same gap definition as MarkovWhistle's
# time_interval). Units - whether a single standalone whistle or a whole
# multi-loop chain - never overlap in time (a new chain only starts after
# a >= MULTILOOP_IWI_THRESHOLD silence, and a standalone whistle is by
# definition a chain of its own), so there's no need for MarkovWhistle's
# overlap/skip-to-next-candidate handling here - just a straightforward
# window check, same as the multi-loop-only version of this script.
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

# node_type + key together are the unique identity of a unit's node (a
# whistle-type key and a composition_pattern key live in disjoint namespaces
# in practice, but the composite avoids relying on that).
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

# Chain-level counterpart of MarkovWhistle()'s NULL_MODEL = "shuffle": a
# global permutation of which node sits on which unit, keeping every
# unit's timing - and therefore `within_window` - fixed. Built once here
# since it's shared by every null iteration inside
# compute_multiloop_type_significance() (name notwithstanding, that
# function is generic over what a "chain" node id represents - see its
# doc comment in R/multiloop_utils.R).
unit_node_ids <- as.integer(node_key_to_id[ordered_chains$node_key])

m <- sum(transition_counts_matrix > 0)
iterations_type <- recommended_iterations(m, GRAPH_ALPHA, margin = 5)
message(
  "Compute combined whistle/multi-loop-type null model (", m, " observed edges, ",
  "GRAPH_ALPHA = ", GRAPH_ALPHA, " -> ", iterations_type,
  " shuffle-null iterations)..."
)
set.seed(SEED)
significance <- compute_multiloop_type_significance(
  unit_node_ids, within_window, n_nodes, transition_counts_matrix,
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
V(gra_type)$node_type <- node_table$node_type
V(gra_type)$sub_category <- node_table$key

# A node with no significant transition in either direction - including a
# significant self-transition, since degree() counts loop edges - has
# nothing to show on this plot; drop it here (only from the plot;
# multiloop_type_nodes.csv above still lists every node). Mirrors
# compute_graph()'s isolate filter (R/graph_utils.R), which also runs
# before self-loop edges are stripped, so a node whose only significant
# transition is a self-loop stays in the plot as an isolated,
# black-bordered dot rather than being dropped.
n_before <- vcount(gra_type)
gra_type_plot <- delete_vertices(gra_type, V(gra_type)[degree(gra_type) == 0])
message(
  vcount(gra_type_plot), " / ", n_before,
  " nodes have at least one significant transition (plotted)"
)

# Black border = this node has a *significant* self-transition (a nonzero
# diagonal entry in selection$matrix) - the same "black node border, no
# drawn loop edge" convention script 01 uses for significant Markov
# self-transitions (see compute_graph_from_selection(), R/graph_utils.R).
# Assigned after the isolate filter (matching diag() up by name, since node
# ordering is unaffected by delete_vertices()) and the self-loop edges
# themselves are dropped last so they never render as arcs.
node_id_plot <- as.integer(V(gra_type_plot)$name)
has_self_loop <- diag(selection$matrix)[node_id_plot] > 0
V(gra_type_plot)$vertex.frame.col <- ifelse(has_self_loop, "black", "gray")
gra_type_plot <- delete_edges(gra_type_plot, which(which_loop(gra_type_plot)))

message("Plotting combined whistle / multi-loop-type network...")
with_pdf_plot(
  plot_path("multiloop_type_network", REPO_ROOT, subdir = "multiloops"),
  width = 9,
  height = 9,
  plot_markov_graph(
    gra_type_plot,
    seed = PLOT_SEED,
    show_legend = FALSE,
    caption = c(
      "each node is either one whistle sub-category produced on its own (as in script 01), or one multi-loop composition (the same distinct sub-categories present, any order/position/repeat count, produced as a loop)",
      "edge width = transition probability between nodes; node color = single main category (as in the whistle-type network); gray = a multi-loop composition mixing more than one main category",
      "black border = node has a significant self-transition (BH-adjusted, not drawn as an edge)",
      paste0(
        "edges require the next unit within TIME_WINDOW (",
        TIME_WINDOW[1], "-", TIME_WINDOW[2], "s) and ",
        MULTIPLE_TESTING_CORRECTION, "-adjusted p < ", GRAPH_ALPHA,
        " vs. a shuffle-null model (", iterations_type, " iterations); ",
        "nodes with no significant transition are omitted"
      )
    )
  )
)

#!/usr/bin/env Rscript

# Draws a multi-loop-only network: same nodes/layout as script 01's
# significant-transition network, but edges are replaced entirely by
# adjacent *different*-type pairs observed within multi-loop chains
# (IWI < MULTILOOP_IWI_THRESHOLD) - not the Markov model's wide-TIME_WINDOW
# transitions. A node gets a black border if it has a same-type multi-loop
# repeat (e.g. A -> A within a chain); gray otherwise.

suppressPackageStartupMessages({
  library(igraph)
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
source(file.path(REPO_ROOT, "R", "plot_io.R"))

ensure_plots_dir(REPO_ROOT, "multiloops")
multiloop_dir <- file.path(REPO_ROOT, OUTPUT_DIR, "multiloop_analysis")
dir.create(multiloop_dir, recursive = TRUE, showWarnings = FALSE)

model_rds <- file.path(REPO_ROOT, OUTPUT_DIR, "markov_model.rds")
if (!file.exists(model_rds)) {
  stop("Run scripts/01_run_markov_model.R first. Missing: ", model_rds, call. = FALSE)
}

saved <- readRDS(model_rds)
gra1 <- saved$gra1

register_fcircle_shape()

# Computed directly from gra1 - the same function/seed/graph script 01 uses
# for markov_significant_network.pdf (via plot_markov_graph(gra1, seed =
# PLOT_SEED) with coordinates = NULL) - rather than a saved CSV, since
# layout_with_fr is a force-directed layout: it depends on the edge set
# it's computed on, so recomputing it later on a graph with *different*
# edges (as gra_multiloop below will have) would silently drift from gra1's
# original layout even with the same seed. Computing it here, once, on the
# untouched gra1, guarantees identical node positions.
coordinates <- compute_markov_layout(gra1, PLOT_SEED)

message("Loading whistles from ", DATA_PATH)
whistles_list <- load_whistles(DATA_PATH)

message(
  "Grouping whistles into multi-loop chains (IWI < ",
  MULTILOOP_IWI_THRESHOLD * 1000, "ms)..."
)
whistle_chains <- compute_multiloop_chains(
  whistles_list, threshold = MULTILOOP_IWI_THRESHOLD, list_names = LIST_NAMES
)

edge_counts <- compute_multiloop_edge_counts(whistle_chains)
self_repeats <- compute_self_repeat_counts(whistle_chains)

node_ids <- as.integer(V(gra1)$name)
sub_category_lookup <- setNames(V(gra1)$sub_category, node_ids)

# Nodes are pinned to script 01's significant-network layout, so a
# multi-loop edge can only be drawn if both its whistle types are part of
# that node set; edges to/from a type outside it are dropped (reported
# below) rather than silently extending the network with new nodes.
in_network <- edge_counts$from_whistle_type %in% node_ids &
  edge_counts$to_whistle_type %in% node_ids
n_dropped <- sum(edge_counts$n[!in_network])
edge_counts <- edge_counts[in_network, ]

message(
  nrow(edge_counts), " distinct multi-loop edges (", sum(edge_counts$n),
  " pair instances); ", n_dropped,
  " pair instances dropped (whistle type outside the significant network)"
)

edge_counts_csv <- file.path(multiloop_dir, "multiloop_edge_counts.csv")
write.csv(
  data.frame(
    from_whistle_type = edge_counts$from_whistle_type,
    from_sub_category = sub_category_lookup[as.character(edge_counts$from_whistle_type)],
    to_whistle_type = edge_counts$to_whistle_type,
    to_sub_category = sub_category_lookup[as.character(edge_counts$to_whistle_type)],
    n = edge_counts$n,
    stringsAsFactors = FALSE
  ),
  edge_counts_csv,
  row.names = FALSE
)
message("Saved ", edge_counts_csv)

has_self_multiloop <- node_ids %in% self_repeats$whistle_type
node_table <- data.frame(
  whistle_type = node_ids,
  sub_category = V(gra1)$sub_category,
  n_self_repeats = self_repeats$n_self_repeats[
    match(node_ids, self_repeats$whistle_type)
  ],
  stringsAsFactors = FALSE
)
node_table$n_self_repeats[is.na(node_table$n_self_repeats)] <- 0
node_table$has_self_multiloop <- has_self_multiloop
node_csv <- file.path(multiloop_dir, "node_multiloop_participation.csv")
write.csv(node_table, node_csv, row.names = FALSE)
message("Saved ", node_csv)
message(sum(has_self_multiloop), " / ", length(node_ids), " nodes have a same-type multi-loop repeat")

# Start from gra1 to keep its vertex set/attributes (position, category
# color, occurrences) and just replace the edges entirely with the
# multi-loop-internal ones computed above.
gra_multiloop <- delete_edges(gra1, E(gra1))
gra_multiloop <- add_edges(
  gra_multiloop,
  edges = as.vector(rbind(
    as.character(edge_counts$from_whistle_type),
    as.character(edge_counts$to_whistle_type)
  ))
)
E(gra_multiloop)$weight <- edge_counts$n
V(gra_multiloop)$vertex.frame.col <- ifelse(has_self_multiloop, "black", "gray")

message("Plotting multi-loop-only network...")
with_pdf_plot(
  plot_path("multiloop_network_overlay", REPO_ROOT, subdir = "multiloops"),
  width = 9,
  height = 9,
  plot_markov_graph(
    gra_multiloop,
    list_names = LIST_NAMES,
    seed = PLOT_SEED,
    coordinates = coordinates,
    show_legend = FALSE,
    caption = c(
      "nodes at the same position as the significant Markov network; edges are multi-loop-internal transitions only",
      "edge width = count of that transition within multi-loop chains; black border = node has a same-type multi-loop repeat"
    )
  )
)

#!/usr/bin/env Rscript

# A third whistle-level network alongside script 01's full one: the same
# node universe (every whistle sub-category MarkovWhistle() would consider,
# not just script 01's already-significant subset) and the exact same
# BH-corrected empirical-shuffle significance test, but the transition
# counts feeding that test are restricted to Markov transition instances
# where both whistles belong to the same multi-loop chain (IWI <
# MULTILOOP_IWI_THRESHOLD; compute_markov_transition_multiloop_overlap()'s
# within_multiloop_chain flag) - i.e. "does the significant Markov
# structure look different if you only look at multi-loop-internal
# transitions?" This is a rebuild of script 01's whole pipeline on a
# filtered edge list, not an overlay on top of script 01's existing
# network (that's scripts/06, which uses raw multi-loop counts with no
# significance test and pins node positions to script 01's layout).

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

# -- transition instances restricted to multi-loop-internal pairs -----------

message(
  "Grouping whistles into multi-loop chains (IWI < ",
  MULTILOOP_IWI_THRESHOLD * 1000, "ms) and finding Markov transition ",
  "instances (TIME_WINDOW ", TIME_WINDOW[1], "-", TIME_WINDOW[2], "s)..."
)
overlap <- compute_markov_transition_multiloop_overlap(
  whistles_list, time_window = TIME_WINDOW,
  threshold = MULTILOOP_IWI_THRESHOLD, list_names = LIST_NAMES
)
filtered <- overlap[overlap$within_multiloop_chain, ]
message(
  nrow(filtered), " / ", nrow(overlap),
  " Markov transition instances are within a multi-loop chain (",
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
# in script 07, one null iteration only needs to permute the label vector
# and re-look-up the fixed set of transition-instance row pairs -
# `filtered$from_row_id`/`to_row_id` - rather than rebuilding chain
# structure per iteration the way MarkovWhistle() rebuilds transitions.
labels <- whistles_list$whistle_type
m <- sum(transition_counts_matrix > 0)
iterations_ml <- recommended_iterations(m, GRAPH_ALPHA, margin = 5)
message(
  "Compute multi-loop-only null model (", m, " observed edges, ",
  "GRAPH_ALPHA = ", GRAPH_ALPHA, " -> ", iterations_ml,
  " shuffle-null iterations)..."
)
set.seed(SEED)
significance <- compute_permutation_transition_significance(
  labels, filtered$from_row_id, filtered$to_row_id, whistle_number,
  transition_counts_matrix, iterations = iterations_ml
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
edge_csv <- file.path(multiloop_dir, "multiloop_only_network_edges.csv")
write.csv(edge_table, edge_csv, row.names = FALSE)
message("Saved ", edge_csv, " (", nrow(edge_table), " edges, ",
        sum(selection$mask), " significant)")

# -- graph + plot -------------------------------------------------------------

# Same node/edge construction path as script 01's gra1 - compute_graph()
# assigns category colors and drops zero-degree nodes, compute_graph_from_
# selection() marks significant self-transitions with a black border and
# drops the self-loop edges themselves - just fed transition_matrix =
# selection$matrix built from the multi-loop-filtered counts above instead
# of MarkovWhistle()'s full transition matrix.
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
gra_multiloop_only <- selected_graphs$no_loops

message(
  vcount(gra_multiloop_only), " / ", whistle_number,
  " whistle sub-categories have at least one significant",
  " multi-loop-only transition (plotted)"
)

# A network built from only the significant, multi-loop-filtered edges
# often isn't fully connected - e.g. a pair of sub-categories that only
# ever transition to each other, nowhere near the rest of the network in
# the data. Forcing compute_markov_layout()'s force-directed layout to
# place every component in one shared coordinate space stretches the whole
# plot out to make room for these outliers, at the cost of legibility for
# the (much larger) main component. Each connected component instead gets
# its own layout and its own panel - the largest as the main panel, every
# other component stacked below it, smallest last - so nothing is dropped
# but nothing distorts the main network's layout either.
components_info <- components(gra_multiloop_only, mode = "weak")
component_order <- order(-components_info$csize)
n_components <- components_info$no
message(
  n_components, " connected component(s) in the plotted network (sizes: ",
  paste(components_info$csize[component_order], collapse = ", "), ")"
)

component_graphs <- lapply(component_order, function(component_id) {
  induced_subgraph(
    gra_multiloop_only,
    V(gra_multiloop_only)[components_info$membership == component_id]
  )
})
# compute_markov_layout()'s default half_width (1.3) is tuned for the
# main, ~80-node network; applied unscaled to a 2-node component it just
# spreads two points across the same wide box, mostly empty space. Scale
# it down for smaller components (relative to the main component's size),
# floored so a 1-2 node component doesn't collapse to a single point.
main_size <- components_info$csize[component_order[1]]
component_coordinates <- lapply(component_graphs, function(g_i) {
  half_width <- max(0.3, 1.3 * sqrt(vcount(g_i) / main_size))
  compute_markov_layout(g_i, seed = PLOT_SEED, half_width = half_width)
})
# plot_markov_graph() defaults weight_reference_range to that call's own
# range(E(gra)$weight), which is degenerate (min == max, so
# rescale_from_reference() errors) whenever a component has only a single
# edge - true of every small disconnected component here. Pass one shared
# range, computed across all edges before the split, to every panel so
# edge widths stay on the same, non-degenerate scale.
weight_reference_range <- range(E(gra_multiloop_only)$weight)

node_table <- do.call(rbind, lapply(seq_along(component_graphs), function(i) {
  g_i <- component_graphs[[i]]
  coords_i <- component_coordinates[[i]]
  data.frame(
    node_id = as.integer(V(g_i)$name),
    sub_category = V(g_i)$sub_category,
    category = V(g_i)$category,
    occurrences = V(g_i)$occurrences,
    has_significant_self_transition = V(g_i)$vertex.frame.col == "black",
    component = i,
    component_size = vcount(g_i),
    x = coords_i[, "x"],
    y = coords_i[, "y"],
    stringsAsFactors = FALSE
  )
}))
node_table <- node_table[order(node_table$component, node_table$node_id), ]
node_csv <- file.path(multiloop_dir, "multiloop_only_network_nodes.csv")
write.csv(node_table, node_csv, row.names = FALSE)
message("Saved ", node_csv)

message("Plotting multi-loop-only significant network...")
main_caption <- c(
  "same node universe and significance test as the full Markov network (script 01)",
  paste0(
    "transition counts restricted to instances where both whistles are in the same multi-loop chain (IWI < ",
    MULTILOOP_IWI_THRESHOLD * 1000, "ms) - ", nrow(filtered), " / ", nrow(overlap), " (",
    round(100 * nrow(filtered) / nrow(overlap), 1), "%) of all Markov transition instances qualify"
  ),
  paste0(
    MULTIPLE_TESTING_CORRECTION, "-adjusted p < ", GRAPH_ALPHA,
    " vs. a shuffle-null model (", iterations_ml, " iterations); ",
    "black border = significant self-transition (not drawn as an edge)"
  )
)
if (n_components > 1) {
  main_caption <- c(main_caption, paste0(
    "largest connected component shown above (", vcount(component_graphs[[1]]),
    " / ", vcount(gra_multiloop_only), " plotted nodes); ",
    n_components - 1L, " smaller disconnected component(s) shown below"
  ))
}

with_pdf_plot(
  plot_path("multiloop_only_significant_network", REPO_ROOT, subdir = "multiloops"),
  width = 9,
  height = 9,
  {
    if (n_components == 1) {
      plot_markov_graph(
        component_graphs[[1]],
        list_names = LIST_NAMES,
        seed = PLOT_SEED,
        coordinates = component_coordinates[[1]],
        weight_reference_range = weight_reference_range,
        show_legend = FALSE,
        caption = main_caption
      )
    } else {
      n_other <- n_components - 1L
      # Row 1 = the main component, spanning every column; row 2 = one
      # column per smaller component, sized well below the main row.
      layout(
        matrix(
          c(rep(1L, n_other), seq_len(n_other) + 1L),
          nrow = 2, byrow = TRUE
        ),
        heights = c(5, 1.4)
      )
      plot_markov_graph(
        component_graphs[[1]],
        list_names = LIST_NAMES,
        seed = PLOT_SEED,
        coordinates = component_coordinates[[1]],
        weight_reference_range = weight_reference_range,
        show_legend = FALSE,
        caption = main_caption
      )
      for (i in seq_len(n_other)) {
        g_i <- component_graphs[[i + 1L]]
        plot_markov_graph(
          g_i,
          list_names = LIST_NAMES,
          seed = PLOT_SEED,
          coordinates = component_coordinates[[i + 1L]],
          weight_reference_range = weight_reference_range,
          show_legend = FALSE,
          # A `main=` title uses plot.igraph()'s default title cex, sized
          # for the large main panel - comically oversized on these small
          # ones. caption's smaller, fixed cex (see plot_markov_graph())
          # reads correctly regardless of panel size.
          caption = paste0(
            "component ", i + 1L, " (", vcount(g_i), " node",
            if (vcount(g_i) > 1) "s" else "", ")"
          )
        )
      }
    }
  }
)

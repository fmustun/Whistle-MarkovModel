#!/usr/bin/env Rscript

# Overlays the most recurrent multi-loop sub-category co-occurrences (see
# scripts/11_multiloop_type_cooccurrence.R) directly on script 01's
# significant Markov transition network - same nodes, layout, and edges as
# markov_significant_network.pdf, unchanged - as two separate plots:
#  - pairs: an orange line for every sub-category pair that co-occurred
#    (order/position ignored) in >= PAIR_MIN_CHAINS multi-loop chains,
#    connecting its 2 nodes' existing positions in the significant network.
#  - triplets: an orange triangle (translucent orange fill, solid orange
#    border) for every sub-category triplet that co-occurred in
#    >= TRIPLET_MIN_CHAINS multi-loop chains, connecting its 3 nodes.
# A pair/triplet with a member sub-category outside the significant network
# (no significant Markov transition at all) can't be placed and is skipped -
# same handling as script 06's multi-loop-only overlay.

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

# User-specified recurrence thresholds (not derived): a triplet needs >= 3
# shared multi-loop chains, a pair needs >= 15, before it's drawn.
TRIPLET_MIN_CHAINS <- 3
PAIR_MIN_CHAINS <- 10
# Same color for both plots, per request; the triplet plot additionally
# fills each triangle with a translucent version of it.
OVERLAY_COLOR <- "darkorange"

model_rds <- file.path(REPO_ROOT, OUTPUT_DIR, "markov_model.rds")
if (!file.exists(model_rds)) {
  stop("Run scripts/01_run_markov_model.R first. Missing: ", model_rds, call. = FALSE)
}

saved <- readRDS(model_rds)
gra1 <- saved$gra1

register_fcircle_shape()

# Same function/seed/graph as script 01 uses for markov_significant_network.pdf
# (and script 06 reuses the same way) - guarantees identical node positions.
coordinates <- compute_markov_layout(gra1, PLOT_SEED)
valid_node_names <- V(gra1)$name

message("Loading whistles from ", DATA_PATH)
whistles_list <- load_whistles(DATA_PATH)

message(
  "Grouping whistles into multi-loop chains (IWI < ",
  MULTILOOP_IWI_THRESHOLD * 1000, "ms)..."
)
whistle_chains <- compute_multiloop_chains(
  whistles_list, threshold = MULTILOOP_IWI_THRESHOLD, list_names = LIST_NAMES
)
membership <- compute_multiloop_type_membership(whistle_chains)

pair_counts <- compute_type_cooccurrence_pairs(membership)
triplet_counts <- compute_type_cooccurrence_triplets(membership)

pairs_recurring <- pair_counts[pair_counts$n_chains >= PAIR_MIN_CHAINS, ]
triplets_recurring <- triplet_counts[triplet_counts$n_chains >= TRIPLET_MIN_CHAINS, ]

# A pair/triplet can only be drawn on this layout if every member
# sub-category is itself a node in the significant network (some
# sub-categories have no significant Markov transition at all, in either
# direction, and were dropped from gra1 entirely - see compute_graph()).
pairs_in_network <- as.character(pairs_recurring$type_a) %in% valid_node_names &
  as.character(pairs_recurring$type_b) %in% valid_node_names
triplets_in_network <- as.character(triplets_recurring$type_a) %in% valid_node_names &
  as.character(triplets_recurring$type_b) %in% valid_node_names &
  as.character(triplets_recurring$type_c) %in% valid_node_names

n_pairs_dropped <- sum(!pairs_in_network)
n_triplets_dropped <- sum(!triplets_in_network)
pairs_plot <- pairs_recurring[pairs_in_network, ]
triplets_plot <- triplets_recurring[triplets_in_network, ]

message(
  nrow(pairs_plot), " / ", nrow(pairs_recurring),
  " recurring pairs (>= ", PAIR_MIN_CHAINS, " chains) plotted (",
  n_pairs_dropped, " dropped: a member sub-category has no significant Markov transition)"
)
message(
  nrow(triplets_plot), " / ", nrow(triplets_recurring),
  " recurring triplets (>= ", TRIPLET_MIN_CHAINS, " chains) plotted (",
  n_triplets_dropped, " dropped: a member sub-category has no significant Markov transition)"
)

# Line width scales with how many multi-loop chains share that pair/triplet,
# so the most recurrent co-occurrences stand out. range() degenerates
# (min == max) whenever only one row clears the threshold; pmax(..., 1)
# below then falls back to a flat mid-range width instead of erroring.
pair_lwd <- if (nrow(pairs_plot) > 1 &&
                 diff(range(pairs_plot$n_chains)) > 0) {
  rescale_from_reference(pairs_plot$n_chains, range(pairs_plot$n_chains), c(2, 6))
} else {
  rep(4, nrow(pairs_plot))
}
triplet_lwd <- if (nrow(triplets_plot) > 1 &&
                    diff(range(triplets_plot$n_chains)) > 0) {
  rescale_from_reference(triplets_plot$n_chains, range(triplets_plot$n_chains), c(1.5, 4))
} else {
  rep(2.5, nrow(triplets_plot))
}

message("Plotting significant network with multi-loop pair overlay...")
with_pdf_plot(
  plot_path("multiloop_recurring_pairs_overlay", REPO_ROOT, subdir = "multiloops"),
  width = 9,
  height = 9,
  {
    plot_markov_graph(
      gra1,
      list_names = LIST_NAMES,
      seed = PLOT_SEED,
      coordinates = coordinates,
      show_legend = FALSE,
      caption = c(
        "background: script 01's significant Markov transition network (same nodes/layout/edges, unchanged)",
        paste0(
          "orange line = multi-loop sub-category pair co-occurring (order/position ignored) in >= ",
          PAIR_MIN_CHAINS, " multi-loop chains; line width ~ number of chains sharing that pair"
        ),
        paste0(
          n_pairs_dropped,
          " pair(s) skipped (a member sub-category has no significant Markov transition)"
        )
      )
    )

    for (i in seq_len(nrow(pairs_plot))) {
      nodes <- as.character(unlist(pairs_plot[i, c("type_a", "type_b")]))
      pts <- coordinates[nodes, , drop = FALSE]
      segments(
        pts[1, "x"], pts[1, "y"], pts[2, "x"], pts[2, "y"],
        col = adjustcolor(OVERLAY_COLOR, alpha.f = 0.75), lwd = pair_lwd[i]
      )
    }

    legend(
      "topleft", inset = 0.01, xpd = NA, bty = "n", cex = 0.7,
      legend = paste0("multi-loop pair (>= ", PAIR_MIN_CHAINS, " chains)"),
      col = OVERLAY_COLOR, lty = 1, lwd = 3
    )
  }
)

message("Plotting significant network with multi-loop triplet overlay...")
with_pdf_plot(
  plot_path("multiloop_recurring_triplets_overlay", REPO_ROOT, subdir = "multiloops"),
  width = 9,
  height = 9,
  {
    plot_markov_graph(
      gra1,
      list_names = LIST_NAMES,
      seed = PLOT_SEED,
      coordinates = coordinates,
      show_legend = FALSE,
      caption = c(
        "background: script 01's significant Markov transition network (same nodes/layout/edges, unchanged)",
        paste0(
          "orange triangle = multi-loop sub-category triplet co-occurring (order/position ignored) in >= ",
          TRIPLET_MIN_CHAINS, " multi-loop chains; border width ~ number of chains sharing that triplet"
        ),
        paste0(
          n_triplets_dropped,
          " triplet(s) skipped (a member sub-category has no significant Markov transition)"
        )
      )
    )

    for (i in seq_len(nrow(triplets_plot))) {
      nodes <- as.character(unlist(triplets_plot[i, c("type_a", "type_b", "type_c")]))
      pts <- coordinates[nodes, , drop = FALSE]
      polygon(
        pts[, "x"], pts[, "y"],
        border = OVERLAY_COLOR, lty = 1, lwd = triplet_lwd[i],
        col = adjustcolor(OVERLAY_COLOR, alpha.f = 0.25)
      )
    }

    legend(
      "topleft", inset = 0.01, xpd = NA, bty = "n", cex = 0.7,
      legend = paste0("multi-loop triplet (>= ", TRIPLET_MIN_CHAINS, " chains)"),
      col = OVERLAY_COLOR, lty = 1, lwd = 3
    )
  }
)

# Remove the earlier combined pairs+triplets overlay this pair/triplet split
# replaces.
for (ext in c("pdf", "png")) {
  stale_path <- plot_path("multiloop_recurring_overlay", REPO_ROOT, ext = ext, subdir = "multiloops")
  if (file.exists(stale_path)) file.remove(stale_path)
}

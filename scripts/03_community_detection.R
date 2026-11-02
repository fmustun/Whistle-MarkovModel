#!/usr/bin/env Rscript
#
# Manuscript Fig. 5a and Methods 4.9: Leiden community detection on the
# significant network (modularity) and per-community / per-category
# compactness (mean within-group over mean between-group transition
# probability).

suppressPackageStartupMessages({
  library(igraph)
  library(scales)
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
source(file.path(REPO_ROOT, "R", "plot_io.R"))

ensure_plots_dir(REPO_ROOT)

rds_path <- file.path(REPO_ROOT, OUTPUT_DIR, "markov_model.rds")
if (!file.exists(rds_path)) {
  stop("Run scripts/01_run_markov_model.R first. Missing: ", rds_path, call. = FALSE)
}

saved <- readRDS(rds_path)
gra1 <- saved$gra1
MarkovModel <- saved$MarkovModel

register_fcircle_shape()

gra1_undirected <- as.undirected(gra1)
set.seed(PLOT_SEED)
cl <- cluster_leiden(
  gra1_undirected,
  objective_function = "modularity",
  resolution_parameter = COMMUNITY_RESOLUTION,
  n_iterations = -1
)
# cluster_leiden(), unlike cluster_louvain(), does not store modularity on
# the returned communities object - modularity(cl) alone would error with
# "Modularity was not calculated" - so it's recomputed from the membership
# vector at the same resolution the partition was found at. weights must be
# passed explicitly: unlike cluster_leiden()/cluster_louvain(), modularity()
# does NOT fall back to the graph's own "weight" edge attribute when weights
# is omitted - it silently computes *unweighted* modularity instead, which
# would silently mismatch what cluster_leiden() actually optimized for.
cl_modularity <- modularity(
  gra1_undirected, membership(cl),
  weights = E(gra1_undirected)$weight, resolution = COMMUNITY_RESOLUTION
)

with_pdf_plot(
  plot_path("community_detection", REPO_ROOT),
  width = 9,
  height = 9,
  {
    set.seed(PLOT_SEED)
    plot(
      cl, gra1,
      col = V(gra1)$color,
      edge.curved = 0.2,
      edge.arrow.size = 1.2,
      vertex.size = 6,
      arrow.width = 0.5,
      edge.arrow.width = 3,
      vertex.shape = "fcircle",
      edge.color = rgb(10 / 255, 10 / 255, 10 / 255, 0.05),
      vertex.frame.color = V(gra1)$vertex.frame.col,
      vertex.frame.width = 2.5,
      asp = 1,
      main = "Community detection"
    )
    a <- legend(-2, 1, legend = unique(V(gra1)$category))
    x <- (a$text$x + a$rect$left) / 2
    y <- a$text$y
    symbols(
      x, y,
      circles = rep(1 / 30, length(unique(V(gra1)$category))),
      inches = FALSE, add = TRUE,
      bg = unique(V(gra1)$color), col = "gray"
    )
  }
)
message(sprintf("Modularity: %.4f", cl_modularity))

# Per-community/per-category compactness: mean within-group transition
# probability over mean between-group transition probability (see
# community_compactness(), R/graph_utils.R), computed with self-loops
# included (a self-transition is always "within" its own node's group, so
# leaving it out - as gra1 does, self-loops are stripped from it - would
# understate every group's within side). That means this needs the
# pre-loop-removal graph (bh_analysis/bh_corrected_network.rds's
# graph_with_loops), not gra1 - same vertex set/order as gra1 (self-loop
# removal only ever drops edges, never vertices), so membership(cl)/
# V(gra1)$category still line up with it vertex-for-vertex.
#
# Reported both unweighted (plain mean over edges - every significant
# transition counts once, regardless of how many times it was actually
# observed) and count-weighted (weighted by each edge's observation count,
# so a probability backed by many observations moves the mean more than one
# backed by few) side by side. The values quoted in the manuscript (SW_Luna
# 6.65, SW_Neo 4.00, SW_Nana 2.04, SW_Nikita 1.91, SW_Yosefa 2.38) are the
# count-weighted, per-main-category ones (category_compactness.csv).
bh_network_rds <- file.path(REPO_ROOT, OUTPUT_DIR, "bh_analysis", "bh_corrected_network.rds")
if (!file.exists(bh_network_rds)) {
  stop("Run scripts/01_run_markov_model.R first. Missing: ", bh_network_rds, call. = FALSE)
}
gra1_with_loops <- readRDS(bh_network_rds)$graph_with_loops
stopifnot(identical(V(gra1_with_loops)$name, V(gra1)$name))

compactness_endpoints <- ends(gra1_with_loops, E(gra1_with_loops), names = TRUE)
compactness_edge_counts <- MarkovModel$transition_counts_matrix[
  cbind(as.integer(compactness_endpoints[, 1L]), as.integer(compactness_endpoints[, 2L]))
]

# Combines an unweighted and a count-weighted community_compactness() call
# (same graph/grouping) into one side-by-side comparison table.
compactness_comparison <- function(gra, membership, edge_counts) {
  unweighted <- community_compactness(gra, membership)
  weighted <- community_compactness(gra, membership, edge_counts = edge_counts)
  merge(
    unweighted, weighted,
    by = c("group", "n_nodes", "n_within_edges", "n_between_edges"),
    suffixes = c("_unweighted", "_weighted")
  )
}

community_compactness_table <- compactness_comparison(
  gra1_with_loops, membership(cl), compactness_edge_counts
)
compactness_csv <- file.path(REPO_ROOT, OUTPUT_DIR, "community_compactness.csv")
write.csv(community_compactness_table, compactness_csv, row.names = FALSE)
message("Saved ", compactness_csv)
print(community_compactness_table)

# Same metric, grouped by main whistle category (V(gra1)$category - SW_Neo,
# SW_Luna, ...) instead of the Leiden-detected community above. A different
# grouping over the same edges, so a category can come out considerably
# more (or less) compact than the community(ies) its nodes ended up in -
# categories aren't required to align with what Leiden found.
category_compactness_table <- compactness_comparison(
  gra1_with_loops, V(gra1_with_loops)$category, compactness_edge_counts
)
category_compactness_csv <- file.path(REPO_ROOT, OUTPUT_DIR, "category_compactness.csv")
write.csv(category_compactness_table, category_compactness_csv, row.names = FALSE)
message("Saved ", category_compactness_csv)
print(category_compactness_table)

with_pdf_plot(
  plot_path("community_subgraphs", REPO_ROOT),
  width = 14,
  height = 9,
  {
    par(mfrow = c(2, 3))
    for (com in seq_along(cl)) {
      set.seed(PLOT_SEED)
      subgra1 <- subgraph(gra1, cl$membership == com)
      plot(
        subgra1,
        layout = layout_with_fr(subgra1, niter = 500),
        edge.curved = 0.2,
        vertex.size = 8,
        edge.color = rgb(180 / 255, 180 / 255, 180 / 255, rescale(E(subgra1)$weight, c(0.2, 1))),
        edge.width = rescale(log(E(subgra1)$weight), c(1, 5)),
        edge.arrow.size = 0.75,
        edge.arrow.width = 1.5,
        vertex.shape = "fcircle",
        vertex.frame.color = V(subgra1)$vertex.frame.col,
        vertex.frame.width = 2.5,
        asp = 1
      )
    }
    par(mfrow = c(1, 1))
  }
)

out_path <- file.path(REPO_ROOT, OUTPUT_DIR, "community_detection.rds")
saveRDS(
  list(
    resolution = COMMUNITY_RESOLUTION,
    membership = setNames(membership(cl), V(gra1)$name),
    community_modularity = cl_modularity,
    community_compactness = community_compactness_table,
    category_compactness = category_compactness_table
  ),
  out_path
)
message("Saved ", out_path)

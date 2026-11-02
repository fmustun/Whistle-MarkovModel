#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo_root <- if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."),
                mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(ccber))
source(file.path(repo_root, "R", "config.R"))
source(file.path(repo_root, "R", "graph_utils.R"))

# --- leiden_resolution_sweep() -----------------------------------------

# Two triangles joined by a single weak edge: an unambiguous two-community
# structure at any reasonable resolution.
g <- graph_from_edgelist(
  rbind(
    c(1, 2), c(2, 3), c(1, 3),
    c(4, 5), c(5, 6), c(4, 6),
    c(3, 4)
  ),
  directed = FALSE
)
E(g)$weight <- c(rep(5, 6), 0.1)

resolutions <- c(0.5, 1)
sweep <- leiden_resolution_sweep(g, resolutions, seed = 1)

stopifnot(
  identical(dim(sweep), c(2L, 3L)),
  identical(names(sweep), c("resolution", "n_communities", "modularity")),
  identical(sweep$resolution, resolutions),
  all(sweep$n_communities >= 1),
  all(is.finite(sweep$modularity))
)

# Reproducible given the same seed, regardless of how many resolutions are
# swept together vs. one at a time.
sweep_again <- leiden_resolution_sweep(g, resolutions, seed = 1)
sweep_one_at_a_time <- do.call(rbind, lapply(resolutions, function(r) {
  leiden_resolution_sweep(g, r, seed = 1)
}))
stopifnot(
  identical(sweep, sweep_again),
  identical(sweep$n_communities, sweep_one_at_a_time$n_communities),
  identical(sweep$modularity, sweep_one_at_a_time$modularity)
)

# The two triangles should end up in separate communities at resolution 1,
# given how weak the bridging edge is relative to the triangle edges.
set.seed(1)
cl_check <- cluster_leiden(
  g, objective_function = "modularity", weights = E(g)$weight,
  resolution_parameter = 1, n_iterations = -1
)
stopifnot(
  membership(cl_check)[1] == membership(cl_check)[2],
  membership(cl_check)[1] == membership(cl_check)[3],
  membership(cl_check)[4] == membership(cl_check)[5],
  membership(cl_check)[4] == membership(cl_check)[6],
  membership(cl_check)[1] != membership(cl_check)[4]
)
stopifnot(isTRUE(all.equal(
  sweep$modularity[sweep$resolution == 1],
  modularity(g, membership(cl_check), weights = E(g)$weight, resolution = 1)
)))

# An explicit weights argument overrides the graph's own "weight" attribute.
uniform_sweep <- leiden_resolution_sweep(
  g, 1, weights = rep(1, ecount(g)), seed = 1
)
stopifnot(!isTRUE(all.equal(uniform_sweep$modularity, sweep$modularity[2])))

# Guards the gotcha that caused a real bug in scripts/03_community_
# detection.R: unlike cluster_leiden()/cluster_louvain(), plain modularity()
# does NOT fall back to the graph's own "weight" edge attribute when
# `weights` is omitted - it silently computes *unweighted* modularity
# instead. leiden_resolution_sweep() must always pass weights explicitly, so
# its result has to match the explicit-weights call and differ from the
# no-weights one for this same weighted graph/membership.
stopifnot(
  isTRUE(all.equal(
    sweep$modularity[sweep$resolution == 1],
    modularity(g, membership(cl_check), weights = E(g)$weight, resolution = 1)
  )),
  !isTRUE(all.equal(
    sweep$modularity[sweep$resolution == 1],
    modularity(g, membership(cl_check), resolution = 1)
  ))
)

# --- reweight_graph_from_null_counts() ----------------------------------

# Template edges: (1,2) and (3,1). Node 3 never emits anything in
# null_counts (row 3 sums to zero) - that edge's null probability must come
# out as 0 (CalcTransitionMatrix's NaN -> 0 rule), not NA.
template <- graph_from_edgelist(rbind(c(1, 2), c(3, 1)), directed = TRUE)
V(template)$name <- as.character(1:3)
E(template)$weight <- c(0.9, 0.9)

null_counts <- matrix(0, nrow = 3, ncol = 3)
null_counts[1, 2] <- 3
null_counts[1, 3] <- 1
null_counts[2, 3] <- 5
null_counts[2, 1] <- 5
# row 3 left all-zero: node 3 never occurs under this null realization.

reweighted <- reweight_graph_from_null_counts(template, null_counts)

stopifnot(
  identical(ecount(reweighted), ecount(template)),
  identical(
    ends(reweighted, E(reweighted), names = TRUE),
    ends(template, E(template), names = TRUE)
  ),
  # (1,2): 3 of node 1's 4 outgoing counts go to node 2 -> 3/4.
  isTRUE(all.equal(E(reweighted)$weight[1], 0.75)),
  # (3,1): node 3 has zero outgoing counts -> NaN -> 0, not NA.
  isTRUE(all.equal(E(reweighted)$weight[2], 0))
)

# --- community_compactness() --------------------------------------------

# Communities: A = {1,2,3}, B = {4,5} (both connected, to each other and
# internally), C = {6} (an isolated singleton, no edges at all), D = {7,8}
# (an isolated *component* - has an internal edge, but nothing connecting
# it to any other community).
compactness_edges <- rbind(
  c(1, 2), c(2, 3), # within A
  c(3, 4), c(5, 1), # between A and B
  c(4, 5), # within B
  c(7, 8) # within D, but D has no edges to anywhere else
)
compactness_graph <- graph_from_edgelist(compactness_edges, directed = TRUE)
E(compactness_graph)$weight <- c(0.6, 0.4, 0.2, 0.4, 0.8, 0.9)
compactness_membership <- c(1, 1, 1, 2, 2, 3, 4, 4)

compactness <- community_compactness(compactness_graph, compactness_membership)
compactness <- compactness[order(compactness$group), ]

stopifnot(
  vcount(compactness_graph) == 8L,
  identical(compactness$group, c(1, 2, 3, 4)),
  identical(compactness$n_nodes, c(3L, 2L, 1L, 2L)),
  identical(compactness$n_within_edges, c(2L, 1L, 0L, 1L)),
  # between is strictly outgoing now: A's between-set is just 3->4 (not
  # also 5->1, which is B's outgoing / A's incoming); B's is just 5->1.
  identical(compactness$n_between_edges, c(1L, 1L, 0L, 0L)),
  # A: within = mean(0.6, 0.4) = 0.5; between = its own outgoing edge
  # 3->4 only (0.2) - 5->1 is incoming (B's outgoing), excluded.
  isTRUE(all.equal(compactness$mean_within_probability[1], 0.5)),
  isTRUE(all.equal(compactness$mean_between_probability[1], 0.2)),
  isTRUE(all.equal(compactness$compactness[1], 0.5 / 0.2)),
  # B: within = 0.8; between = its own outgoing edge 5->1 only (0.4) -
  # 3->4 is incoming (A's outgoing), excluded.
  isTRUE(all.equal(compactness$mean_within_probability[2], 0.8)),
  isTRUE(all.equal(compactness$mean_between_probability[2], 0.4)),
  isTRUE(all.equal(compactness$compactness[2], 0.8 / 0.4)),
  # C: a singleton with no edges at all - both means, and compactness,
  # are NA (not 0 or dropped).
  is.na(compactness$mean_within_probability[3]),
  is.na(compactness$mean_between_probability[3]),
  is.na(compactness$compactness[3]),
  # D: has an internal edge (mean_within defined) but no outgoing edges to
  # any other community, so mean_between and compactness are NA.
  isTRUE(all.equal(compactness$mean_within_probability[4], 0.9)),
  is.na(compactness$mean_between_probability[4]),
  is.na(compactness$compactness[4])
)

# edge_counts weights each edge's contribution by its observation count
# instead of a plain mean - edges in the same order as compactness_edges
# above (1->2, 2->3, 3->4, 5->1, 4->5, 7->8).
compactness_counts <- c(10, 30, 5, 15, 1, 100)
weighted_compactness <- community_compactness(
  compactness_graph, compactness_membership, edge_counts = compactness_counts
)
weighted_compactness <- weighted_compactness[order(weighted_compactness$group), ]

stopifnot(
  # A: within = (0.6*10 + 0.4*30) / (10+30) = 0.45;
  #    between = its own outgoing edge 3->4 only (weight 0.2, count 5) = 0.2.
  isTRUE(all.equal(weighted_compactness$mean_within_probability[1], 0.45)),
  isTRUE(all.equal(weighted_compactness$mean_between_probability[1], 0.2)),
  isTRUE(all.equal(weighted_compactness$compactness[1], 0.45 / 0.2)),
  # B: within = 0.8 (single edge, its count doesn't matter on its own);
  #    between = its own outgoing edge 5->1 only (weight 0.4, count 15) = 0.4.
  isTRUE(all.equal(weighted_compactness$mean_within_probability[2], 0.8)),
  isTRUE(all.equal(weighted_compactness$mean_between_probability[2], 0.4)),
  isTRUE(all.equal(weighted_compactness$compactness[2], 0.8 / 0.4)),
  # Unweighted and weighted differ here (A's within isn't a straight 0.5/0.5
  # split of counts) - confirms edge_counts actually changes the result
  # rather than silently being ignored.
  !isTRUE(all.equal(
    weighted_compactness$mean_within_probability[1],
    compactness$mean_within_probability[1]
  ))
)

message("All graph_utils.R tests passed")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(scales)
  library(tnet)
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

rds_path <- file.path(REPO_ROOT, OUTPUT_DIR, "markov_model.rds")
if (!file.exists(rds_path)) {
  stop("Run scripts/01_run_markov_model.R first. Missing: ", rds_path, call. = FALSE)
}

saved <- readRDS(rds_path)
gra1 <- saved$gra1

register_fcircle_shape()

message("Generating ", NULL_NETWORK_ITERATIONS, " degree-preserving random networks...")
netObs <- cbind(get.edgelist(gra1, names = FALSE), E(gra1)$weight)

rdm_path_length <- rep(NA_real_, NULL_NETWORK_ITERATIONS)
rdm_clustering_coefficients <- rep(NA_real_, NULL_NETWORK_ITERATIONS)
rdm_diameter <- rep(NA_real_, NULL_NETWORK_ITERATIONS)
rdm_reciprocity <- rep(NA_real_, NULL_NETWORK_ITERATIONS)
rdm_nb_community <- rep(NA_real_, NULL_NETWORK_ITERATIONS)
rdm_modularity <- rep(NA_real_, NULL_NETWORK_ITERATIONS)

for (i in seq_len(NULL_NETWORK_ITERATIONS)) {
  netRdm <- rg_reshuffling_w(netObs, option = "links", directed = TRUE, seed = i)
  netRdm <- tnet_igraph(netRdm)
  rdm_clustering_coefficients[[i]] <- transitivity(netRdm, type = "global")
  rdm_path_length[[i]] <- mean_distance(netRdm)
  rdm_diameter[[i]] <- diameter(netRdm)
  rdm_reciprocity[[i]] <- reciprocity(netRdm)
  cl_rdm <- cluster_louvain(as.undirected(netRdm))
  rdm_nb_community[i] <- length(cl_rdm)
  rdm_modularity[i] <- modularity(cl_rdm)
}

obs_clustering_coefficients <- transitivity(gra1, type = "global")

aclust <- rdm_clustering_coefficients[rdm_clustering_coefficients > obs_clustering_coefficients]
p1clust <- length(aclust) / length(rdm_clustering_coefficients)
bclust <- rdm_clustering_coefficients[rdm_clustering_coefficients < obs_clustering_coefficients]
p0clust <- length(bclust) / length(rdm_clustering_coefficients)
p2clust <- 2 * min(p0clust, p1clust)

mean_CCrd <- mean(rdm_clustering_coefficients)
sd_CCrd <- sd(rdm_clustering_coefficients)
Z_score_CC <- (obs_clustering_coefficients - mean_CCrd) / sd_CCrd

print(
  ggplot(data = data.frame(rdm_clustering_coefficients), aes(x = rdm_clustering_coefficients)) +
    geom_histogram(aes(y = after_stat(density)), bins = 100, fill = "gray", position = "identity", linewidth = 0.8) +
    geom_vline(aes(xintercept = obs_clustering_coefficients, color = "Observed"), linewidth = 0.6) +
    labs(
      title = "Distribution of clustering coefficients on random networks",
      caption = sprintf("P-value = %.2f. Z-score = %.3f.", p2clust, Z_score_CC),
      x = "clustering coefficient",
      color = "Legend"
    )
)

SWC <- (transitivity(gra1, type = "global") / mean(rdm_clustering_coefficients)) /
  (mean_distance(gra1) / mean(rdm_path_length))
message(sprintf("Small-world coefficient: %.4f", SWC))

set.seed(PLOT_SEED)
cl <- cluster_louvain(as.undirected(gra1), resolution = 0.7)

set.seed(PLOT_SEED)
plot(
  cl, gra1,
  col = V(gra1)$color,
  edge.curved = 0.2,
  edge.arrow.size = 0.28,
  vertex.size = 6,
  arrow.width = 0.5,
  edge.arrow.width = 1,
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
message(sprintf("Modularity: %.4f", modularity(cl)))

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
    edge.arrow.size = 0.3,
    edge.arrow.width = 0.45,
    vertex.shape = "fcircle",
    vertex.frame.color = V(subgra1)$vertex.frame.col,
    vertex.frame.width = 2.5,
    asp = 1,
    main = paste("Community", com)
  )
}
par(mfrow = c(1, 1))

out_path <- file.path(REPO_ROOT, OUTPUT_DIR, "global_null_metrics.rds")
saveRDS(
  list(
    rdm_clustering_coefficients = rdm_clustering_coefficients,
    rdm_path_length = rdm_path_length,
    rdm_diameter = rdm_diameter,
    rdm_reciprocity = rdm_reciprocity,
    rdm_nb_community = rdm_nb_community,
    rdm_modularity = rdm_modularity,
    obs_clustering_coefficients = obs_clustering_coefficients,
    p2clust = p2clust,
    Z_score_CC = Z_score_CC,
    small_world_coefficient = SWC,
    community_modularity = modularity(cl)
  ),
  out_path
)
message("Saved ", out_path)

#!/usr/bin/env Rscript
#
# Leiden community-detection resolution sweep for the network built in
# scripts/03_global_network_measures.R: scans COMMUNITY_RESOLUTION_RANGE
# (R/config.R) looking for a plateau in community count / modularity, and
# benchmarks the observed network's modularity at each resolution against
# a subsample of the shuffle-null count matrices already generated for the
# edge-significance test (MarkovModel$random_whistle_occurrences_matrices -
# see MarkovWhistle's null_model = "shuffle", R/markov_whistle.R). Each null
# realization reweights gra1's own significant edges with that null's
# transition probabilities, holding topology fixed (see
# reweight_graph_from_null_counts(), R/graph_utils.R) so the comparison
# isolates the observed weight pattern rather than differences in network
# density/size.

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(ccber)
  library(gridExtra)
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

ensure_plots_dir(REPO_ROOT, "network_measures")

rds_path <- file.path(REPO_ROOT, OUTPUT_DIR, "markov_model.rds")
if (!file.exists(rds_path)) {
  stop("Run scripts/01_run_markov_model.R first. Missing: ", rds_path, call. = FALSE)
}

saved <- readRDS(rds_path)
gra1 <- saved$gra1
MarkovModel <- saved$MarkovModel

n_res <- length(COMMUNITY_RESOLUTION_RANGE)
message(
  "Sweeping Leiden resolution_parameter over ", n_res, " values (",
  min(COMMUNITY_RESOLUTION_RANGE), " to ", max(COMMUNITY_RESOLUTION_RANGE),
  ", step ", COMMUNITY_RESOLUTION_RANGE[2] - COMMUNITY_RESOLUTION_RANGE[1], ")..."
)
gra1_undirected <- as.undirected(gra1)
real_sweep <- leiden_resolution_sweep(
  gra1_undirected, COMMUNITY_RESOLUTION_RANGE, seed = PLOT_SEED
)

message(
  "Benchmarking against ", COMMUNITY_NULL_REALIZATIONS, " of ",
  MarkovModel$iterations_number,
  " saved shuffle-null realizations (edge topology held fixed at gra1's",
  " significant edges; weights drawn from each null's transition",
  " probabilities)..."
)
set.seed(PLOT_SEED)
null_indices <- sample(
  seq_len(MarkovModel$iterations_number), COMMUNITY_NULL_REALIZATIONS
)

null_modularity_mat <- matrix(
  NA_real_, nrow = COMMUNITY_NULL_REALIZATIONS, ncol = n_res
)
null_communities_mat <- matrix(
  NA_integer_, nrow = COMMUNITY_NULL_REALIZATIONS, ncol = n_res
)
for (k in seq_along(null_indices)) {
  null_gra <- reweight_graph_from_null_counts(
    gra1, MarkovModel$random_whistle_occurrences_matrices[[null_indices[k]]]
  )
  null_gra_undirected <- as.undirected(null_gra)
  sweep_k <- leiden_resolution_sweep(
    null_gra_undirected, COMMUNITY_RESOLUTION_RANGE, seed = PLOT_SEED
  )
  null_modularity_mat[k, ] <- sweep_k$modularity
  null_communities_mat[k, ] <- sweep_k$n_communities
  if (k %% 20 == 0 || k == length(null_indices)) {
    message("  null realization ", k, " / ", length(null_indices))
  }
}

# Two-tailed empirical p, same construction as p2clust earlier in
# scripts/03_global_network_measures.R (2 * min(one-tailed p in each
# direction)) - real modularity can plausibly land in *either* tail: below
# the null (a coarser/less-modular partition than random reweighting would
# give, seen at low resolution) as well as above it (seen at high
# resolution), and a one-tailed "is observed higher" test would report the
# below-null case as merely non-significant (p near 1) rather than flagging
# it as significant in the other direction. Each one-tailed p keeps the same
# Laplace-smoothed, tie-inclusive formula used elsewhere in this codebase -
# see empirical_enrichment_p_values(), R/statistical_utils.R.
observed_modularity_mat <- matrix(
  real_sweep$modularity, nrow = COMMUNITY_NULL_REALIZATIONS, ncol = n_res,
  byrow = TRUE
)
exceedance_count_high <- colSums(null_modularity_mat >= observed_modularity_mat)
exceedance_count_low <- colSums(null_modularity_mat <= observed_modularity_mat)
p_modularity_high <- (1 + exceedance_count_high) / (COMMUNITY_NULL_REALIZATIONS + 1)
p_modularity_low <- (1 + exceedance_count_low) / (COMMUNITY_NULL_REALIZATIONS + 1)
p_modularity_two_tailed <- pmin(1, 2 * pmin(p_modularity_high, p_modularity_low))

sweep_summary <- data.frame(
  resolution = COMMUNITY_RESOLUTION_RANGE,
  n_communities_real = real_sweep$n_communities,
  n_communities_null_median = apply(null_communities_mat, 2, median),
  n_communities_null_lo = apply(null_communities_mat, 2, quantile, probs = 0.025, names = FALSE),
  n_communities_null_hi = apply(null_communities_mat, 2, quantile, probs = 0.975, names = FALSE),
  modularity_real = real_sweep$modularity,
  modularity_null_median = apply(null_modularity_mat, 2, median),
  modularity_null_lo = apply(null_modularity_mat, 2, quantile, probs = 0.025, names = FALSE),
  modularity_null_hi = apply(null_modularity_mat, 2, quantile, probs = 0.975, names = FALSE),
  empirical_p_modularity_high = p_modularity_high,
  empirical_p_modularity_low = p_modularity_low,
  empirical_p_modularity_two_tailed = p_modularity_two_tailed
)
sweep_summary$modularity_significant <- sweep_summary$empirical_p_modularity_two_tailed < 0.05

message(sprintf(
  "Observed modularity peak: %.4f at resolution %.2f",
  max(real_sweep$modularity), real_sweep$resolution[which.max(real_sweep$modularity)]
))
message(sprintf(
  "Community count range across sweep: %d - %d",
  min(real_sweep$n_communities), max(real_sweep$n_communities)
))
selected_row <- sweep_summary[
  which.min(abs(sweep_summary$resolution - COMMUNITY_RESOLUTION)),
]
message(sprintf(
  "At the selected resolution (%.2f): %d communities, modularity = %.4f (null median %.4f), two-tailed p = %.4f (%s)",
  selected_row$resolution, selected_row$n_communities_real,
  selected_row$modularity_real, selected_row$modularity_null_median,
  selected_row$empirical_p_modularity_two_tailed,
  if (selected_row$modularity_real < selected_row$modularity_null_median) {
    "below null"
  } else {
    "above null"
  }
))

p_communities <- ggplot(sweep_summary, aes(x = resolution)) +
  geom_ribbon(
    aes(ymin = n_communities_null_lo, ymax = n_communities_null_hi),
    fill = "gray70", alpha = 0.5
  ) +
  geom_line(aes(y = n_communities_null_median), color = "gray40", linetype = "dashed") +
  geom_line(aes(y = n_communities_real), color = "steelblue", linewidth = 0.9) +
  geom_point(aes(y = n_communities_real), color = "steelblue", size = 1.2) +
  geom_vline(xintercept = COMMUNITY_RESOLUTION, linetype = "dotted", color = "black") +
  theme_minimal() +
  labs(
    title = "Leiden community count vs. resolution",
    subtitle = paste0(
      "blue = observed network; gray band = ", COMMUNITY_NULL_REALIZATIONS,
      " shuffle-null realizations (95% band, dashed = median);\n",
      "dotted line = resolution used in the community_detection figure"
    ),
    x = "resolution_parameter",
    y = "number of communities"
  )

p_modularity <- ggplot(sweep_summary, aes(x = resolution)) +
  geom_ribbon(
    aes(ymin = modularity_null_lo, ymax = modularity_null_hi),
    fill = "gray70", alpha = 0.5
  ) +
  geom_line(aes(y = modularity_null_median), color = "gray40", linetype = "dashed") +
  geom_line(aes(y = modularity_real), color = "firebrick", linewidth = 0.9) +
  geom_point(aes(y = modularity_real), color = "firebrick", size = 1.2) +
  geom_point(
    data = subset(sweep_summary, modularity_significant),
    aes(y = modularity_real), shape = 8, size = 2.6, color = "black"
  ) +
  geom_vline(xintercept = COMMUNITY_RESOLUTION, linetype = "dotted", color = "black") +
  theme_minimal() +
  labs(
    title = "Leiden modularity vs. resolution",
    subtitle = paste0(
      "red = observed network; gray band = shuffle-null 95% envelope (dashed = median);\n",
      "black asterisk = two-tailed empirical p < 0.05 (either above or below the null)"
    ),
    x = "resolution_parameter",
    y = "modularity"
  )

save_ggplot(
  gridExtra::arrangeGrob(p_communities, p_modularity, ncol = 1),
  "leiden_resolution_sweep",
  REPO_ROOT,
  width = 8,
  height = 9,
  subdir = "network_measures"
)

csv_path <- file.path(REPO_ROOT, OUTPUT_DIR, "leiden_resolution_sweep.csv")
write.csv(sweep_summary, csv_path, row.names = FALSE)
message("Saved ", csv_path)

out_path <- file.path(REPO_ROOT, OUTPUT_DIR, "leiden_resolution_sweep.rds")
saveRDS(
  list(
    resolutions = COMMUNITY_RESOLUTION_RANGE,
    real_sweep = real_sweep,
    null_indices = null_indices,
    null_modularity = null_modularity_mat,
    null_n_communities = null_communities_mat,
    summary = sweep_summary
  ),
  out_path
)
message("Saved ", out_path)

#!/usr/bin/env Rscript
#
# Manuscript Fig. S4 and Methods 4.9: node strength (self-transitions
# included) against emission frequency, tested against the shuffle null
# that preserves each node's occurrence count; BH-corrected per-node
# empirical p-values define the hubs.

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
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
source(file.path(REPO_ROOT, "R", "statistical_utils.R"))
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

# Node strength here is the raw total transition count in + out of a node
# (self-transitions included, as in Methods 4.9), not the row-normalized
# probability strength of the plotted network - see the comment on
# strength_null_distribution() in R/statistical_utils.R for why raw counts
# are the quantity that can be nulled meaningfully.
message(
  "Computing frequency-based null for node strength (",
  MarkovModel$iterations_number, " shuffle iterations)..."
)
strength_null_incl_self <- strength_null_distribution(
  MarkovModel$transition_counts_matrix,
  MarkovModel$random_whistle_occurrences_matrices,
  exclude_self = FALSE
)

n_nodes <- nrow(MarkovModel$transition_counts_matrix)
node_ids <- as.character(seq_len(n_nodes))

# assign_vertex_category_color() only needs vertex names, so a bare
# n-vertex graph gets every whistle sub-category (86) a category/color
# without going through compute_graph() - which drops any node with zero
# transitions (there is exactly one: whistle type 76, occurrences = 4)
# and would silently lose it from this frequency-driven analysis.
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

# Empirical p-value per node (rank-based: (1 + count(null >= observed)) /
# (iterations + 1), same formula as empirical_enrichment_p_values() used
# for edges) rather than reading significance off the z-score directly -
# the z-score assumes an approximately normal null, which is a much rougher
# approximation for low-count nodes (see strength_null_distribution()'s
# p_value in R/statistical_utils.R). BH-corrected across the family of all
# 86 node-level tests at GRAPH_ALPHA (R/config.R), the same alpha and
# correction already used for the edge-level significant network in
# scripts/01_run_markov_model.R - so "significant" here means the same
# thing it means everywhere else in this pipeline.
bh_q_value_incl_self <- stats::p.adjust(strength_null_incl_self$p_value, method = "BH")

df <- data.frame(
  nodes = factor(node_ids, levels = node_ids),
  sub_category = unname(whistle_chr_map[node_ids]),
  category = V(node_graph)$category,
  color = V(node_graph)$color,
  occurrences = as.numeric(MarkovModel$whistle_occurrence[node_ids]),
  strength_incl_self = strength_null_incl_self$observed_strength,
  null_mean_strength_incl_self = strength_null_incl_self$null_mean,
  null_sd_strength_incl_self = strength_null_incl_self$null_sd,
  z_score_incl_self = strength_null_incl_self$z_score,
  p_value_incl_self = strength_null_incl_self$p_value,
  bh_q_value_incl_self = bh_q_value_incl_self
)

csv_path <- file.path(REPO_ROOT, OUTPUT_DIR, "node_strength_frequency_null.csv")
write.csv(df, csv_path, row.names = FALSE, quote = FALSE)
message("Saved ", csv_path)

# The two panels of Fig. S4 are built by these two helpers, which return a
# ggplot object rather than saving directly so they can be combined into a
# single two-panel figure below.
#
# The gray band is the same BH-corrected significance boundary as the
# z-score panel's dashed line (z_threshold, computed below by
# z_score_bh_threshold()), converted back to strength units per node as
# null_mean + z_threshold * null_sd - not an arbitrary +/-2 SD. It only has
# an upper edge (ymin fixed at 0) because strength_null_distribution()'s
# p-value is one-sided (tests enrichment above the null, never depletion
# below it - see R/statistical_utils.R), so there is no tested boundary for
# a symmetric lower band. A node is "in the gray zone" here exactly when its
# z-score falls below the dashed line in the paired z-score plot.
build_strength_vs_frequency_null_plot <- function(
    df, strength_col, null_mean_col, null_sd_col, z_threshold, y_label,
    x_label = "emission frequency (occurrences)", caption = NULL, legend = TRUE
) {
  df_null_curve <- df[order(df$occurrences), ]
  p <- ggplot(df) +
    aes(x = occurrences, y = .data[[strength_col]]) +
    (if (!is.na(z_threshold)) {
      geom_ribbon(
        data = df_null_curve,
        aes(
          x = occurrences,
          y = .data[[null_mean_col]],
          ymin = 0,
          ymax = .data[[null_mean_col]] + z_threshold * .data[[null_sd_col]]
        ),
        fill = "gray70", alpha = 0.4, inherit.aes = FALSE
      )
    }) +
    geom_line(
      data = df_null_curve,
      aes(x = occurrences, y = .data[[null_mean_col]]),
      color = "gray30", linetype = "dashed", inherit.aes = FALSE
    ) +
    geom_point(aes(fill = category), colour = "black", pch = 21, size = 5) +
    scale_fill_manual(breaks = unique(df$category), values = unique(df$color)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 13)
    ) +
    labs(x = x_label, y = y_label, caption = caption) +
    geom_text_repel(aes(label = nodes), size = 4.5)
  if (!legend) p <- p + theme(legend.position = "none")
  p
}

# z_threshold replaces an arbitrary "z = 2" convention with the boundary
# actually implied by BH-corrected significance: the smallest z-score among
# nodes whose empirical p-value survives BH correction at GRAPH_ALPHA (NA,
# drawing no line, if nothing survives). Shared with the paired
# strength-vs-frequency plot's gray band above, so the two panels agree on
# exactly which nodes count as significant.
build_z_score_vs_frequency_plot <- function(
    df, z_col, z_threshold, bh_q_col, x_label = "emission frequency (occurrences)",
    caption = NULL, legend = TRUE, annotate_threshold = FALSE
) {
  p <- ggplot(df) +
    aes(x = occurrences, y = .data[[z_col]]) +
    geom_hline(yintercept = 0, color = "gray50", linewidth = 0.4) +
    (if (!is.na(z_threshold)) {
      geom_hline(yintercept = z_threshold, color = "black", linetype = "dashed", linewidth = 0.6)
    }) +
    (if (!is.na(z_threshold) && annotate_threshold) {
      annotate(
        "text", x = Inf, y = z_threshold,
        label = sprintf("z = %.2f", z_threshold),
        hjust = 1.05, vjust = 1.5, size = 6
      )
    }) +
    geom_point(aes(fill = category), colour = "black", pch = 21, size = 5) +
    scale_fill_manual(breaks = unique(df$category), values = unique(df$color)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 13)
    ) +
    labs(
      x = x_label,
      y = "standardized strength (z-score vs. shuffle null)",
      caption = caption
    ) +
    # Only the BH-significant nodes get a label - labeling all 86 was
    # illegible (ggrepel was already dropping most of them as "too many
    # overlaps"). Uses bh_q_col directly (the exact, distribution-free test)
    # rather than z_col >= z_threshold: the dashed line is a single z-value
    # shared across all nodes, which only marks the same nodes as
    # BH-significant if every node's null distribution is close enough to
    # Gaussian for a common z to imply a common p-value - true here (checked
    # against the per-node empirical quantiles directly) but not guaranteed
    # on a future rerun with different data, so labeling shouldn't depend on
    # it.
    geom_text_repel(aes(
      label = ifelse(
        !is.na(.data[[bh_q_col]]) & .data[[bh_q_col]] < GRAPH_ALPHA,
        as.character(nodes), NA_character_
      )
    ), size = 4.5)
  if (!legend) p <- p + theme(legend.position = "none")
  p
}

# Smallest z-score among BH-significant nodes - the z-score value at which
# the dashed line above actually separates BH-significant from
# not-significant, IF z-score and empirical-p rank are concordant. That's
# only guaranteed when every node's null distribution is close enough to
# Gaussian that one shared z maps to a shared p-value; it isn't guaranteed in
# general (low-count nodes' null strength is bounded at 0, so their null can
# skew), which is why the label logic in build_z_score_vs_frequency_plot()
# uses bh_q directly rather than this threshold. This value only drives the
# dashed line / shaded band's *position*; warn_if_z_threshold_discordant()
# below checks that position is still trustworthy for the actual data.
z_score_bh_threshold <- function(z_score, bh_q, alpha = GRAPH_ALPHA) {
  significant <- !is.na(bh_q) & bh_q < alpha
  if (!any(significant)) return(NA_real_)
  min(z_score[significant], na.rm = TRUE)
}

# Confirms the single z_threshold line doesn't visually mislead: warns if any
# node would look BH-significant from the z-score plot (z >= z_threshold)
# without actually being BH-significant (bh_q >= alpha), which is exactly
# the failure mode a skewed per-node null could cause.
warn_if_z_threshold_discordant <- function(z_score, bh_q, z_threshold, label,
                                            alpha = GRAPH_ALPHA) {
  if (is.na(z_threshold)) return(invisible(NULL))
  looks_significant <- !is.na(z_score) & z_score >= z_threshold
  is_significant <- !is.na(bh_q) & bh_q < alpha
  discordant <- looks_significant & !is_significant
  if (any(discordant)) {
    warning(
      sprintf(
        "%s: z_threshold = %.3f visually flags node(s) %s as significant but BH q-value does not (null distribution likely skewed for these nodes) - dashed line/band position is misleading for them.",
        label, z_threshold, paste(which(discordant), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

z_threshold_incl_self <- z_score_bh_threshold(df$z_score_incl_self, df$bh_q_value_incl_self)
warn_if_z_threshold_discordant(
  df$z_score_incl_self, df$bh_q_value_incl_self, z_threshold_incl_self, "incl_self"
)

# One two-panel figure (Fig. S4): no per-panel caption or legend (the dashed
# BH-threshold line in the z-score panel is labeled directly instead).
strength_panel_incl_self <- build_strength_vs_frequency_null_plot(
  df, "strength_incl_self", "null_mean_strength_incl_self", "null_sd_strength_incl_self",
  z_threshold_incl_self,
  y_label = "strength (self-transitions included)",
  x_label = "number of occurrences",
  legend = FALSE
)
z_score_panel_incl_self <- build_z_score_vs_frequency_plot(
  df, "z_score_incl_self", z_threshold_incl_self, "bh_q_value_incl_self",
  x_label = "number of occurrences",
  legend = FALSE,
  annotate_threshold = TRUE
)
save_ggplot(
  gridExtra::arrangeGrob(strength_panel_incl_self, z_score_panel_incl_self, ncol = 2),
  "strength_vs_frequency_null_incl_self",
  REPO_ROOT,
  width = 14,
  height = 6,
  subdir = "network_measures"
)

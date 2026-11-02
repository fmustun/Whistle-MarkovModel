#!/usr/bin/env Rscript
#
# Manuscript Fig. 5b-c and Methods 4.9: node-level measures of the
# significant network - turn-taking (probability that the next whistle
# belongs to a different main category, Fig. 5b) and standardized strength
# (z-score vs. the shuffle null, self-transitions included) against
# betweenness (Fig. 5c). Writes outputs/node_measures.csv, read by
# scripts/11_multiloop_position_role.R.

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(scales)
  library(dplyr)
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

ensure_plots_dir(REPO_ROOT)
ensure_plots_dir(REPO_ROOT, "network_measures")

rds_path <- file.path(REPO_ROOT, OUTPUT_DIR, "markov_model.rds")
if (!file.exists(rds_path)) {
  stop("Run scripts/01_run_markov_model.R first. Missing: ", rds_path, call. = FALSE)
}

saved <- readRDS(rds_path)
MarkovModel <- saved$MarkovModel
whistles_list <- saved$whistles_list
gra1 <- saved$gra1
list_names <- saved$config$list_names

deg <- degree(gra1, mode = "total")
deg_in <- degree(gra1, mode = "in")
deg_out <- degree(gra1, mode = "out")
nodes_betweenness <- betweenness(gra1)
nodes_strength <- strength(gra1, mode = "all")
nodes_strength_in <- strength(gra1, mode = "in")
nodes_strength_out <- strength(gra1, mode = "out")
nodes_closeness <- closeness(gra1)

seq_category <- build_seq_category(whistles_list, list_names)

col_same_cat_outbound <- c()
col_diff_cat_outbound <- c()
col_same_cat_inbound <- c()
col_diff_cat_inbound <- c()

for (node in V(gra1)$name) {
  node_i <- as.integer(node)
  i <- 1
  while (!node_i %in% seq_category[[list_names[i]]]) {
    i <- i + 1
  }
  sub_network <- list_names[i]
  subv <- seq_category[[sub_network]]

  same_cat_occ_outbound <- sum(MarkovModel$transition_counts_matrix[node_i, subv])
  total_occ_outbound <- sum(MarkovModel$transition_counts_matrix[node_i, ])
  diff_cat_occ_outbound <- total_occ_outbound - same_cat_occ_outbound

  same_cat_occ_inbound <- sum(MarkovModel$transition_counts_matrix[subv, node_i])
  total_occ_inbound <- sum(MarkovModel$transition_counts_matrix[, node_i])
  diff_cat_occ_inbound <- total_occ_inbound - same_cat_occ_inbound

  col_same_cat_outbound <- c(col_same_cat_outbound, same_cat_occ_outbound)
  col_diff_cat_outbound <- c(col_diff_cat_outbound, diff_cat_occ_outbound)
  col_same_cat_inbound <- c(col_same_cat_inbound, same_cat_occ_inbound)
  col_diff_cat_inbound <- c(col_diff_cat_inbound, diff_cat_occ_inbound)
}

# Strength with self-transitions included, frequency-normalized: the same
# raw-count total-degree quantity used for the null model in
# scripts/05_node_strength_frequency_null.R (rowSums + colSums of
# transition_counts_matrix, self-transitions left in), divided by each
# node's own occurrence count so common and rare whistle types are
# comparable on the same axis - "transitions per whistle produced" rather
# than a raw total that mostly just tracks frequency.
node_ids_int <- as.integer(V(gra1)$name)
strength_incl_self <- rowSums(MarkovModel$transition_counts_matrix)[node_ids_int] +
  colSums(MarkovModel$transition_counts_matrix)[node_ids_int]

# Same quantity, normalized a different way: instead of dividing by
# occurrences (a simple rate), standardize against the frequency-based
# shuffle null from scripts/05_node_strength_frequency_null.R (mean/SD of
# the same raw-count strength, self-transitions included, recomputed across
# all 150k null iterations already saved on MarkovModel) - a z-score
# already accounts for the nonlinear relationship between frequency and
# expected strength, rather than assuming a flat per-occurrence rate. This
# is the strength axis of Fig. 5c.
message(
  "Computing frequency-based null for strength (self-transitions included, ",
  MarkovModel$iterations_number, " shuffle iterations)..."
)
strength_null_incl_self <- strength_null_distribution(
  MarkovModel$transition_counts_matrix,
  MarkovModel$random_whistle_occurrences_matrices,
  exclude_self = FALSE
)
z_score_incl_self <- strength_null_incl_self$z_score[node_ids_int]

df3 <- data.frame(
  nodes = V(gra1)$name,
  category = V(gra1)$category,
  sub_category = V(gra1)$sub_category,
  color = V(gra1)$color,
  degree = deg,
  degree_in = deg_in,
  degree_out = deg_out,
  betweenness = nodes_betweenness,
  strength = nodes_strength,
  strength_in = nodes_strength_in,
  strength_out = nodes_strength_out,
  strength_incl_self_per_occurrence = strength_incl_self / V(gra1)$occurrences,
  z_score_incl_self = z_score_incl_self,
  diff_degree = deg_in - deg_out,
  closeness = nodes_closeness,
  occurrences = V(gra1)$occurrences,
  prob_produced = V(gra1)$occurrences / sum(V(gra1)$occurrences),
  prob_same_cat_outbound = col_same_cat_outbound /
    (col_same_cat_outbound + col_diff_cat_outbound),
  prob_diff_cat_outbound = col_diff_cat_outbound /
    (col_same_cat_outbound + col_diff_cat_outbound),
  prob_same_cat_inbound = col_same_cat_inbound /
    (col_same_cat_inbound + col_diff_cat_inbound),
  prob_diff_cat_inbound = col_diff_cat_inbound /
    (col_same_cat_inbound + col_diff_cat_inbound)
)
df3$nodes <- as.factor(df3$nodes)

# Fig. 5b: turn-taking, P(following node in a different category).
df_order <- order_nodes_by_category(df3, list_names, "prob_diff_cat_outbound")

g1 <- ggplot(df_order) +
  aes(x = nodes, y = prob_diff_cat_outbound, fill = category) +
  geom_col(color = "black", width = 0.75) +
  scale_fill_manual(breaks = unique(df3$category), values = unique(df3$color)) +
  theme(
    axis.text.x = element_text(angle = 45),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  labs(
    title = "P(following node in a different category)",
    y = "probability"
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.75)

# One row per category group, at the mean/min/max position of that group's
# bars in df_order's category-grouped factor levels (order_nodes_by_
# category() above) - x_mid prints each group's name once, centered under
# its bars, below the per-node id labels; x_min/x_max delimit the bracket
# bar drawn under each group further down. The three NSW_* categories are
# merged into one "NSW" group here (LIST_NAMES keeps them contiguous, so
# their combined bar block is still one contiguous run) - each individually
# only spans 2-3 bars, too narrow for even a short horizontal label without
# colliding with its NSW neighbor.
category_label_positions <- df_order %>%
  dplyr::mutate(
    x_pos = as.numeric(nodes),
    label_category = ifelse(grepl("^NSW", category), "NSW", as.character(category))
  ) %>%
  dplyr::group_by(label_category) %>%
  dplyr::summarise(
    x_mid = mean(x_pos), x_min = min(x_pos), x_max = max(x_pos), .groups = "drop"
  ) %>%
  dplyr::rename(category = label_category)
category_label_positions <- category_label_positions[order(category_label_positions$x_min), ]

# Bracket delimiting each group: a horizontal bar from its first bar's left
# edge to its last bar's right edge (half of geom_col's own width = 0.75 out
# past x_min/x_max, so it lines up with the actual bar edges, not their
# centers), with a short vertical tick at each end - "|----|" under each
# group's numbers. bracket_y = -0.38 sits in the empty gap between the
# angled node-id row (whose 2-digit labels extend down to about y = -0.3)
# and the category-name row (moved to y = -0.5 below to keep clear of the
# bracket) - checked directly against a rendered proof rather than guessed,
# since axis text position isn't in data units at all (it's placed by the
# theme, in the plot margin) while annotate() layers are, so the two only
# align by comparing pixels.
bracket_y <- -0.38
bracket_tick <- 0.035
bar_half_width <- 0.75 / 2
bracket_start <- category_label_positions$x_min - bar_half_width
bracket_end <- category_label_positions$x_max + bar_half_width

# Wide and short, with its legend dropped (the same category/color mapping
# is already shown by the standalone legend, outputs/plots/MM_legend.png).
save_ggplot(
  g1 +
    theme(
      legend.position = "none",
      text = element_text(size = 14),
      # Overrides the size axis.text.x would otherwise inherit from the
      # bumped-up "text" above - at 82 nodes across a 10in-wide, 2.2in-tall
      # panel, that inherited size made the angled node-id labels overlap.
      axis.text.x = element_text(angle = 45, size = 9),
      # Redundant now that every category is separately labeled below (and
      # its default position collides with that new row).
      axis.title.x = element_blank(),
      # Room below the angled node-id row for the bracket bar and the
      # horizontal category-name row added below it (clip = "off" lets that
      # annotation draw past y = 0 into this margin instead of being cut
      # off at the panel border).
      plot.margin = margin(t = 5.5, r = 5.5, b = 42, l = 5.5)
    ) +
    # ylim pinned explicitly (matching the ~5% auto-expansion c(0, 1) would
    # otherwise get by default) rather than left to auto-range: annotate()
    # below adds a layer at y = -0.5, and without an explicit ylim here
    # coord_cartesian() would auto-expand the visible panel to include it,
    # pulling the "below the axis" annotation back inside the panel instead
    # of out into the margin clip = "off" opens up.
    coord_cartesian(ylim = c(-0.05, 1.05), clip = "off") +
    # Horizontal, each centered under its (possibly NSW-merged, see
    # category_label_positions above) group of bars.
    annotate(
      "text",
      x = category_label_positions$x_mid, y = -0.5,
      label = category_label_positions$category, size = 3.2,
      angle = 0, hjust = 0.5, vjust = 1
    ) +
    # Bracket bar under each group (see bracket_start/bracket_end above),
    # sitting in the gap between the angled node-id row above and the
    # category-name row below - the horizontal span plus its two end ticks.
    annotate(
      "segment",
      x = bracket_start, xend = bracket_end, y = bracket_y, yend = bracket_y,
      linewidth = 0.4, color = "grey30"
    ) +
    annotate(
      "segment",
      x = c(bracket_start, bracket_end), xend = c(bracket_start, bracket_end),
      y = bracket_y - bracket_tick, yend = bracket_y + bracket_tick,
      linewidth = 0.4, color = "grey30"
    ),
  "turn_taking",
  REPO_ROOT,
  width = 12,
  height = 2.2,
  subdir = "network_measures"
)

# Fig. 5c: betweenness vs. standardized strength (z_score_incl_self above),
# with marginal histograms, 90th-percentile betweenness line and outlier
# labels. z-scores go negative, so plot_betweenness_strength_scatter()
# min-max normalizes this axis to [0, 1] over its own observed range instead
# of dividing by max (which only works for non-negative values), and its
# marginal fit switches to a normal density instead of an exponential - see
# the "minmax" strength_scale_mode branch in R/graph_utils.R. df3 is built
# directly from V(gra1) attributes, so its row order already matches gra1's
# vertex order, which is what strength_values must align to.
plot_betweenness_strength_scatter(
  gra1,
  plot_path("markov_betweenness_vs_strength_zscore", REPO_ROOT),
  width = 15,
  strength_values = df3$z_score_incl_self,
  strength_label = "standardized strength (z-score vs. shuffle null)",
  show_caption = FALSE,
  text_cex = 1.5,
  axis_padding = 0.05
)

csv_path <- file.path(REPO_ROOT, OUTPUT_DIR, "node_measures.csv")
write.csv(df3, csv_path, row.names = FALSE, quote = FALSE)
message("Saved ", csv_path)

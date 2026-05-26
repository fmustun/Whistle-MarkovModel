#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(scales)
  library(dplyr)
  library(gridExtra)
  library(ggrepel)
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

df <- data.frame(
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
  diff_degree = deg_in - deg_out,
  closeness = nodes_closeness,
  occurrences = V(gra1)$occurrences
)
df$nodes <- as.factor(df$nodes)

df_order <- order_nodes_by_category(df, list_names, "occurrences")
save_ggplot(
  ggplot(df_order) +
    aes(x = nodes, y = occurrences, fill = category) +
    geom_col(width = 0.75, color = "black") +
    scale_fill_manual(breaks = unique(df$category), values = unique(df$color)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45)) +
    labs(title = "Whistle occurrences by node"),
  "whistle_occurrences",
  REPO_ROOT
)

seq_category <- build_seq_category(whistles_list, list_names)

col_same_cat_outbound <- c()
col_diff_cat_outbound <- c()
col_ind_prob_outbound <- c()
col_same_cat_inbound <- c()
col_diff_cat_inbound <- c()
col_ind_prob_inbound <- c()

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
  ind_prob_outbound <- (same_cat_occ_outbound - diff_cat_occ_outbound) /
    (same_cat_occ_outbound + diff_cat_occ_outbound)

  same_cat_occ_inbound <- sum(MarkovModel$transition_counts_matrix[subv, node_i])
  total_occ_inbound <- sum(MarkovModel$transition_counts_matrix[, node_i])
  diff_cat_occ_inbound <- total_occ_inbound - same_cat_occ_inbound
  ind_prob_inbound <- (same_cat_occ_inbound - diff_cat_occ_inbound) /
    (same_cat_occ_inbound + diff_cat_occ_inbound)

  col_same_cat_outbound <- c(col_same_cat_outbound, same_cat_occ_outbound)
  col_diff_cat_outbound <- c(col_diff_cat_outbound, diff_cat_occ_outbound)
  col_ind_prob_outbound <- c(col_ind_prob_outbound, ind_prob_outbound)
  col_same_cat_inbound <- c(col_same_cat_inbound, same_cat_occ_inbound)
  col_diff_cat_inbound <- c(col_diff_cat_inbound, diff_cat_occ_inbound)
  col_ind_prob_inbound <- c(col_ind_prob_inbound, ind_prob_inbound)
}

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

df2_order <- order_nodes_by_category(
  data.frame(
    nodes = V(gra1)$name,
    category = V(gra1)$category,
    color = V(gra1)$color,
    prob_same_cat_outbound = df3$prob_same_cat_outbound
  ),
  list_names,
  "prob_same_cat_outbound"
)
save_ggplot(
  ggplot(df2_order) +
    aes(x = nodes, y = prob_same_cat_outbound, fill = category) +
    geom_col(color = "black", width = 0.75) +
    scale_fill_manual(breaks = unique(df3$category), values = unique(df3$color)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45)) +
    labs(
      title = "Probability of next node in same category",
      y = "probability"
    ) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 0.5),
  "prob_same_category_outbound",
  REPO_ROOT
)

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

g2 <- ggplot(df_order) +
  aes(x = nodes, y = prob_diff_cat_inbound, fill = category) +
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
    title = "P(previous node in a different category)",
    y = "probability"
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.75)

g3 <- ggplot(df_order) +
  aes(x = nodes, y = strength_in, fill = category) +
  geom_col(color = "black", width = 0.75) +
  scale_fill_manual(breaks = unique(df3$category), values = unique(df3$color)) +
  labs(y = "strength in") +
  theme(
    axis.text.x = element_text(angle = 45),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

g4 <- ggplot(df_order) +
  aes(x = nodes, y = betweenness, fill = category) +
  geom_col(color = "black", width = 0.75) +
  scale_fill_manual(breaks = unique(df3$category), values = unique(df3$color)) +
  labs(x = "whistle category") +
  theme(
    axis.text.x = element_text(angle = 45),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

save_ggplot(
  gridExtra::arrangeGrob(g1, g3, g4, ncol = 1),
  "turn_taking_and_centrality",
  REPO_ROOT,
  height = 14
)

save_ggplot(
  ggplot(df3, aes(x = betweenness, y = strength_in)) +
    geom_point(aes(fill = category), colour = "black", pch = 21, size = 5) +
    scale_fill_manual(breaks = unique(df3$category), values = unique(df3$color)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    labs(y = "strength in") +
    geom_text_repel(aes(label = nodes)),
  "betweenness_vs_strength_in",
  REPO_ROOT
)

save_ggplot(
  ggplot(df3, aes(x = prob_diff_cat_outbound, y = betweenness)) +
    geom_point(aes(fill = category), colour = "black", pch = 21, size = 5) +
    scale_fill_manual(breaks = unique(df3$category), values = unique(df3$color)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    geom_text_repel(aes(label = nodes)),
  "prob_diff_outbound_vs_betweenness",
  REPO_ROOT
)

save_ggplot(
  ggplot(df3, aes(x = prob_diff_cat_inbound, y = prob_diff_cat_outbound)) +
    geom_point(aes(fill = category), colour = "black", pch = 21, size = 5) +
    scale_fill_manual(breaks = unique(df3$category), values = unique(df3$color)) +
    geom_segment(aes(x = 0, y = 0.5, xend = 0.5, yend = 0.5), linetype = "dashed") +
    geom_segment(aes(x = 0.5, y = 0, xend = 0.5, yend = 0.5), linetype = "dashed") +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed") +
    labs(
      x = "P(previous node in different category)",
      y = "P(next node in different category)"
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    geom_text_repel(aes(label = nodes)),
  "prob_diff_inbound_vs_outbound",
  REPO_ROOT
)

bet_str <- df3 %>%
  group_by(category) %>%
  summarise(
    betweenness = mean(betweenness),
    strength_in = mean(strength_in),
    .groups = "drop"
  )

save_ggplot(
  ggplot(bet_str, aes(x = betweenness, y = strength_in)) +
    geom_point() +
    theme_minimal() +
    theme(axis.line = element_line(colour = "black")) +
    geom_text_repel(aes(label = category)),
  "category_mean_betweenness_vs_strength_in",
  REPO_ROOT
)

csv_path <- file.path(REPO_ROOT, OUTPUT_DIR, "node_measures.csv")
write.csv(df3, csv_path, row.names = FALSE, quote = FALSE)
message("Saved ", csv_path)

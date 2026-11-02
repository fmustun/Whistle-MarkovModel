#!/usr/bin/env Rscript

# Tests whether a whistle's position within its multi-loop chain
# (compute_multiloop_chains(), R/multiloop_utils.R) is associated with its
# sub-category (whistle_type_chr) - i.e. do some sub-categories tend to lead,
# trail, or sit in the middle of a multi-loop, rather than positions being
# interchangeable across sub-categories?
#
# Position is summarized two ways:
#  - position_category: "first" / "middle" / "last" within the chain. Using
#    the categorical first/middle/last split (rather than the raw
#    position_in_chain integer) avoids confounding position with chain_length
#    - position_in_chain == 4 only exists in chains of length >= 4, so a raw
#    position x sub-category table would partly just reflect which
#    sub-categories happen to occur in longer chains.
#  - relative_position: (position_in_chain - 1) / (chain_length - 1), in
#    [0, 1] (0 = first, 1 = last), comparable across chains of any length.
#
# Whistles within a chain are not independent draws (a chain's sub-category
# composition and length are fixed; only which whistle sits at which
# position is at stake), so a textbook chi-square test's asymptotic p-value
# would not be calibrated here. Instead, both test statistics below are
# assessed against a null built by permuting whistle_type_chr WITHIN each
# chain_id (keeping chain membership, chain length, and each chain's own
# sub-category multiset fixed, and shuffling only which position gets which
# label) - the permutation analogue of "does order matter, given what's in
# the chain", and the same shuffle-null logic already used for transition
# significance elsewhere in this repo (see NULL_MODEL, R/config.R;
# compute_permutation_transition_significance(), R/multiloop_utils.R).
#
# Iteration count: this is a single omnibus test per statistic (not a
# multiple-testing family), so recommended_iterations() (R/statistical_utils.R,
# which sizes iterations from a family size m and BH alpha) does not apply
# here. 9999 permutations gives an empirical p-value floor of 1 / 10000,
# comfortably below the 0.05 significance level used to interpret it, and is
# the conventional default iteration count for a single permutation test.

suppressPackageStartupMessages({
  library(ggplot2)
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
source(file.path(REPO_ROOT, "R", "load_whistles.R"))
source(file.path(REPO_ROOT, "R", "graph_utils.R"))
source(file.path(REPO_ROOT, "R", "multiloop_utils.R"))
source(file.path(REPO_ROOT, "R", "plot_io.R"))

ensure_plots_dir(REPO_ROOT, "multiloops")
multiloop_dir <- file.path(REPO_ROOT, OUTPUT_DIR, "multiloop_analysis")
dir.create(multiloop_dir, recursive = TRUE, showWarnings = FALSE)

iterations_position <- 9999L
min_n_for_plot <- 10

message("Loading whistles from ", DATA_PATH)
whistles_list <- load_whistles(DATA_PATH)

message(
  "Grouping whistles into multi-loop chains (IWI < ",
  MULTILOOP_IWI_THRESHOLD * 1000, "ms)..."
)
whistle_chains <- compute_multiloop_chains(
  whistles_list, threshold = MULTILOOP_IWI_THRESHOLD, list_names = LIST_NAMES
)
ml <- whistle_chains[whistle_chains$chain_length >= 2, ]
ml$position_category <- ifelse(
  ml$position_in_chain == 1L, "first",
  ifelse(ml$position_in_chain == ml$chain_length, "last", "middle")
)
ml$relative_position <- (ml$position_in_chain - 1) / (ml$chain_length - 1)

n_whistles <- nrow(ml)
n_chains <- length(unique(ml$chain_id))
message(
  n_whistles, " whistles across ", n_chains,
  " multi-loop chains (chain length >= 2)"
)

whistles_csv <- file.path(multiloop_dir, "multiloop_position_subcategory_whistles.csv")
write.csv(ml, whistles_csv, row.names = FALSE)
message("Saved ", whistles_csv)

# -- multi-loop vs. isolation proportion by sub-category ---------------------
# Unlike ml above (chain_length >= 2 only), this uses whistle_chains
# unfiltered, since isolated whistles (chain_length == 1) are the other side
# of this comparison.

isolation_sub_levels <- sort(unique(whistle_chains$whistle_type_chr))
isolation_sub_code <- match(whistle_chains$whistle_type_chr, isolation_sub_levels)
in_multiloop <- whistle_chains$chain_length >= 2

n_total_per_isolation_sub <- tabulate(
  isolation_sub_code, nbins = length(isolation_sub_levels)
)
n_multiloop_per_isolation_sub <- tabulate(
  isolation_sub_code[in_multiloop], nbins = length(isolation_sub_levels)
)
n_isolated_per_isolation_sub <-
  n_total_per_isolation_sub - n_multiloop_per_isolation_sub

# node_id matches the vertex numbering scripts/01_run_markov_model.R builds
# the Markov network with (V(gra1)$name <- as.character(whistle_type); see
# compute_graph(), R/graph_utils.R) - i.e. the integer whistle_type code
# whistle_type_chr is derived from, not a re-ranking specific to this script.
node_id_per_isolation_sub <- as.integer(tapply(
  whistle_chains$whistle_type, whistle_chains$whistle_type_chr,
  function(x) x[1]
)[isolation_sub_levels])
category_per_isolation_sub <- tapply(
  whistle_chains$category, whistle_chains$whistle_type_chr, function(x) x[1]
)[isolation_sub_levels]

isolation_summary <- data.frame(
  whistle_type_chr = isolation_sub_levels,
  node_id = node_id_per_isolation_sub,
  category = unname(category_per_isolation_sub),
  n_isolated = n_isolated_per_isolation_sub,
  n_multiloop = n_multiloop_per_isolation_sub,
  n_total = n_total_per_isolation_sub,
  pct_multiloop = round(
    100 * n_multiloop_per_isolation_sub / n_total_per_isolation_sub, 1
  ),
  low_n = n_total_per_isolation_sub < min_n_for_plot,
  stringsAsFactors = FALSE
)
isolation_summary <- isolation_summary[order(isolation_summary$pct_multiloop), ]

n_total_all <- nrow(whistle_chains)
n_multiloop_total <- sum(in_multiloop)
message(
  n_multiloop_total, " of ", n_total_all, " whistles (",
  round(100 * n_multiloop_total / n_total_all, 1),
  "%) occur within a multi-loop chain (chain_length >= 2)"
)

isolation_summary_csv <- file.path(
  multiloop_dir, "multiloop_isolation_subcategory_summary.csv"
)
write.csv(isolation_summary, isolation_summary_csv, row.names = FALSE)
message("Saved ", isolation_summary_csv)

# -- observed statistics ----------------------------------------------------

pos_levels <- c("first", "middle", "last")
sub_levels <- sort(unique(ml$whistle_type_chr))
n_pos <- length(pos_levels)
n_sub <- length(sub_levels)

pos_code <- match(ml$position_category, pos_levels)
sub_code <- match(ml$whistle_type_chr, sub_levels)
chain_code <- as.integer(factor(ml$chain_id))

observed_counts <- tabulate(pos_code + (sub_code - 1L) * n_pos, nbins = n_pos * n_sub)
observed_table <- matrix(
  observed_counts, nrow = n_pos, ncol = n_sub,
  dimnames = list(pos_levels, sub_levels)
)
row_totals <- rowSums(observed_table)
col_totals <- colSums(observed_table)
expected_table <- outer(row_totals, col_totals) / n_whistles
chisq_stat <- function(observed) sum((observed - expected_table)^2 / expected_table)
observed_chisq <- chisq_stat(observed_table)

grand_mean_relpos <- mean(ml$relative_position)
sum_relpos_per_sub <- as.vector(rowsum(ml$relative_position, sub_code))
n_per_sub <- col_totals
mean_relpos_per_sub <- sum_relpos_per_sub / n_per_sub
between_ss <- function(sum_per_group, n_per_group) {
  mean_per_group <- sum_per_group / n_per_group
  sum(n_per_group * (mean_per_group - grand_mean_relpos)^2)
}
observed_between_ss <- between_ss(sum_relpos_per_sub, n_per_sub)

# Standardized (Pearson) residuals from the observed table only - descriptive,
# to flag which sub-category x position cells are over/under-represented;
# not itself the significance test (that's the permutation p-values below).
stdres_table <- suppressWarnings(chisq.test(observed_table, correct = FALSE)$stdres)

# -- permutation null: shuffle whistle_type_chr within chain_id -------------

message(
  "Running ", iterations_position,
  " within-chain permutations for the position x sub-category null model..."
)
sort_idx <- order(chain_code)
chain_code_s <- chain_code[sort_idx]
pos_code_s <- pos_code[sort_idx]
sub_code_s <- sub_code[sort_idx]
relpos_s <- ml$relative_position[sort_idx]

set.seed(SEED)
exceed_chisq <- 0L
exceed_ss <- 0L
progress_every <- 2000L
for (iteration in seq_len(iterations_position)) {
  ord <- order(chain_code_s, runif(n_whistles))
  shuffled_sub_s <- sub_code_s[ord]

  null_counts <- tabulate(
    pos_code_s + (shuffled_sub_s - 1L) * n_pos, nbins = n_pos * n_sub
  )
  null_table <- matrix(null_counts, nrow = n_pos, ncol = n_sub)
  exceed_chisq <- exceed_chisq + (chisq_stat(null_table) >= observed_chisq)

  null_sum_relpos <- as.vector(rowsum(relpos_s, shuffled_sub_s))
  exceed_ss <- exceed_ss +
    (between_ss(null_sum_relpos, n_per_sub) >= observed_between_ss)

  if (iteration %% progress_every == 0L || iteration == iterations_position) {
    message("Permutation pass: ", iteration, " / ", iterations_position)
  }
}
p_chisq <- (1 + exceed_chisq) / (iterations_position + 1)
p_between_ss <- (1 + exceed_ss) / (iterations_position + 1)

message(
  "position_category x sub-category association: chi-sq = ",
  round(observed_chisq, 1), ", permutation p = ", signif(p_chisq, 3)
)
message(
  "relative_position ~ sub-category (between-group SS): stat = ",
  round(observed_between_ss, 3), ", permutation p = ", signif(p_between_ss, 3)
)

# -- per-sub-category descriptive summary -----------------------------------

category_per_sub <- tapply(ml$category, ml$whistle_type_chr, function(x) x[1])
pct_table <- 100 * t(observed_table) / n_per_sub
subcat_summary <- data.frame(
  whistle_type_chr = sub_levels,
  category = unname(category_per_sub[sub_levels]),
  n = as.integer(n_per_sub),
  mean_relative_position = round(mean_relpos_per_sub, 3),
  pct_first = round(pct_table[, "first"], 1),
  pct_middle = round(pct_table[, "middle"], 1),
  pct_last = round(pct_table[, "last"], 1),
  low_n = n_per_sub < min_n_for_plot,
  stringsAsFactors = FALSE
)
subcat_summary <- subcat_summary[order(subcat_summary$mean_relative_position), ]

summary_csv <- file.path(multiloop_dir, "multiloop_position_subcategory_summary.csv")
write.csv(subcat_summary, summary_csv, row.names = FALSE)
message("Saved ", summary_csv)

residual_df <- as.data.frame(as.table(stdres_table))
names(residual_df) <- c("position_category", "whistle_type_chr", "stdres")
residual_df <- residual_df[order(-abs(residual_df$stdres)), ]
residual_csv <- file.path(multiloop_dir, "multiloop_position_subcategory_residuals.csv")
write.csv(residual_df, residual_csv, row.names = FALSE)
message("Saved ", residual_csv)

# -- text summary -------------------------------------------------------------

top_early <- head(subcat_summary[!subcat_summary$low_n, ], 5)
top_late <- tail(subcat_summary[!subcat_summary$low_n, ], 5)
top_residuals <- head(residual_df[abs(residual_df$stdres) >= 2, ], 10)

summary_lines <- c(
  "Multi-loop position x sub-category association",
  paste("IWI threshold:", MULTILOOP_IWI_THRESHOLD * 1000, "ms"),
  paste("whistles in multi-loop chains (chain length >= 2):", n_whistles),
  paste("multi-loop chains:", n_chains),
  paste("distinct sub-categories (whistle_type_chr):", n_sub),
  "",
  "position_category defined as first (position_in_chain == 1), last",
  "(position_in_chain == chain_length), middle (otherwise); relative_position",
  "= (position_in_chain - 1) / (chain_length - 1), in [0, 1].",
  "",
  "Null model: whistle_type_chr shuffled within chain_id only (chain",
  "membership, length, and each chain's own sub-category multiset held",
  "fixed) - isolates whether ORDER within the chain carries information",
  "about identity, as opposed to which sub-categories co-occur.",
  paste0("permutations: ", iterations_position),
  "",
  paste0(
    "Test 1 - position_category (first/middle/last) x sub-category, ",
    "chi-square statistic:"
  ),
  paste0("  observed statistic: ", round(observed_chisq, 1)),
  paste0("  permutation p-value: ", signif(p_chisq, 3)),
  "",
  paste0(
    "Test 2 - relative_position by sub-category, between-group sum of ",
    "squares:"
  ),
  paste0("  observed statistic: ", round(observed_between_ss, 3)),
  paste0("  permutation p-value: ", signif(p_between_ss, 3)),
  "",
  paste0(
    "Descriptive only (not multiple-testing corrected) - sub-categories ",
    "with n >= ", min_n_for_plot, " ranked by mean relative_position:"
  ),
  "  earliest (lowest mean relative_position):",
  paste0(
    "    ", top_early$whistle_type_chr, " (", top_early$category, "): n=",
    top_early$n, ", mean_relative_position=", top_early$mean_relative_position
  ),
  "  latest (highest mean relative_position):",
  paste0(
    "    ", top_late$whistle_type_chr, " (", top_late$category, "): n=",
    top_late$n, ", mean_relative_position=", top_late$mean_relative_position
  ),
  "",
  "Largest standardized residuals (observed table vs. independence,",
  "|stdres| >= 2; descriptive, not permutation-tested individually):",
  paste0(
    "  ", top_residuals$whistle_type_chr, " x ", top_residuals$position_category,
    ": stdres=", round(top_residuals$stdres, 2)
  ),
  "",
  "full per-sub-category breakdown: multiloop_position_subcategory_summary.csv",
  "full residual table: multiloop_position_subcategory_residuals.csv"
)
writeLines(
  summary_lines,
  file.path(multiloop_dir, "multiloop_position_subcategory_summary.txt")
)
message(
  "Saved ",
  file.path(multiloop_dir, "multiloop_position_subcategory_summary.txt")
)

# -- plots --------------------------------------------------------------------

null_chisq_sample <- local({
  set.seed(SEED)
  vapply(seq_len(2000L), function(i) {
    ord <- order(chain_code_s, runif(n_whistles))
    null_counts <- tabulate(
      pos_code_s + (sub_code_s[ord] - 1L) * n_pos, nbins = n_pos * n_sub
    )
    chisq_stat(matrix(null_counts, nrow = n_pos, ncol = n_sub))
  }, numeric(1))
})

panel_null <- ggplot(data.frame(stat = null_chisq_sample), aes(x = stat)) +
  geom_histogram(bins = 40, fill = "#A9A9A9", color = "black") +
  geom_vline(xintercept = observed_chisq, color = "#ee0000", linewidth = 1) +
  theme_minimal() +
  labs(
    title = paste0(
      "Null distribution: position_category x sub-category chi-square ",
      "(within-chain shuffle, 2000 of ", iterations_position, " draws shown)"
    ),
    subtitle = paste0(
      "observed = ", round(observed_chisq, 1),
      " (red line); permutation p = ", signif(p_chisq, 3)
    ),
    x = "chi-square statistic under the null", y = "count"
  )

plot_subcats <- subcat_summary$whistle_type_chr[!subcat_summary$low_n]
order_levels <- subcat_summary$whistle_type_chr[!subcat_summary$low_n]
ml_plot <- ml[ml$whistle_type_chr %in% plot_subcats, ]
ml_plot$whistle_type_chr <- factor(ml_plot$whistle_type_chr, levels = order_levels)

panel_box <- ggplot(
  ml_plot, aes(x = whistle_type_chr, y = relative_position, fill = category)
) +
  geom_boxplot(outlier.size = 0.6) +
  coord_flip() +
  scale_fill_manual(breaks = LIST_NAMES, values = LIST_COLORS[seq_along(LIST_NAMES)]) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6)) +
  labs(
    title = paste0(
      "Relative position in multi-loop by sub-category (n >= ",
      min_n_for_plot, "; ordered by mean)"
    ),
    subtitle = paste0(
      "0 = first in chain, 1 = last; between-group SS permutation p = ",
      signif(p_between_ss, 3)
    ),
    x = NULL, y = "relative position", fill = "category"
  )

pct_long <- data.frame(
  whistle_type_chr = factor(rep(sub_levels, n_pos), levels = order_levels),
  position_category = factor(rep(pos_levels, each = n_sub), levels = pos_levels),
  pct = as.vector(t(pct_table))
)
pct_long <- pct_long[pct_long$whistle_type_chr %in% plot_subcats, ]

panel_bar <- ggplot(
  pct_long, aes(x = whistle_type_chr, y = pct, fill = position_category)
) +
  geom_col(position = "fill", color = "black", width = 0.75) +
  coord_flip() +
  scale_fill_manual(
    breaks = pos_levels, values = c("#61D04F", "#A9A9A9", "#2297E6")
  ) +
  theme_minimal() +
  theme(axis.text.y = element_blank()) +
  labs(
    title = "Position-in-chain composition by sub-category",
    x = NULL, y = "proportion of occurrences", fill = "position in chain"
  )

layout_matrix <- rbind(c(1, 1), c(2, 3))
combined_plot <- gridExtra::arrangeGrob(
  panel_null, panel_box, panel_bar, layout_matrix = layout_matrix,
  widths = c(1.3, 1)
)
save_ggplot(
  combined_plot, "multiloop_position_subcategory", REPO_ROOT,
  width = 13, height = 15, subdir = "multiloops"
)

isolation_order <- isolation_summary$whistle_type_chr[!isolation_summary$low_n]
isolation_plotted <- isolation_summary[isolation_summary$whistle_type_chr %in% isolation_order, ]
# Axis breaks are still keyed by whistle_type_chr (unique per node_id) so
# the fill/color layers and axis stay in sync; only the printed tick text
# is swapped for node_id via scale_x_discrete(labels = ...) below.
isolation_node_labels <- isolation_plotted$node_id[
  match(isolation_order, isolation_plotted$whistle_type_chr)
]

isolation_long <- data.frame(
  whistle_type_chr = factor(
    rep(isolation_summary$whistle_type_chr, 2), levels = isolation_order
  ),
  status = factor(
    rep(c("multi-loop", "isolated"), each = nrow(isolation_summary)),
    levels = c("isolated", "multi-loop")
  ),
  n = c(isolation_summary$n_multiloop, isolation_summary$n_isolated)
)
isolation_long <- isolation_long[isolation_long$whistle_type_chr %in% isolation_order, ]

panel_isolation <- ggplot(
  isolation_long, aes(x = whistle_type_chr, y = n, fill = status)
) +
  geom_col(position = "fill", color = "black", width = 0.75) +
  # Small square past the bar's end flags each node's main category -
  # a distinct `color` scale from the bars' `fill`, so both keep their own
  # legend without needing a third-party multi-scale package.
  geom_point(
    data = isolation_plotted,
    aes(x = whistle_type_chr, y = 1.05, color = category),
    inherit.aes = FALSE, shape = 15, size = 1.8
  ) +
  # ylim (not scale_y_continuous(limits=...)) so the extra room for the
  # category marker only clips the viewport - scale limits would instead
  # censor raw counts to NA before geom_col's position="fill" normalizes
  # them, silently dropping almost every bar.
  coord_flip(clip = "off", ylim = c(0, 1.08)) +
  scale_fill_manual(
    name = NULL, breaks = c("isolated", "multi-loop"),
    values = c("#A9A9A9", "#2297E6")
  ) +
  scale_color_manual(
    name = "category", breaks = LIST_NAMES,
    values = LIST_COLORS[seq_along(LIST_NAMES)]
  ) +
  scale_x_discrete(breaks = isolation_order, labels = isolation_node_labels) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6)) +
  labs(
    title = paste0(
      "Multi-loop vs. isolated occurrence by sub-category (n >= ",
      min_n_for_plot, "; ordered by % multi-loop)"
    ),
    subtitle = paste0(
      n_multiloop_total, " of ", n_total_all, " whistles overall (",
      round(100 * n_multiloop_total / n_total_all, 1), "%) occur in a chain; ",
      "y-axis = Markov model node number; square = main category"
    ),
    x = NULL, y = "proportion of occurrences", fill = NULL
  )

save_ggplot(
  panel_isolation, "multiloop_isolation_subcategory", REPO_ROOT,
  width = 10, height = 11, subdir = "multiloops"
)

# -- multi-loop rate vs. occurrence frequency scatter -------------------------
# Same n >= min_n_for_plot restriction as panel_isolation above (proportion
# is unstable at low n); node_id/category match the Markov model node
# numbering scripts/01_run_markov_model.R builds the network with (see the
# node_id_per_isolation_sub comment above).

scatter_data <- isolation_summary[!isolation_summary$low_n, ]
overall_pct_multiloop <- 100 * n_multiloop_total / n_total_all

panel_occurrence_scatter <- ggplot(
  scatter_data, aes(x = n_total, y = pct_multiloop)
) +
  geom_hline(
    yintercept = overall_pct_multiloop, linetype = "dashed", color = "grey40"
  ) +
  geom_point(aes(fill = category), shape = 21, color = "black", size = 2.5) +
  geom_text_repel(
    aes(label = node_id), size = 2.2, max.overlaps = Inf,
    segment.size = 0.2, segment.color = "grey50"
  ) +
  scale_x_log10() +
  scale_fill_manual(
    name = "category", breaks = LIST_NAMES,
    values = LIST_COLORS[seq_along(LIST_NAMES)]
  ) +
  theme_minimal() +
  labs(
    title = paste0(
      "Multi-loop rate vs. occurrence frequency by sub-category (n >= ",
      min_n_for_plot, ")"
    ),
    subtitle = paste0(
      "dashed line = overall rate (", round(overall_pct_multiloop, 1),
      "%); label = Markov model node number; point color = main category"
    ),
    x = "number of occurrences (log scale)",
    y = "% of occurrences in a multi-loop chain"
  )

save_ggplot(
  panel_occurrence_scatter, "multiloop_occurrence_rate_scatter", REPO_ROOT,
  width = 9, height = 7, subdir = "multiloops"
)

# -- same rate, against whistle-level network measures instead of frequency --
# turn-taking (P(next node in a different category)), betweenness, and
# standardized strength (z-score vs. shuffle null, self-transitions
# included) all come from the BH-significant whistle-level network
# (scripts/01_run_markov_model.R, scripts/02_node_network_measures.R),
# joined here by node_id - the same Markov model node numbering used
# throughout this script. Unlike occurrence frequency, none of these three
# are meaningful on a log scale (betweenness is routinely 0; the z-score can
# be negative), so all three keep a linear x-axis.

node_measures_csv <- file.path(REPO_ROOT, OUTPUT_DIR, "node_measures.csv")
if (!file.exists(node_measures_csv)) {
  stop(
    "Run scripts/02_node_network_measures.R first. Missing: ", node_measures_csv,
    call. = FALSE
  )
}
node_measures <- read.csv(node_measures_csv, stringsAsFactors = FALSE)
node_measures$node_id <- as.integer(node_measures$nodes)

network_measure_data <- merge(
  scatter_data,
  node_measures[, c("node_id", "betweenness", "prob_diff_cat_outbound", "z_score_incl_self")],
  by = "node_id"
)
message(
  nrow(scatter_data) - nrow(network_measure_data), " / ", nrow(scatter_data),
  " sub-categories (n >= ", min_n_for_plot, ") have no significant transition in the ",
  "whistle-level network (scripts/01/02) and are excluded from the network-measure scatter plots"
)

plot_occurrence_rate_vs_measure <- function(x_col, x_label, name) {
  save_ggplot(
    ggplot(network_measure_data, aes(x = .data[[x_col]], y = pct_multiloop)) +
      geom_hline(
        yintercept = overall_pct_multiloop, linetype = "dashed", color = "grey40"
      ) +
      geom_point(aes(fill = category), shape = 21, color = "black", size = 2.5) +
      geom_text_repel(
        aes(label = node_id), size = 2.2, max.overlaps = Inf,
        segment.size = 0.2, segment.color = "grey50"
      ) +
      scale_fill_manual(
        name = "category", breaks = LIST_NAMES,
        values = LIST_COLORS[seq_along(LIST_NAMES)]
      ) +
      theme_minimal() +
      labs(
        title = paste0(
          "Multi-loop rate vs. ", x_label, " by sub-category (n >= ",
          min_n_for_plot, ")"
        ),
        subtitle = paste0(
          "dashed line = overall rate (", round(overall_pct_multiloop, 1),
          "%); label = Markov model node number; point color = main category; ",
          "restricted to sub-categories in the BH-significant network (scripts/01/02)"
        ),
        x = x_label,
        y = "% of occurrences in a multi-loop chain"
      ),
    name,
    REPO_ROOT,
    width = 9, height = 7, subdir = "multiloops"
  )
}

plot_occurrence_rate_vs_measure(
  "prob_diff_cat_outbound", "turn-taking (P(next node in a different category))",
  "multiloop_occurrence_rate_vs_turn_taking"
)
plot_occurrence_rate_vs_measure(
  "betweenness", "betweenness",
  "multiloop_occurrence_rate_vs_betweenness"
)
plot_occurrence_rate_vs_measure(
  "z_score_incl_self", "standardized strength (z-score vs. shuffle null)",
  "multiloop_occurrence_rate_vs_strength_zscore"
)

# -- same, but P(first|in chain) / P(last|in chain) vs. chance expectation --
# "% of occurrences in the first/last position" (n_first / n_total) is
# capped by pct_multiloop, so a low rate can just mean the sub-category is
# rarely in a chain at all - it confounds "how often in a chain" with
# "where in the chain". Conditioning on being in a chain (denominator
# n_multiloop instead of n_total) removes that: P(first | in chain) =
# n_first / n_multiloop, P(last | in chain) = n_last / n_multiloop.
#
# Even so, these aren't comparable to a single flat rate across
# sub-categories, because chance itself depends on each sub-category's own
# chain-length mix: in a chain of length L, only 1 of L slots is "first"
# (same for "last"), so a sub-category that mostly sits in 2-whistle chains
# has a much higher chance level than one that mostly sits in 8-whistle
# chains, with no position preference at all. The chance level is exactly
# the expectation of this statistic under the within-chain label-shuffle
# null already used above (see "permutation null: shuffle whistle_type_chr
# within chain_id"): a uniformly random permutation of a chain's L labels
# sends any given token to the first slot (or the last slot) with
# probability 1/L, by symmetry of a random permutation over L items -
# regardless of how many chains the sub-category appears in or how many
# times it repeats within one chain. So the chance level for a
# sub-category is the average of 1/chain_length over its own in-chain
# occurrences (an exact analytic expectation, not a Monte Carlo estimate).

first_code <- match(
  ml$whistle_type_chr[ml$position_category == "first"], isolation_sub_levels
)
last_code <- match(
  ml$whistle_type_chr[ml$position_category == "last"], isolation_sub_levels
)
sum_inv_chain_length <- tapply(1 / ml$chain_length, ml$whistle_type_chr, sum)
sum_inv_chain_length <- unname(sum_inv_chain_length[isolation_sub_levels])
sum_inv_chain_length[is.na(sum_inv_chain_length)] <- 0

position_end_counts <- data.frame(
  whistle_type_chr = isolation_sub_levels,
  n_first = tabulate(first_code, nbins = length(isolation_sub_levels)),
  n_last = tabulate(last_code, nbins = length(isolation_sub_levels)),
  sum_inv_chain_length = sum_inv_chain_length,
  stringsAsFactors = FALSE
)

scatter_data_positions <- merge(scatter_data, position_end_counts, by = "whistle_type_chr")

# n_total >= min_n_for_plot (already true of scatter_data) no longer implies
# a stable in-chain proportion, since n_multiloop can be much smaller than
# n_total - so apply the same minimum sample-size threshold to the actual
# denominator these conditional probabilities are estimated from.
scatter_data_positions <- scatter_data_positions[
  scatter_data_positions$n_multiloop >= min_n_for_plot,
]
scatter_data_positions$p_first_given_chain <- round(
  scatter_data_positions$n_first / scatter_data_positions$n_multiloop, 3
)
scatter_data_positions$p_last_given_chain <- round(
  scatter_data_positions$n_last / scatter_data_positions$n_multiloop, 3
)
scatter_data_positions$p_chance <- round(
  scatter_data_positions$sum_inv_chain_length / scatter_data_positions$n_multiloop, 3
)

network_measure_data_positions <- merge(
  scatter_data_positions,
  node_measures[, c("node_id", "betweenness", "prob_diff_cat_outbound", "z_score_incl_self")],
  by = "node_id"
)

# turn-taking (prob_diff_cat_outbound) is already a probability, in [0, 1].
# betweenness (>= 0, routinely 0, long right tail) and the strength z-score
# (can be negative) are not, so min-max scale each to [0, 1] over this same
# set of sub-categories - a linear rescale (preserves relative spacing,
# unlike e.g. a rank/percentile transform) so it changes only the axis
# units, not the shape of the relationship being plotted.
min_max_scale <- function(x) (x - min(x)) / (max(x) - min(x))
network_measure_data_positions$betweenness_scaled <-
  min_max_scale(network_measure_data_positions$betweenness)
network_measure_data_positions$z_score_scaled <-
  min_max_scale(network_measure_data_positions$z_score_incl_self)

position_chance_csv <- file.path(multiloop_dir, "multiloop_position_chance_summary.csv")
write.csv(
  scatter_data_positions[, c(
    "whistle_type_chr", "node_id", "category", "n_total", "n_multiloop",
    "n_first", "n_last", "p_first_given_chain", "p_last_given_chain", "p_chance"
  )],
  position_chance_csv, row.names = FALSE
)
message("Saved ", position_chance_csv)

positions_subdir <- file.path("multiloops", "positions_vs_network_measure")
ensure_plots_dir(REPO_ROOT, positions_subdir)

# Long form: one row per (sub-category, position) so first and last can
# share a single plot as two facet panels instead of two separate files.
long_position_data <- rbind(
  transform(
    network_measure_data_positions,
    position = factor("first", levels = c("first", "last")),
    p_observed = p_first_given_chain
  ),
  transform(
    network_measure_data_positions,
    position = factor("last", levels = c("first", "last")),
    p_observed = p_last_given_chain
  )
)

position_measures <- list(
  list(
    col = "prob_diff_cat_outbound",
    label = "turn-taking (P(next node in a different category))",
    suffix = "turn_taking"
  ),
  list(
    col = "betweenness_scaled",
    label = "betweenness (normalized)",
    suffix = "betweenness"
  ),
  list(
    col = "z_score_scaled",
    label = "standardized strength (z-score vs. shuffle null, normalized)",
    suffix = "strength_zscore"
  )
)

# -- Spearman correlation: position rate vs. each network measure -----------
# Rank-based, so it's unaffected by the min-max scaling of betweenness/
# z-score above (a monotonic transform leaves ranks unchanged) - the same
# rho would come out of the raw, unscaled measure. Ties (e.g. betweenness
# is routinely 0) mean the exact null distribution doesn't apply, so
# cor.test() falls back to its asymptotic p-value approximation; this is
# the standard behavior for Spearman with ties, not a special case here.
position_correlations <- do.call(rbind, lapply(position_measures, function(measure) {
  do.call(rbind, lapply(c("first", "last"), function(pos) {
    subset_data <- long_position_data[long_position_data$position == pos, ]
    test <- suppressWarnings(cor.test(
      subset_data[[measure$col]], subset_data$p_observed, method = "spearman"
    ))
    data.frame(
      measure = measure$suffix,
      measure_label = measure$label,
      position = pos,
      n = nrow(subset_data),
      rho = round(unname(test$estimate), 3),
      p_value = signif(unname(test$p.value), 3),
      stringsAsFactors = FALSE
    )
  }))
}))

# BH correction across all 6 tests as one family (3 network measures x
# first/last), the same stats::p.adjust(method = "BH") used for the
# transition-significance network in scripts/01_run_markov_model.R
# (adjust_transition_p_values(), R/statistical_utils.R) - these 6 are
# examined together here, so that's the family a reader is implicitly
# asking "how many of these are real" over, not each test in isolation.
position_correlations$p_value_bh <- signif(
  stats::p.adjust(position_correlations$p_value, method = "BH"), 3
)
position_correlations$significant_bh <- position_correlations$p_value_bh < GRAPH_ALPHA

position_correlations_csv <- file.path(
  multiloop_dir, "multiloop_position_correlation_summary.csv"
)
write.csv(position_correlations, position_correlations_csv, row.names = FALSE)
message("Saved ", position_correlations_csv)
for (i in seq_len(nrow(position_correlations))) {
  row <- position_correlations[i, ]
  message(
    "Spearman correlation, P(", row$position, " | in loop) vs. ", row$measure_label,
    ": rho = ", row$rho, ", p = ", row$p_value, ", BH-adjusted p = ", row$p_value_bh,
    " (n = ", row$n, ")"
  )
}

format_p <- function(p) {
  ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, format = "f", digits = 3)))
}
position_correlations$annotation <- paste0(
  "rho = ", formatC(position_correlations$rho, format = "f", digits = 2),
  ", BH-adj. ", format_p(position_correlations$p_value_bh)
)

# Point = observed P(position | in loop); "x" mark + connecting segment =
# that sub-category's own chance level (see derivation above) - so each
# point is compared directly to its own chance expectation rather than to a
# single population-wide reference line. Per-facet text = Spearman
# correlation between the observed rate and the x-axis measure, computed
# above (position_correlations).
plot_position_rate_vs_measure <- function(x_col, x_label, name, measure_suffix) {
  ann <- position_correlations[position_correlations$measure == measure_suffix, ]
  save_ggplot(
    ggplot(long_position_data, aes(x = .data[[x_col]])) +
      geom_segment(
        aes(xend = .data[[x_col]], y = p_chance, yend = p_observed),
        color = "grey65", linewidth = 0.4
      ) +
      geom_point(
        aes(y = p_chance), shape = 4, color = "grey40", size = 1.6, stroke = 0.8
      ) +
      geom_point(aes(y = p_observed, fill = category), shape = 21, color = "black", size = 2.5) +
      geom_text_repel(
        aes(y = p_observed, label = node_id), size = 3.2, max.overlaps = Inf,
        segment.size = 0.2, segment.color = "grey50"
      ) +
      geom_text(
        data = ann, aes(x = -Inf, y = Inf, label = annotation),
        inherit.aes = FALSE, hjust = -0.08, vjust = 1.6, size = 4.5, fontface = "italic"
      ) +
      facet_wrap(~position, ncol = 2) +
      scale_fill_manual(
        name = "category", breaks = LIST_NAMES, values = LIST_COLORS[seq_along(LIST_NAMES)]
      ) +
      theme_minimal(base_size = 15) +
      theme(
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 16),
        strip.background = element_rect(fill = "grey85", color = "grey50"),
        panel.spacing = unit(3, "lines"),
        panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.6),
        plot.title = element_text(size = 17)
      ) +
      labs(
        title = paste0("P(position | in loop) vs. ", x_label, " by sub-category"),
        x = x_label,
        y = "P(position | in loop)"
      ),
    name,
    REPO_ROOT,
    width = 14, height = 7, subdir = positions_subdir
  )
}

for (measure in position_measures) {
  plot_position_rate_vs_measure(
    measure$col, measure$label,
    paste0("multiloop_position_rate_vs_", measure$suffix),
    measure$suffix
  )
}

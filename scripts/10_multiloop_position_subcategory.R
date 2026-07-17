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

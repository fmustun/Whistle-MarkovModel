#!/usr/bin/env Rscript

# Tests whether individual whistle sub-categories (whistle_type_chr) are
# over- or under-represented among whistles that end up inside a multi-loop
# chain (compute_multiloop_chains(), chain_length >= 2 - "in a loop"),
# compared to whistles that don't join one (chain_length == 1 - the same
# "singleton" whistles script 04 already reports; here called "inter-event"
# to match how the user framed the question: not part of a tight IWI burst,
# but still an individual whistle in the recording's overall sequence).
# This is a whistle-level loop/not-loop split - a different unit from
# script 09's "inter-event" Markov TRANSITION instances (a pair of
# whistles not sharing a chain); those aren't the same thing, so results
# from the two scripts describe different questions and shouldn't be read
# as directly comparable.
#
# Two results per sub-category, both from one shared permutation null:
#  - a single omnibus test: is sub-category identity associated with
#    loop/inter-event status at all, across every sub-category jointly
#    (chi-square statistic on the 2 x n_sub table).
#  - a per-sub-category enrichment test: for THIS sub-category specifically,
#    is its count among loop whistles higher (or lower) than a shuffle-null
#    would produce, BH-corrected across the n_sub sub-categories tested
#    (the same restricted-family idea as elsewhere in this repo - see
#    R/statistical_utils.R).
#
# Null model: whistle_type_chr shuffled globally across every whistle
# (NULL_MODEL = "shuffle", R/config.R - the same null MarkovWhistle() uses).
# Unlike scripts/10's within-chain shuffle, no chain-level restriction is
# needed here: loop/inter-event status is fixed by timing alone (chain_id
# assignment never depends on the label), so a whistle's chain membership
# is exactly the same under any relabeling - only which category ends up
# "in a loop how often" is being tested.
#
# Iterations: this per-sub-category test IS a multiple-testing family (one
# hypothesis per sub-category, m = n_sub), so recommended_iterations()
# (R/statistical_utils.R) applies directly and is used to derive the
# permutation count from m and GRAPH_ALPHA, rather than picking a number
# by hand.

suppressPackageStartupMessages({
  library(ggplot2)
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
source(file.path(REPO_ROOT, "R", "load_whistles.R"))
source(file.path(REPO_ROOT, "R", "graph_utils.R"))
source(file.path(REPO_ROOT, "R", "multiloop_utils.R"))
source(file.path(REPO_ROOT, "R", "statistical_utils.R"))
source(file.path(REPO_ROOT, "R", "plot_io.R"))

ensure_plots_dir(REPO_ROOT, "multiloops")
multiloop_dir <- file.path(REPO_ROOT, OUTPUT_DIR, "multiloop_analysis")
dir.create(multiloop_dir, recursive = TRUE, showWarnings = FALSE)

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
is_loop <- whistle_chains$chain_length >= 2
n_total <- nrow(whistle_chains)
n_loop <- sum(is_loop)
n_interevent <- n_total - n_loop
message(
  n_loop, " / ", n_total, " whistles (", round(100 * n_loop / n_total, 1),
  "%) are in a multi-loop chain (chain_length >= 2); ", n_interevent,
  " are inter-event (chain_length == 1)"
)

# -- observed statistics ------------------------------------------------------

sub_levels <- sort(unique(whistle_chains$whistle_type_chr))
n_sub <- length(sub_levels)
sub_code <- match(whistle_chains$whistle_type_chr, sub_levels)
loop_idx <- which(is_loop)

observed_loop_counts <- tabulate(sub_code[loop_idx], nbins = n_sub)
observed_total_counts <- tabulate(sub_code, nbins = n_sub)
observed_interevent_counts <- observed_total_counts - observed_loop_counts

observed_table <- rbind(observed_loop_counts, observed_interevent_counts)
dimnames(observed_table) <- list(c("loop", "inter_event"), sub_levels)
row_totals <- c(n_loop, n_interevent)
expected_table <- outer(row_totals, observed_total_counts) / n_total
chisq_stat <- function(observed) sum((observed - expected_table)^2 / expected_table)
observed_chisq <- chisq_stat(observed_table)

overall_loop_rate <- n_loop / n_total
category_loop_rate <- observed_loop_counts / observed_total_counts
enrichment_ratio <- category_loop_rate / overall_loop_rate

# -- permutation null: shuffle whistle_type_chr globally ---------------------

m <- n_sub
iterations_enrichment <- recommended_iterations(m, GRAPH_ALPHA, margin = 5)
message(
  "Running ", iterations_enrichment, " global-shuffle permutations (m = ",
  m, " sub-categories, GRAPH_ALPHA = ", GRAPH_ALPHA, ")..."
)

set.seed(SEED)
exceed_chisq <- 0L
exceed_more <- integer(n_sub)
exceed_less <- integer(n_sub)
progress_every <- 2000L
for (iteration in seq_len(iterations_enrichment)) {
  shuffled_sub_code <- sample(sub_code)
  null_loop_counts <- tabulate(shuffled_sub_code[loop_idx], nbins = n_sub)
  null_total_counts <- tabulate(shuffled_sub_code, nbins = n_sub)
  null_table <- rbind(null_loop_counts, null_total_counts - null_loop_counts)

  exceed_chisq <- exceed_chisq + (chisq_stat(null_table) >= observed_chisq)
  exceed_more <- exceed_more + (null_loop_counts >= observed_loop_counts)
  exceed_less <- exceed_less + (null_loop_counts <= observed_loop_counts)

  if (iteration %% progress_every == 0L || iteration == iterations_enrichment) {
    message("Permutation pass: ", iteration, " / ", iterations_enrichment)
  }
}

p_chisq <- (1 + exceed_chisq) / (iterations_enrichment + 1)
p_more_in_loop <- (1 + exceed_more) / (iterations_enrichment + 1)
p_more_interevent <- (1 + exceed_less) / (iterations_enrichment + 1)
p_two_sided <- pmin(1, 2 * pmin(p_more_in_loop, p_more_interevent))
bh_q <- stats::p.adjust(p_two_sided, method = "BH")

message(
  "Omnibus (loop/inter-event x sub-category) chi-sq = ",
  round(observed_chisq, 1), ", permutation p = ", signif(p_chisq, 3)
)
n_significant <- sum(bh_q < GRAPH_ALPHA)
message(
  n_significant, " / ", n_sub,
  " sub-categories significantly enriched in one direction (BH q < ",
  GRAPH_ALPHA, ")"
)

# -- per-sub-category table ---------------------------------------------------

category_of_sub <- tapply(
  whistle_chains$category, whistle_chains$whistle_type_chr, function(x) x[1]
)
direction <- ifelse(
  bh_q >= GRAPH_ALPHA, "not_significant",
  ifelse(category_loop_rate > overall_loop_rate, "more_in_loops", "more_inter_event")
)

enrichment_table <- data.frame(
  whistle_type_chr = sub_levels,
  category = unname(category_of_sub[sub_levels]),
  n_total = observed_total_counts,
  n_loop = observed_loop_counts,
  n_inter_event = observed_interevent_counts,
  pct_in_loop = round(100 * category_loop_rate, 1),
  enrichment_ratio = round(enrichment_ratio, 3),
  p_more_in_loop = signif(p_more_in_loop, 3),
  p_more_inter_event = signif(p_more_interevent, 3),
  bh_q = signif(bh_q, 3),
  direction = direction,
  stringsAsFactors = FALSE
)
enrichment_table <- enrichment_table[order(enrichment_table$bh_q, -enrichment_table$enrichment_ratio), ]

enrichment_csv <- file.path(multiloop_dir, "multiloop_subcategory_enrichment.csv")
write.csv(enrichment_table, enrichment_csv, row.names = FALSE)
message("Saved ", enrichment_csv)

# -- text summary -------------------------------------------------------------

more_in_loops <- enrichment_table[enrichment_table$direction == "more_in_loops", ]
more_interevent <- enrichment_table[enrichment_table$direction == "more_inter_event", ]

summary_lines <- c(
  "Multi-loop vs. inter-event sub-category enrichment",
  paste("IWI threshold:", MULTILOOP_IWI_THRESHOLD * 1000, "ms"),
  paste0(
    "whistles: ", n_total, " total, ", n_loop, " in a multi-loop chain (",
    round(100 * overall_loop_rate, 1), "%), ", n_interevent, " inter-event (singleton)"
  ),
  paste("distinct sub-categories (whistle_type_chr):", n_sub),
  "",
  "Null model: whistle_type_chr shuffled globally across all whistles",
  "(loop/inter-event status is fixed by timing alone, unaffected by the",
  "shuffle) - tests whether sub-category identity predicts loop membership.",
  paste0("permutations: ", iterations_enrichment,
         " (= recommended_iterations(m = ", m, ", alpha = ", GRAPH_ALPHA, ", margin = 5))"),
  "",
  "Omnibus test - loop/inter-event x sub-category, chi-square statistic:",
  paste0("  observed statistic: ", round(observed_chisq, 1)),
  paste0("  permutation p-value: ", signif(p_chisq, 3)),
  "",
  paste0(
    "Per-sub-category enrichment (two-sided, BH-adjusted across all ",
    n_sub, " sub-categories, alpha = ", GRAPH_ALPHA, "):"
  ),
  paste0(
    "  ", n_significant, " / ", n_sub, " sub-categories significant (",
    nrow(more_in_loops), " more in loops, ", nrow(more_interevent), " more inter-event)"
  ),
  "  significantly MORE in loops than expected:",
  if (nrow(more_in_loops) > 0) {
    paste0(
      "    ", more_in_loops$whistle_type_chr, " (", more_in_loops$category,
      "): ", more_in_loops$n_loop, "/", more_in_loops$n_total, " in loops (",
      more_in_loops$pct_in_loop, "%), enrichment ratio=", more_in_loops$enrichment_ratio,
      ", q=", more_in_loops$bh_q
    )
  } else {
    "    (none)"
  },
  "  significantly MORE inter-event (LESS in loops) than expected:",
  if (nrow(more_interevent) > 0) {
    paste0(
      "    ", more_interevent$whistle_type_chr, " (", more_interevent$category,
      "): ", more_interevent$n_loop, "/", more_interevent$n_total, " in loops (",
      more_interevent$pct_in_loop, "%), enrichment ratio=", more_interevent$enrichment_ratio,
      ", q=", more_interevent$bh_q
    )
  } else {
    "    (none)"
  },
  "",
  "full per-sub-category breakdown: multiloop_subcategory_enrichment.csv"
)
writeLines(
  summary_lines,
  file.path(multiloop_dir, "multiloop_subcategory_enrichment_summary.txt")
)
message(
  "Saved ", file.path(multiloop_dir, "multiloop_subcategory_enrichment_summary.txt")
)

# -- plots --------------------------------------------------------------------

null_chisq_sample <- local({
  set.seed(SEED)
  vapply(seq_len(2000L), function(i) {
    shuffled <- sample(sub_code)
    null_loop <- tabulate(shuffled[loop_idx], nbins = n_sub)
    null_tot <- tabulate(shuffled, nbins = n_sub)
    chisq_stat(rbind(null_loop, null_tot - null_loop))
  }, numeric(1))
})

panel_null <- ggplot(data.frame(stat = null_chisq_sample), aes(x = stat)) +
  geom_histogram(bins = 40, fill = "#A9A9A9", color = "black") +
  geom_vline(xintercept = observed_chisq, color = "#ee0000", linewidth = 1) +
  theme_minimal() +
  labs(
    title = paste0(
      "Null distribution: loop x sub-category chi-square (global shuffle, ",
      "2000 of ", iterations_enrichment, " draws shown)"
    ),
    subtitle = paste0(
      "observed = ", round(observed_chisq, 1),
      " (red line); permutation p = ", signif(p_chisq, 3)
    ),
    x = "chi-square statistic under the null", y = "count"
  )

volcano_df <- enrichment_table[enrichment_table$n_total >= min_n_for_plot, ]
# log2(0) is undefined for a sub-category never seen in a loop; floored to
# just below the smallest non-zero ratio actually observed so it still
# plots (at the far left) instead of being dropped.
floor_ratio <- min(volcano_df$enrichment_ratio[volcano_df$enrichment_ratio > 0]) / 2
volcano_df$log2_ratio <- log2(pmax(volcano_df$enrichment_ratio, floor_ratio))
volcano_df$neg_log10_q <- -log10(pmax(volcano_df$bh_q, 1e-300))

panel_volcano <- ggplot(
  volcano_df, aes(x = log2_ratio, y = neg_log10_q, fill = category)
) +
  geom_hline(yintercept = -log10(GRAPH_ALPHA), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(shape = 21, size = 2.6, color = "black") +
  ggrepel::geom_text_repel(
    data = volcano_df[volcano_df$direction != "not_significant", ],
    aes(label = whistle_type_chr), size = 2.6, max.overlaps = 30
  ) +
  scale_fill_manual(breaks = LIST_NAMES, values = LIST_COLORS[seq_along(LIST_NAMES)]) +
  theme_minimal() +
  labs(
    title = paste0(
      "Sub-category enrichment in multi-loop chains vs. inter-event (n >= ",
      min_n_for_plot, ")"
    ),
    subtitle = paste0(
      "dashed horizontal line = BH q = ", GRAPH_ALPHA,
      "; labeled points are significant; left = more inter-event, right = more in loops"
    ),
    x = "log2(enrichment ratio): loop rate / overall loop rate",
    y = "-log10(BH q)",
    fill = "category"
  )

layout_matrix <- rbind(c(1), c(2))
combined_plot <- gridExtra::arrangeGrob(
  panel_null, panel_volcano, layout_matrix = layout_matrix, heights = c(1, 1.6)
)
save_ggplot(
  combined_plot, "multiloop_subcategory_enrichment", REPO_ROOT,
  width = 10, height = 12, subdir = "multiloops"
)

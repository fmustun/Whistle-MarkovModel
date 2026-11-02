#!/usr/bin/env Rscript

# Regenerates panels b-e of outputs/plots/SFigure_null_model.png (the code
# that produced that PNG isn't in this repo - only the rendered image was
# committed). Panels a and f of that figure are out of scope here (interval
# distribution fit and the full network plot, respectively); this script
# only rebuilds the four example null-distribution histograms. Output goes
# to outputs/plots/null_model_examples/: the combined 4-panel column (as in
# the original figure), each panel as its own file, and a stats CSV.
#
# Each panel shows the shuffle-null distribution of one transition
# probability P(from -> to) (MarkovModel$random_whistle_occurrences_matrices,
# R/markov_whistle.R) against the observed probability (red dashed line),
# and is annotated with the current BH-corrected q-value
# (outputs/bh_analysis/bh_q_values.rds) rather than the uncorrected
# direct-randomization p the original legend quoted - see
# [[project-bh-correction-methodology]].
#
# Pair selection: b (8->2) is unchanged from the original figure and is
# still BH-significant. c's original pair (16->26) was dropped: its raw
# empirical p under the current 150,000-iteration run (~0.0023) doesn't
# match the original legend's p=0.004, not because anything is wrong but
# because the two numbers were estimated with very different Monte Carlo
# precision (git history shows ITERATIONS was 1,000 until commit d82eaaf,
# 2026-07-14, when it became 150,000 as part of the BH-correction work - the
# original figure predates that, given its legend only quotes a raw
# direct-randomization p with no BH/q-value language). Rather than carry
# that explanation into the figure, c now uses 40->37: BH-significant
# (q ~= 0.0005) with a similarly sparse null distribution (23.9x spike/rest
# ratio, needing the same axis break as e below).
#
# d's original pair (50->21) is no longer usable - it has zero observed
# occurrences in the current outputs/markov_model.rds (count = 0), so it's
# trivially excluded from the BH hypothesis family (family_mask requires
# observed_probability > 0) and can't be shown as a *tested* non-significant
# example. It's replaced with 26->21: still non-significant (BH q ~= 0.09),
# and keeps node 21 as the target, which is the single most common whistle
# sub-category (1312/8566 occurrences) - its high emission frequency is
# exactly why the shuffle null also lands high probability there by chance,
# so even a fairly high observed probability into it isn't significant.
# e (16->24) is unchanged.

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(gtable)
  library(grid)
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
source(file.path(REPO_ROOT, "R", "statistical_utils.R"))
source(file.path(REPO_ROOT, "R", "plot_io.R"))

PLOT_SUBDIR <- "null_model_examples"
ensure_plots_dir(REPO_ROOT, PLOT_SUBDIR)

rds_path <- file.path(REPO_ROOT, OUTPUT_DIR, "markov_model.rds")
if (!file.exists(rds_path)) {
  stop("Run scripts/01_run_markov_model.R first. Missing: ", rds_path, call. = FALSE)
}
bh_q_path <- file.path(REPO_ROOT, OUTPUT_DIR, "bh_analysis", "bh_q_values.rds")
if (!file.exists(bh_q_path)) {
  stop("Run scripts/01_run_markov_model.R first. Missing: ", bh_q_path, call. = FALSE)
}

saved <- readRDS(rds_path)
MarkovModel <- saved$MarkovModel
bh_q_values <- readRDS(bh_q_path)

# Null count matrices are row counts, not row-normalized probabilities
# (MarkovWhistle() only pre-normalizes the observed matrix). Re-deriving
# P(from -> to) per iteration from just the needed row - rather than calling
# CalcTransitionMatrix() (a full 86x86 divide) on all
# MarkovModel$iterations_number matrices - avoids doing 86x that work for
# four cells. Ties CalcTransitionMatrix's own NaN -> 0 rule (row_total == 0
# yields 0, not NaN) so the two stay equivalent for these cells.
null_transition_probability <- function(null_count_matrices, from, to) {
  vapply(null_count_matrices, function(counts) {
    row_total <- sum(counts[from, ])
    if (row_total == 0) 0 else counts[from, to] / row_total
  }, numeric(1))
}

panels <- list(
  list(tag = "b", from = 8, to = 2, broken = FALSE),
  list(tag = "c", from = 40, to = 37, broken = TRUE),
  list(tag = "d", from = 26, to = 21, broken = FALSE),
  list(tag = "e", from = 16, to = 24, broken = TRUE)
)

panel_stats <- lapply(panels, function(p) {
  observed <- MarkovModel$transition_probabilities_matrix_all[p$from, p$to]
  null_probability <- null_transition_probability(
    MarkovModel$random_whistle_occurrences_matrices, p$from, p$to
  )
  c(
    p,
    list(
      observed = observed,
      count = MarkovModel$transition_counts_matrix[p$from, p$to],
      null_probability = null_probability,
      empirical_p = MarkovModel$p_value_matrix[p$from, p$to],
      bh_q = bh_q_values[p$from, p$to]
    )
  )
})

shared_theme <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  plot.title = element_text(size = 13, hjust = 0),
  plot.tag = element_text(size = 16, face = "bold"),
  plot.tag.position = c(0.02, 0.98),
  legend.background = element_blank(),
  legend.key = element_blank(),
  legend.text = element_text(size = 10)
)

# Panels without a dominant near-zero bin: one histogram, full y-range.
build_null_panel <- function(stats, show_legend = FALSE, show_tag = TRUE) {
  df <- data.frame(probability = stats$null_probability)
  ggplot(df, aes(x = probability)) +
    geom_histogram(aes(fill = "Null model"), bins = 40, colour = NA) +
    geom_vline(
      aes(xintercept = stats$observed, colour = "Observed value"),
      linetype = "dashed", linewidth = 0.8
    ) +
    scale_fill_manual(name = NULL, values = c("Null model" = "grey60")) +
    scale_colour_manual(name = NULL, values = c("Observed value" = "red")) +
    labs(
      title = sprintf("P(%d to %d)", stats$from, stats$to),
      x = "Probability", y = "# of events",
      tag = if (show_tag) stats$tag else NULL
    ) +
    shared_theme +
    theme(
      legend.position = if (show_legend) "inside" else "none",
      legend.position.inside = c(0.72, 0.85)
    )
}

# Break bounds for the two-panel "broken axis" trick used below: for these
# pairs the near-zero bin is 6-9x taller than the next-tallest bin (see the
# per-panel ratios this script prints), which would otherwise swamp the rest
# of the null distribution into a flat line at the bottom of a single
# shared y-scale. bottom_max covers all of the non-spike mass with headroom;
# top_min/top_max bracket the spike bin itself with headroom on both sides,
# so no bar is clipped by the break - only the visually empty gap between
# them is compressed out. Derived from the actual binned counts each run
# rather than pinned, so it stays correct if the pair selection changes.
compute_break_bounds <- function(bin_counts) {
  spike <- max(bin_counts)
  rest <- bin_counts[bin_counts < spike]
  rest_max <- if (length(rest)) max(rest) else 0
  list(bottom_max = rest_max * 1.15, top_min = spike * 0.85, top_max = spike * 1.08)
}

# Two short parallel diagonal segments marking the break - drawn once as a
# grob spanning the exact seam row of the combined gtable (see
# build_broken_null_panel() below), not positioned in either sub-panel's own
# data coordinates. An earlier version placed the marks via each sub-panel's
# own annotate("segment", ...) in data space, which needed
# coord_cartesian(clip = "off") to keep them from being cut off - but that
# also unclipped the histogram bars themselves (which extend far past each
# sub-panel's narrow ylim) and let them bleed into neighbouring panels.
# Anchoring to the gtable row instead sidesteps both problems and guarantees
# the marks sit exactly on the seam regardless of each panel's y-range.
seam_break_grob <- function() {
  segmentsGrob(
    x0 = unit(c(-0.025, 0.045), "npc"), x1 = unit(c(0.025, 0.095), "npc"),
    y0 = unit(c(-0.7, -0.7), "npc"), y1 = unit(c(0.7, 0.7), "npc"),
    gp = gpar(lwd = 1.6, lineend = "butt")
  )
}

# Panels with a dominant near-zero bin: a short "spike" sub-panel stacked on
# a taller "rest of the distribution" sub-panel. Built as two separate
# ggplots, then combined via gtable rbind (rather than gridExtra::arrangeGrob,
# used everywhere else in this script) because arrangeGrob lays each grob out
# in its own independently-sized cell - if the two panels' y-axis labels have
# different widths (e.g. "120000" on top vs "5000" on bottom), their plot
# areas end up different widths too, so the histogram bars don't line up
# vertically between the two sub-panels. gtable::rbind(..., size = "max")
# unions the column widths first, which keeps both panels' bars aligned to
# the same x scale.
build_broken_null_panel <- function(stats, show_tag = TRUE) {
  breaks <- seq(0, max(stats$null_probability) * 1.02, length.out = 41)
  x_range <- range(breaks)
  bin_counts <- hist(stats$null_probability, breaks = breaks, plot = FALSE)$counts
  bounds <- compute_break_bounds(bin_counts)
  message(sprintf(
    "  P(%d to %d) break: near-zero bin = %d, next-tallest bin = %d (%.1fx)",
    stats$from, stats$to, max(bin_counts), max(bin_counts[bin_counts < max(bin_counts)]),
    max(bin_counts) / max(bin_counts[bin_counts < max(bin_counts)])
  ))

  df <- data.frame(probability = stats$null_probability)
  base <- ggplot(df, aes(x = probability)) +
    geom_histogram(aes(fill = "Null model"), breaks = breaks, colour = NA) +
    scale_fill_manual(name = NULL, values = c("Null model" = "grey60")) +
    shared_theme +
    theme(legend.position = "none")

  top <- base +
    coord_cartesian(ylim = c(bounds$top_min, bounds$top_max), xlim = x_range) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 2)) +
    labs(
      title = sprintf("P(%d to %d)", stats$from, stats$to),
      tag = if (show_tag) stats$tag else NULL, y = NULL, x = NULL
    ) +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.title.x = element_blank(), axis.line.x = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 1, l = 5)
    )

  bottom <- base +
    geom_vline(
      aes(xintercept = stats$observed, colour = "Observed value"),
      linetype = "dashed", linewidth = 0.8
    ) +
    scale_colour_manual(name = NULL, values = c("Observed value" = "red")) +
    coord_cartesian(ylim = c(0, bounds$bottom_max), xlim = x_range) +
    labs(x = "Probability", y = "# of events") +
    theme(plot.margin = margin(t = 1, r = 5, b = 5, l = 5))

  top_grob <- ggplotGrob(top)
  bottom_grob <- ggplotGrob(bottom)
  common_widths <- unit.pmax(top_grob$widths, bottom_grob$widths)
  top_grob$widths <- common_widths
  bottom_grob$widths <- common_widths
  # A 1:3.3 null-unit split, tall enough for the top strip's title, tag and
  # 2-break y-axis to fit without crowding, short enough that it reads as a
  # compressed "just the spike" strip rather than a second full panel.
  top_grob$heights[top_grob$layout$t[top_grob$layout$name == "panel"]] <- unit(1, "null")
  bottom_grob$heights[bottom_grob$layout$t[bottom_grob$layout$name == "panel"]] <- unit(3.3, "null")

  seam_row <- nrow(top_grob)
  panel_col <- top_grob$layout$l[top_grob$layout$name == "panel"]
  combined <- rbind(top_grob, bottom_grob, size = "max")
  gtable_add_grob(
    combined, seam_break_grob(),
    t = seam_row, b = seam_row + 1, l = panel_col, clip = "off", name = "break_marks"
  )
}

build_panel <- function(stats, show_legend = FALSE, show_tag = TRUE) {
  if (isTRUE(stats$broken)) {
    build_broken_null_panel(stats, show_tag)
  } else {
    build_null_panel(stats, show_legend, show_tag)
  }
}

message("Broken-axis bounds (derived from each panel's own null histogram):")
combined_panels <- Map(
  build_panel, panel_stats, show_legend = c(TRUE, FALSE, FALSE, FALSE)
)
save_ggplot(
  gridExtra::arrangeGrob(grobs = combined_panels, ncol = 1),
  "SFigure_null_model_examples", REPO_ROOT,
  width = 3.6, height = 10, subdir = PLOT_SUBDIR
)

# Separate per-panel files: no panel tag (each is its own standalone figure,
# not a lettered piece of a multi-panel one), and only panel b carries the
# legend (matching the combined figure's convention of labelling it once).
standalone_names <- c(b = "null_model_b_8to2", c = "null_model_c_40to37",
                       d = "null_model_d_26to21", e = "null_model_e_16to24")
for (stats in panel_stats) {
  save_ggplot(
    build_panel(stats, show_legend = stats$tag == "b", show_tag = FALSE),
    standalone_names[[stats$tag]], REPO_ROOT,
    width = 3.6, height = 3, subdir = PLOT_SUBDIR
  )
}

stats_table <- do.call(rbind, lapply(panel_stats, function(s) {
  data.frame(
    panel = s$tag, from = s$from, to = s$to,
    observed_probability = s$observed, observed_count = s$count,
    iterations = MarkovModel$iterations_number,
    empirical_p = s$empirical_p, bh_q = s$bh_q,
    significant_at_0.05 = s$bh_q < GRAPH_ALPHA
  )
}))
stats_csv <- file.path(plots_dir(REPO_ROOT, PLOT_SUBDIR), "SFigure_null_model_examples_stats.csv")
write.csv(stats_table, stats_csv, row.names = FALSE)
message("Saved ", stats_csv)
message(
  "Updated legend values (BH-corrected q, ", MarkovModel$iterations_number,
  " shuffle iterations, direct randomization test):"
)
for (s in panel_stats) {
  message(sprintf(
    "  panel %s: P(%d to %d) = %.4f, BH q %s",
    s$tag, s$from, s$to, s$observed,
    if (s$bh_q < 0.001) "< 0.001" else sprintf("= %.3f", s$bh_q)
  ))
}

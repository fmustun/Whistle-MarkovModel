#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(igraph)
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

ensure_plots_dir(REPO_ROOT)
multiloop_dir <- file.path(REPO_ROOT, OUTPUT_DIR, "multiloop_analysis")
dir.create(multiloop_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading whistles from ", DATA_PATH)
whistles_list <- load_whistles(DATA_PATH)
message("Rows after filtering Unknown: ", nrow(whistles_list))

message(
  "Grouping whistles into multi-loop chains (IWI < ",
  MULTILOOP_IWI_THRESHOLD * 1000, "ms)..."
)
whistle_chains <- compute_multiloop_chains(
  whistles_list,
  threshold = MULTILOOP_IWI_THRESHOLD,
  list_names = LIST_NAMES
)
chains <- summarize_multiloop_chains(whistle_chains)

n_whistles <- nrow(whistle_chains)
n_chains <- nrow(chains)
n_singleton_whistles <- sum(chains$chain_length == 1)
n_multiloop_whistles <- n_whistles - n_singleton_whistles

# A chain of length 1 is a single whistle with no neighbor within the IWI
# threshold - not a "multi-loop" by definition. It's reported above for
# context, but excluded from every multi-loop measure/plot below so it
# doesn't dilute stats like mean chain length.
multiloop_chains <- chains[chains$chain_length >= 2, ]
n_multiloop_chains <- nrow(multiloop_chains)

message(
  n_chains, " chains from ", n_whistles, " whistles (",
  n_singleton_whistles, " singleton, ",
  n_multiloop_whistles, " whistles across ",
  n_multiloop_chains, " multi-loop chains)"
)

whistles_csv <- file.path(multiloop_dir, "multiloop_whistles.csv")
chains_csv <- file.path(multiloop_dir, "multiloop_chains.csv")
write.csv(whistle_chains, whistles_csv, row.names = FALSE)
write.csv(chains, chains_csv, row.names = FALSE)
message("Saved ", whistles_csv)
message("Saved ", chains_csv)

length_counts <- table(multiloop_chains$chain_length)

# Whether a multi-loop chain mixes whistles from more than one main category
# (e.g. "SW_Neo" vs "SW_Luna") or stays within a single one - main category
# here means the LIST_NAMES-level grouping (V(gra)$category elsewhere in
# this repo), not the finer whistle_type_chr sub-category.
whistle_chains_ml <- whistle_chains[
  whistle_chains$chain_id %in% multiloop_chains$chain_id,
]
chain_category_counts <- aggregate(
  category ~ chain_id, data = whistle_chains_ml,
  FUN = function(x) length(unique(x))
)
names(chain_category_counts)[2] <- "n_distinct_categories"
chain_category_counts$composition <- ifelse(
  chain_category_counts$n_distinct_categories == 1L,
  "same category", "different category"
)
composition_counts <- as.data.frame(table(chain_category_counts$composition))
names(composition_counts) <- c("composition", "n_chains")
composition_counts$pct <- 100 * composition_counts$n_chains / sum(composition_counts$n_chains)
n_same_category_chains <- composition_counts$n_chains[
  composition_counts$composition == "same category"
]
n_diff_category_chains <- composition_counts$n_chains[
  composition_counts$composition == "different category"
]

message(
  n_same_category_chains, " / ", n_multiloop_chains,
  " multi-loop chains (",
  round(100 * n_same_category_chains / n_multiloop_chains, 1),
  "%) stay within a single main category"
)

# How much of the Markov model's transition pairs (script 01, TIME_WINDOW)
# connect two whistles that are themselves part of the same tight-IWI
# multi-loop chain (MULTILOOP_IWI_THRESHOLD)? These are two independent
# "next whistle" definitions - this cross-references them at the level of
# individual recorded transition instances (see
# compute_markov_transition_multiloop_overlap(), R/multiloop_utils.R).
message("Cross-referencing Markov transition pairs against multi-loop chains...")
transition_overlap <- compute_markov_transition_multiloop_overlap(
  whistles_list,
  time_window = TIME_WINDOW,
  threshold = MULTILOOP_IWI_THRESHOLD,
  list_names = LIST_NAMES
)

n_transitions <- nrow(transition_overlap)
n_transitions_within_chain <- sum(transition_overlap$within_multiloop_chain)
pct_within_chain <- 100 * n_transitions_within_chain / n_transitions

is_self_transition <- transition_overlap$from_whistle_type ==
  transition_overlap$to_whistle_type
n_self <- sum(is_self_transition)
n_self_within_chain <- sum(transition_overlap$within_multiloop_chain[is_self_transition])
n_diff <- n_transitions - n_self
n_diff_within_chain <- n_transitions_within_chain - n_self_within_chain

message(
  n_transitions_within_chain, " / ", n_transitions,
  " Markov transition instances (", round(pct_within_chain, 1),
  "%) connect two whistles in the same multi-loop chain"
)

overlap_csv <- file.path(
  multiloop_dir, "markov_transition_multiloop_overlap.csv"
)
write.csv(transition_overlap, overlap_csv, row.names = FALSE)
message("Saved ", overlap_csv)

# -- write the combined text summary --------------------------------------

summary_lines <- c(
  "Multi-loop whistle population summary",
  paste("IWI threshold:", MULTILOOP_IWI_THRESHOLD * 1000, "ms",
        "(gap = start_time[i+1] - end_time[i])"),
  paste("total whistles:", n_whistles),
  paste("total chains (including singletons):", n_chains),
  paste0(
    "singleton whistles (chain length 1, not a multi-loop): ",
    n_singleton_whistles,
    " (", round(100 * n_singleton_whistles / n_whistles, 1), "%)"
  ),
  "",
  "-- measures below cover multi-loop chains only (chain length >= 2) --",
  paste0(
    "whistles in multi-loop chains: ", n_multiloop_whistles,
    " (", round(100 * n_multiloop_whistles / n_whistles, 1), "%)"
  ),
  paste("multi-loop chains:", n_multiloop_chains),
  paste("mean chain length:", round(mean(multiloop_chains$chain_length), 2)),
  paste("median chain length:", median(multiloop_chains$chain_length)),
  paste("max chain length:", max(multiloop_chains$chain_length)),
  "chain-length counts:",
  paste(" ", names(length_counts), ":", as.integer(length_counts)),
  "",
  "main-category composition of multi-loop chains:",
  paste0(
    "  same category: ", n_same_category_chains, " / ", n_multiloop_chains,
    " (", round(100 * n_same_category_chains / n_multiloop_chains, 1), "%)"
  ),
  paste0(
    "  different category: ", n_diff_category_chains, " / ", n_multiloop_chains,
    " (", round(100 * n_diff_category_chains / n_multiloop_chains, 1), "%)"
  ),
  "",
  "Markov model transition / multi-loop chain overlap",
  paste(
    "TIME_WINDOW:", TIME_WINDOW[1], "-", TIME_WINDOW[2], "s",
    "(script 01's Markov transition pair selection)"
  ),
  paste("total Markov transition instances:", n_transitions),
  paste0(
    "within the same multi-loop chain: ", n_transitions_within_chain,
    " (", round(pct_within_chain, 1), "%)"
  ),
  paste0(
    "  same whistle type (self-transitions): ", n_self_within_chain,
    " / ", n_self,
    if (n_self > 0) {
      paste0(" (", round(100 * n_self_within_chain / n_self, 1), "%)")
    } else {
      ""
    }
  ),
  paste0(
    "  different whistle type: ", n_diff_within_chain, " / ", n_diff,
    if (n_diff > 0) {
      paste0(" (", round(100 * n_diff_within_chain / n_diff, 1), "%)")
    } else {
      ""
    }
  )
)
writeLines(summary_lines, file.path(multiloop_dir, "multiloop_summary.txt"))
message("Saved ", file.path(multiloop_dir, "multiloop_summary.txt"))

# -- combined figure: three panels -----------------------------------------

panel_length_by_category <- ggplot(
  multiloop_chains, aes(x = factor(chain_length), fill = leading_category)
) +
  geom_bar(color = "black", width = 0.75) +
  scale_fill_manual(
    breaks = LIST_NAMES, values = LIST_COLORS[seq_along(LIST_NAMES)]
  ) +
  theme_minimal() +
  labs(
    title = paste0(
      "Multi-loop chain length distribution by leading category (IWI < ",
      MULTILOOP_IWI_THRESHOLD * 1000, "ms; singletons excluded)"
    ),
    x = "chain length (whistles)",
    y = "number of chains",
    fill = "leading category"
  )

panel_category_composition <- ggplot(
  composition_counts, aes(x = "", y = n_chains, fill = composition)
) +
  geom_col(color = "black", width = 1) +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = paste0(n_chains, "\n(", round(pct, 1), "%)")),
    position = position_stack(vjust = 0.5)
  ) +
  scale_fill_manual(
    breaks = c("same category", "different category"),
    values = c("#61D04F", "#A9A9A9")
  ) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    title = "Multi-loop chains: same vs. different main category",
    fill = "composition"
  )

overlap_rates <- data.frame(
  group = factor(
    c("all transitions", "self-transitions", "different-type transitions"),
    levels = c("all transitions", "self-transitions", "different-type transitions")
  ),
  pct_within_chain = c(
    pct_within_chain,
    100 * n_self_within_chain / n_self,
    100 * n_diff_within_chain / n_diff
  ),
  n = c(n_transitions, n_self, n_diff),
  n_within = c(n_transitions_within_chain, n_self_within_chain, n_diff_within_chain)
)

panel_markov_overlap <- ggplot(
  overlap_rates, aes(x = group, y = pct_within_chain)
) +
  geom_col(fill = "#2297E6", color = "black", width = 0.6) +
  geom_text(
    aes(label = paste0(round(pct_within_chain, 1), "%\n(", n_within, "/", n, ")")),
    vjust = -0.15, size = 3.2
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(
    title = "Markov transition pairs within a multi-loop chain",
    x = NULL,
    y = "% of transition instances"
  ) +
  ylim(0, 100)

layout_matrix <- rbind(c(1, 1), c(2, 3))
combined_plot <- gridExtra::arrangeGrob(
  panel_length_by_category, panel_category_composition, panel_markov_overlap,
  layout_matrix = layout_matrix
)
save_ggplot(combined_plot, "multiloop_combined_summary", REPO_ROOT, width = 11, height = 12)

# Remove the earlier per-panel plots this combined figure replaces.
for (name in c("multiloop_length_distribution", "multiloop_length_distribution_by_category")) {
  for (ext in c("pdf", "png")) {
    stale_path <- plot_path(name, REPO_ROOT, ext = ext)
    if (file.exists(stale_path)) file.remove(stale_path)
  }
}

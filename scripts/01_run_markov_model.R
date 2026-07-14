#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(igraph)
  library(scales)
  library(ccber)
  library(foreach)
  library(doSNOW)
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
source(file.path(REPO_ROOT, "R", "markov_whistle.R"))
source(file.path(REPO_ROOT, "R", "statistical_utils.R"))
source(file.path(REPO_ROOT, "R", "graph_utils.R"))
source(file.path(REPO_ROOT, "R", "plot_io.R"))

dir.create(file.path(REPO_ROOT, OUTPUT_DIR), showWarnings = FALSE, recursive = TRUE)
ensure_plots_dir(REPO_ROOT)

message("Loading whistles from ", DATA_PATH)
whistles_list <- load_whistles(DATA_PATH)
message("Rows after filtering Unknown: ", nrow(whistles_list))

set.seed(SEED)
message("Running MarkovWhistle (iterations = ", ITERATIONS, ")...")
MarkovModel <- MarkovWhistle(
  whistles_list,
  null_model = NULL_MODEL,
  min_shift = MIN_SHIFT,
  max_shift = MAX_SHIFT,
  time_window = TIME_WINDOW,
  iterations = ITERATIONS,
  threshold = THRESHOLD,
  list_names = LIST_NAMES
)

register_fcircle_shape()
# Never-observed transitions can never pass select_transition_matrix's own
# observed_probability > 0 requirement, so they are excluded from the
# multiple-testing family rather than diluting it (see
# adjust_transition_p_values' family_mask in R/statistical_utils.R).
family_mask <- MarkovModel$transition_probabilities_matrix_all > 0
adjusted_p_value <- adjust_transition_p_values(
  MarkovModel$p_value_matrix,
  correction = MULTIPLE_TESTING_CORRECTION,
  family_mask = family_mask
)
selection <- select_transition_matrix(
  MarkovModel$transition_probabilities_matrix_all,
  adjusted_p_value,
  alpha = GRAPH_ALPHA
)
selected_graphs <- compute_graph_from_selection(
  transition_matrix = selection$matrix,
  selection_values = adjusted_p_value,
  Markov_Model = MarkovModel,
  Whistles_List = whistles_list,
  list_names = LIST_NAMES,
  list_colors = LIST_COLORS,
  edge_attribute_name = if (MULTIPLE_TESTING_CORRECTION == "BH") {
    "bh_q"
  } else {
    "empirical_p"
  },
  raw_p_values = MarkovModel$p_value_matrix
)
gra1 <- selected_graphs$no_loops

counts <- transition_selection_counts(selection$mask)
# No pinned regression check here: the previous pin (150,000 iterations;
# retained_nodes=66 etc.) was tied to a hypothesis family that included
# never-observed transitions and to an iteration count that was tuned to
# reproduce a target result rather than derived from a resolution
# requirement (see recommended_iterations() in R/statistical_utils.R). Once
# you've run the derived iteration count once and trust the result, use
# assert_transition_selection_counts(counts, c(...)) to pin a new baseline.
message(paste(names(counts), counts, sep = "=", collapse = "; "))

message(
  "Plotting significant-transition network (",
  MULTIPLE_TESTING_CORRECTION, " adjusted p < ", GRAPH_ALPHA, ")"
)
with_pdf_plot(
  plot_path("markov_significant_network", REPO_ROOT),
  width = 8,
  height = 8,
  plot_markov_graph(gra1, list_names = LIST_NAMES, seed = PLOT_SEED)
)

out_path <- file.path(REPO_ROOT, OUTPUT_DIR, "markov_model.rds")
saveRDS(
  list(
    MarkovModel = MarkovModel,
    whistles_list = whistles_list,
    gra1 = gra1,
    time_window = TIME_WINDOW,
    graph_p_value = GRAPH_P_VALUE,
    config = list(
      iterations = ITERATIONS,
      null_model = NULL_MODEL,
      time_window = TIME_WINDOW,
      graph_p_value = GRAPH_P_VALUE,
      multiple_testing_correction = MULTIPLE_TESTING_CORRECTION,
      graph_alpha = GRAPH_ALPHA,
      list_names = LIST_NAMES
    )
  ),
  out_path
)
message("Saved ", out_path)

bh_dir <- file.path(REPO_ROOT, OUTPUT_DIR, "bh_analysis")
dir.create(bh_dir, recursive = TRUE, showWarnings = FALSE)

label_map <- tapply(
  whistles_list$whistle_type_chr,
  whistles_list$whistle_type,
  function(value) as.character(value[1L])
)
coordinates <- compute_markov_layout(gra1, PLOT_SEED)

coordinate_table <- data.frame(
  node_id = as.integer(rownames(coordinates)),
  label = unname(label_map[rownames(coordinates)]),
  category = V(gra1)$category,
  x = coordinates[, "x"],
  y = coordinates[, "y"],
  stringsAsFactors = FALSE
)

retained_index <- which(selection$mask, arr.ind = TRUE)
edge_table <- data.frame(
  from_node_id = retained_index[, 1L],
  to_node_id = retained_index[, 2L],
  from_label = unname(label_map[as.character(retained_index[, 1L])]),
  to_label = unname(label_map[as.character(retained_index[, 2L])]),
  transition_probability = MarkovModel$transition_probabilities_matrix_all[retained_index],
  empirical_p = MarkovModel$p_value_matrix[retained_index],
  adjusted_p = adjusted_p_value[retained_index],
  correction = MULTIPLE_TESTING_CORRECTION,
  is_auto_loop = retained_index[, 1L] == retained_index[, 2L],
  stringsAsFactors = FALSE
)
if (MULTIPLE_TESTING_CORRECTION == "BH") {
  names(edge_table)[names(edge_table) == "adjusted_p"] <- "bh_q"
}

write.csv(MarkovModel$p_value_matrix,
          file.path(bh_dir, "empirical_p_values.csv"))
write.csv(adjusted_p_value, file.path(bh_dir,
          if (MULTIPLE_TESTING_CORRECTION == "BH")
            "bh_q_values.csv" else "uncorrected_p_values.csv"))
saveRDS(MarkovModel$p_value_matrix,
        file.path(bh_dir, "empirical_p_values.rds"))
saveRDS(adjusted_p_value, file.path(bh_dir,
        if (MULTIPLE_TESTING_CORRECTION == "BH")
          "bh_q_values.rds" else "uncorrected_p_values.rds"))
saveRDS(MarkovModel$empirical_exceedance_count_matrix,
        file.path(bh_dir, "empirical_exceedance_counts.rds"))
write.csv(coordinate_table, file.path(bh_dir, "node_layout_coordinates.csv"),
          row.names = FALSE)
write.csv(edge_table, file.path(bh_dir, "bh_retained_edges.csv"),
          row.names = FALSE)
write.csv(data.frame(metric = names(counts), value = as.integer(counts)),
          file.path(bh_dir, "bh_validation_counts.csv"), row.names = FALSE)
saveRDS(
  list(
    iterations = ITERATIONS,
    hypothesis_count = sum(family_mask),
    empirical_p_formula = "(1 + count(null >= observed)) / (B + 1)",
    zero_observed_rule = "p = 1",
    correction = MULTIPLE_TESTING_CORRECTION,
    alpha = GRAPH_ALPHA,
    empirical_p_values = MarkovModel$p_value_matrix,
    adjusted_p_values = adjusted_p_value,
    significant_transition_matrix = selection$matrix,
    significant_mask = selection$mask,
    selection_counts = counts,
    graph_with_loops = selected_graphs$with_loops,
    graph = gra1,
    significant_auto_loop_ids = selected_graphs$significant_loop_ids,
    layout_coordinates = coordinates
  ),
  file.path(bh_dir, "bh_corrected_network.rds")
)

writeLines(
  c(
    "BH-corrected significant network",
    paste("iterations:", ITERATIONS),
    paste(
      "hypothesis family:", sum(family_mask),
      "transitions observed at least once (loops included);",
      "never-observed transitions are excluded (see family_mask in",
      "adjust_transition_p_values, R/statistical_utils.R)"
    ),
    "empirical p: (1 + count(null >= observed)) / (B + 1); zero-observed p = 1",
    paste("correction:", MULTIPLE_TESTING_CORRECTION),
    paste("alpha:", GRAPH_ALPHA),
    paste(names(counts), counts, sep = ": "),
    "loops are represented by black node borders and omitted as drawn edges"
  ),
  file.path(bh_dir, "bh_analysis_summary.txt")
)
message("BH analysis outputs written to ", bh_dir)

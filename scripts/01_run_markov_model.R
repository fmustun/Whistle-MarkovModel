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
adjusted_p_value <- adjust_transition_p_values(
  MarkovModel$p_value_matrix,
  correction = MULTIPLE_TESTING_CORRECTION
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

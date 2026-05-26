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
source(file.path(REPO_ROOT, "R", "graph_utils.R"))

dir.create(file.path(REPO_ROOT, OUTPUT_DIR), showWarnings = FALSE, recursive = TRUE)

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
gra1 <- compute_graph_p_value_significant(
  GRAPH_P_VALUE,
  Markov_Model = MarkovModel,
  Whistles_List = whistles_list,
  sub_network = NULL,
  sub_division = TRUE,
  list_names = LIST_NAMES
)

message("Plotting significant-transition network (p < ", GRAPH_P_VALUE, ")")
plot_markov_graph(gra1, list_names = LIST_NAMES, seed = PLOT_SEED)

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
      list_names = LIST_NAMES
    )
  ),
  out_path
)
message("Saved ", out_path)

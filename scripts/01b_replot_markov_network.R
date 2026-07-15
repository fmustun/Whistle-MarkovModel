#!/usr/bin/env Rscript

# Regenerates the significant-transition network plot from the already-saved
# outputs/markov_model.rds, without rerunning MarkovWhistle's null-model
# simulation (see scripts/01_run_markov_model.R for the full pipeline).

suppressPackageStartupMessages({
  library(igraph)
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
gra1 <- saved$gra1

register_fcircle_shape()

message(
  "Plotting significant-transition network (",
  saved$config$multiple_testing_correction, " adjusted p < ",
  saved$config$graph_alpha, ")"
)
with_pdf_plot(
  plot_path("markov_significant_network", REPO_ROOT),
  width = 8,
  height = 8,
  plot_markov_graph(
    gra1, list_names = saved$config$list_names, seed = PLOT_SEED,
    show_legend = FALSE
  )
)

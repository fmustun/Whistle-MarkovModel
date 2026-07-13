#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(igraph)
  library(ccber)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
REPO_ROOT <- if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."),
                mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
setwd(REPO_ROOT)

source(file.path(REPO_ROOT, "R", "config.R"))
source(file.path(REPO_ROOT, "R", "statistical_utils.R"))
source(file.path(REPO_ROOT, "R", "graph_utils.R"))

parse_options <- function(arguments) {
  defaults <- list(
    input = file.path(REPO_ROOT, OUTPUT_DIR, "markov_model.rds"),
    output_dir = file.path(REPO_ROOT, OUTPUT_DIR, "bh_analysis"),
    correction = MULTIPLE_TESTING_CORRECTION,
    alpha = GRAPH_ALPHA,
    precomputed_p = ""
  )
  for (argument in arguments) {
    match <- regexec("^--([^=]+)=(.*)$", argument)
    parts <- regmatches(argument, match)[[1L]]
    if (!length(parts)) stop("Arguments must use --name=value: ", argument)
    name <- gsub("-", "_", parts[2L], fixed = TRUE)
    if (!name %in% names(defaults)) stop("Unknown option: --", parts[2L])
    defaults[[name]] <- parts[3L]
  }
  defaults$alpha <- as.numeric(defaults$alpha)
  defaults
}

options <- parse_options(commandArgs(trailingOnly = TRUE))
options$correction <- match.arg(options$correction, c("BH", "none"))
if (!is.finite(options$alpha) || options$alpha <= 0 || options$alpha >= 1) {
  stop("alpha must be between zero and one", call. = FALSE)
}
if (!file.exists(options$input)) stop("Missing input: ", options$input)

if (dir.exists(options$output_dir) &&
    length(list.files(options$output_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Refusing to overwrite non-empty output directory: ",
       options$output_dir, call. = FALSE)
}
dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading saved null model: ", options$input)
saved <- readRDS(options$input)
validate_saved_markov_model(saved, expected_states = 86L)
model <- saved$MarkovModel
observed <- name_transition_matrix(model$transition_probabilities_matrix_all)
state_ids <- rownames(observed)

if (nzchar(options$precomputed_p)) {
  message("Loading precomputed empirical p-values: ", options$precomputed_p)
  empirical_p <- name_transition_matrix(readRDS(options$precomputed_p))
  if (!identical(dim(empirical_p), dim(observed)) ||
      any(!is.finite(empirical_p)) || any(empirical_p <= 0) ||
      any(empirical_p > 1) || any(empirical_p[observed == 0] != 1)) {
    stop("Precomputed empirical p-value matrix failed validation", call. = FALSE)
  }
  exceedance_count <- NULL
} else {
  result <- empirical_enrichment_p_values(
    observed,
    model$random_whistle_occurrences_matrices,
    ccber::CalcTransitionMatrix
  )
  empirical_p <- result$p_value
  exceedance_count <- result$exceedance_count
}

adjusted_p <- adjust_transition_p_values(
  empirical_p, correction = options$correction
)
selection <- select_transition_matrix(observed, adjusted_p, options$alpha)
counts <- transition_selection_counts(selection$mask)

if (identical(dim(observed), c(86L, 86L)) &&
    identical(as.integer(model$iterations_number), 150000L) &&
    options$correction == "BH" && identical(options$alpha, 0.05)) {
  assert_transition_selection_counts(
    counts,
    c(
      retained_nodes = 66L,
      significant_transitions_including_loops = 128L,
      significant_auto_loops = 31L,
      loop_free_edges = 97L
    )
  )
}
message(paste(names(counts), counts, sep = "=", collapse = "; "))

graphs <- compute_graph_from_selection(
  selection$matrix,
  adjusted_p,
  Markov_Model = model,
  Whistles_List = saved$whistles_list,
  list_names = saved$config$list_names,
  list_colors = LIST_COLORS,
  edge_attribute_name = if (options$correction == "BH") "bh_q" else "empirical_p",
  raw_p_values = empirical_p
)
bh_graph <- graphs$no_loops
original_graph <- saved$gra1
register_fcircle_shape()

original_coordinates <- compute_markov_layout(original_graph, PLOT_SEED)
bh_coordinates <- layout_for_graph(original_coordinates, bh_graph)
weight_reference_range <- range(E(original_graph)$weight)

label_map <- tapply(
  saved$whistles_list$whistle_type_chr,
  saved$whistles_list$whistle_type,
  function(value) as.character(value[1L])
)
coordinate_table <- data.frame(
  node_id = as.integer(rownames(original_coordinates)),
  label = unname(label_map[rownames(original_coordinates)]),
  category = V(original_graph)$category,
  x = original_coordinates[, "x"],
  y = original_coordinates[, "y"],
  bh_retained = rownames(original_coordinates) %in% V(bh_graph)$name,
  stringsAsFactors = FALSE
)

retained_index <- which(selection$mask, arr.ind = TRUE)
edge_table <- data.frame(
  from_node_id = retained_index[, 1L],
  to_node_id = retained_index[, 2L],
  from_label = unname(label_map[as.character(retained_index[, 1L])]),
  to_label = unname(label_map[as.character(retained_index[, 2L])]),
  transition_probability = observed[retained_index],
  empirical_p = empirical_p[retained_index],
  adjusted_p = adjusted_p[retained_index],
  correction = options$correction,
  is_auto_loop = retained_index[, 1L] == retained_index[, 2L],
  stringsAsFactors = FALSE
)
if (options$correction == "BH") {
  names(edge_table)[names(edge_table) == "adjusted_p"] <- "bh_q"
}

write.csv(empirical_p, file.path(options$output_dir, "empirical_p_values.csv"))
write.csv(adjusted_p, file.path(options$output_dir,
                                if (options$correction == "BH")
                                  "bh_q_values.csv" else "uncorrected_p_values.csv"))
saveRDS(empirical_p, file.path(options$output_dir, "empirical_p_values.rds"))
saveRDS(adjusted_p, file.path(options$output_dir,
                              if (options$correction == "BH")
                                "bh_q_values.rds" else "uncorrected_p_values.rds"))
if (!is.null(exceedance_count)) {
  saveRDS(exceedance_count,
          file.path(options$output_dir, "empirical_exceedance_counts.rds"))
}
write.csv(coordinate_table,
          file.path(options$output_dir, "figure4_original_layout_coordinates.csv"),
          row.names = FALSE)
write.csv(edge_table, file.path(options$output_dir, "bh_retained_edges.csv"),
          row.names = FALSE)
write.csv(data.frame(metric = names(counts), value = as.integer(counts)),
          file.path(options$output_dir, "bh_validation_counts.csv"),
          row.names = FALSE)
saveRDS(
  list(
    source_model = normalizePath(options$input),
    source_iterations = model$iterations_number,
    hypothesis_count = length(empirical_p),
    empirical_p_formula = "(1 + count(null >= observed)) / (B + 1)",
    zero_observed_rule = "p = 1",
    correction = options$correction,
    alpha = options$alpha,
    empirical_p_values = empirical_p,
    adjusted_p_values = adjusted_p,
    significant_transition_matrix = selection$matrix,
    significant_mask = selection$mask,
    selection_counts = counts,
    graph_with_loops = graphs$with_loops,
    graph = graphs$no_loops,
    significant_auto_loop_ids = graphs$significant_loop_ids,
    original_layout_coordinates = original_coordinates,
    graph_layout_coordinates = bh_coordinates
  ),
  file.path(options$output_dir, "bh_corrected_network.rds")
)

draw_network <- function(graph, title = NULL, legend = TRUE) {
  par(mar = c(5.2, 1.2, if (is.null(title)) 1.2 else 2.5, 1.2))
  plot_markov_graph(
    graph,
    list_names = saved$config$list_names,
    list_colors = LIST_COLORS,
    coordinates = original_coordinates,
    weight_reference_range = weight_reference_range,
    show_legend = legend,
    main = title
  )
}

pdf_path <- file.path(options$output_dir, "figure4_bh_original_style.pdf")
pdf(pdf_path, width = 8, height = 8, useDingbats = FALSE)
draw_network(bh_graph)
dev.off()

png_path <- file.path(options$output_dir, "figure4_bh_original_style.png")
png(png_path, width = 3200, height = 3200, res = 400,
    type = "cairo", bg = "white")
draw_network(bh_graph)
dev.off()

comparison_pdf <- file.path(options$output_dir, "figure4_original_vs_bh.pdf")
pdf(comparison_pdf, width = 16, height = 8, useDingbats = FALSE)
par(mfrow = c(1, 2))
draw_network(original_graph, "Original: raw p < 0.05", legend = FALSE)
draw_network(bh_graph, "BH q < 0.05", legend = TRUE)
dev.off()

comparison_png <- file.path(options$output_dir, "figure4_original_vs_bh.png")
png(comparison_png, width = 6400, height = 3200, res = 400,
    type = "cairo", bg = "white")
par(mfrow = c(1, 2))
draw_network(original_graph, "Original: raw p < 0.05", legend = FALSE)
draw_network(bh_graph, "BH q < 0.05", legend = TRUE)
dev.off()

writeLines(
  c(
    "BH-corrected Figure 4 network",
    paste("input:", normalizePath(options$input)),
    paste("iterations:", model$iterations_number),
    "hypothesis family: all 86 x 86 ordered transitions, including loops and zero-observed pairs",
    "empirical p: (1 + count(null >= observed)) / (B + 1); zero-observed p = 1",
    paste("correction:", options$correction),
    paste("alpha:", options$alpha),
    paste(names(counts), counts, sep = ": "),
    "layout: seed-22 original 85-node FR coordinates; BH uses the retained rows",
    paste("absolute edge-weight reference range:",
          paste(signif(weight_reference_range, 8), collapse = " to ")),
    "loops are represented by black node borders and omitted as drawn edges"
  ),
  file.path(options$output_dir, "bh_analysis_summary.txt")
)
message("BH analysis outputs written to ", options$output_dir)

validate_saved_markov_model <- function(saved, expected_states = 86L,
                                        expected_iterations = NULL) {
  required_saved <- c("MarkovModel", "whistles_list", "config")
  if (!is.list(saved) || !all(required_saved %in% names(saved))) {
    stop("Saved model lacks: ",
         paste(setdiff(required_saved, names(saved)), collapse = ", "),
         call. = FALSE)
  }

  model <- saved$MarkovModel
  required_model <- c(
    "transition_probabilities_matrix_all", "transition_counts_matrix",
    "random_whistle_occurrences_matrices", "iterations_number",
    "whistle_occurrence"
  )
  if (!all(required_model %in% names(model))) {
    stop("MarkovModel lacks: ",
         paste(setdiff(required_model, names(model)), collapse = ", "),
         call. = FALSE)
  }

  expected_dim <- c(as.integer(expected_states), as.integer(expected_states))
  matrices <- list(
    observed_probability = model$transition_probabilities_matrix_all,
    observed_counts = model$transition_counts_matrix
  )
  for (matrix_name in names(matrices)) {
    value <- matrices[[matrix_name]]
    if (!is.matrix(value) || !identical(dim(value), expected_dim) ||
        !is.numeric(value) || any(!is.finite(value)) || any(value < 0)) {
      stop(matrix_name, " is not a finite non-negative ", expected_states,
           " x ", expected_states, " numeric matrix", call. = FALSE)
    }
  }

  null_counts <- model$random_whistle_occurrences_matrices
  if (!is.null(expected_iterations) &&
      !identical(as.integer(model$iterations_number),
                 as.integer(expected_iterations))) {
    stop("Expected ", expected_iterations, " iterations; found ",
         model$iterations_number, call. = FALSE)
  }
  if (length(null_counts) != as.integer(model$iterations_number)) {
    stop("Saved null-matrix count differs from iterations_number", call. = FALSE)
  }
  null_valid <- vapply(
    null_counts,
    function(value) is.matrix(value) && identical(dim(value), expected_dim) &&
      is.numeric(value) && all(is.finite(value)) && all(value >= 0),
    logical(1)
  )
  if (!all(null_valid)) {
    stop("Invalid saved null count matrix at iteration ", which(!null_valid)[1L],
         call. = FALSE)
  }
  invisible(saved)
}

name_transition_matrix <- function(value) {
  if (!is.matrix(value) || nrow(value) != ncol(value)) {
    stop("Transition value must be a square matrix", call. = FALSE)
  }
  state_ids <- as.character(seq_len(nrow(value)))
  dimnames(value) <- list(state_ids, state_ids)
  value
}

empirical_enrichment_p_values <- function(
    observed_probability,
    null_count_matrices,
    transition_probability_function,
    progress_every = 5000L
) {
  observed_probability <- name_transition_matrix(observed_probability)
  iterations <- length(null_count_matrices)
  if (iterations < 1L) stop("At least one null matrix is required", call. = FALSE)

  # Hypothesis family: every ordered state pair, including auto-loops and
  # zero-observed transitions. Enrichment p = (1 + count(null >= observed)) /
  # (B + 1), with ties included. Zero-observed transitions remain in the family
  # but are assigned p = 1 after the exceedance pass.
  exceedances <- matrix(
    0L, nrow = nrow(observed_probability), ncol = ncol(observed_probability),
    dimnames = dimnames(observed_probability)
  )
  for (iteration in seq_len(iterations)) {
    null_probability <- transition_probability_function(
      null_count_matrices[[iteration]]
    )
    if (!identical(dim(null_probability), dim(observed_probability))) {
      stop("Null transition matrix has unexpected dimensions at iteration ",
           iteration, call. = FALSE)
    }
    exceedances <- exceedances + (null_probability >= observed_probability)
    if (progress_every > 0L &&
        (iteration %% progress_every == 0L || iteration == iterations)) {
      message("Empirical p-value pass: ", iteration, " / ", iterations)
    }
  }

  empirical_p <- (1 + exceedances) / (iterations + 1)
  empirical_p[observed_probability == 0] <- 1
  empirical_p <- name_transition_matrix(empirical_p)
  if (any(!is.finite(empirical_p)) || any(empirical_p <= 0) ||
      any(empirical_p > 1)) {
    stop("Invalid empirical p-values produced", call. = FALSE)
  }
  list(p_value = empirical_p, exceedance_count = exceedances)
}

adjust_transition_p_values <- function(empirical_p, correction = c("BH", "none")) {
  correction <- match.arg(correction)
  empirical_p <- name_transition_matrix(empirical_p)
  adjusted <- if (correction == "none") {
    empirical_p
  } else {
    matrix(
      stats::p.adjust(as.vector(empirical_p), method = correction),
      nrow = nrow(empirical_p), ncol = ncol(empirical_p),
      dimnames = dimnames(empirical_p)
    )
  }
  name_transition_matrix(adjusted)
}

select_transition_matrix <- function(observed_probability, adjusted_p,
                                     alpha = 0.05) {
  observed_probability <- name_transition_matrix(observed_probability)
  adjusted_p <- name_transition_matrix(adjusted_p)
  if (!identical(dim(observed_probability), dim(adjusted_p))) {
    stop("Observed and adjusted-p matrices have different dimensions",
         call. = FALSE)
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be one finite number strictly between zero and one",
         call. = FALSE)
  }
  selected <- observed_probability > 0 & adjusted_p < alpha
  transition_matrix <- observed_probability
  transition_matrix[!selected] <- 0
  list(matrix = transition_matrix, mask = selected)
}

transition_selection_counts <- function(mask) {
  active <- which(rowSums(mask) + colSums(mask) > 0)
  loops <- sum(diag(mask))
  c(
    retained_nodes = length(active),
    significant_transitions_including_loops = sum(mask),
    significant_auto_loops = loops,
    loop_free_edges = sum(mask) - loops
  )
}

assert_transition_selection_counts <- function(counts, expected) {
  shared <- intersect(names(expected), names(counts))
  if (!identical(as.integer(counts[shared]), as.integer(expected[shared]))) {
    stop(
      "Transition-selection validation failed. Observed: ",
      paste(names(counts), counts, sep = "=", collapse = ", "),
      "; expected: ",
      paste(names(expected), expected, sep = "=", collapse = ", "),
      call. = FALSE
    )
  }
  invisible(counts)
}

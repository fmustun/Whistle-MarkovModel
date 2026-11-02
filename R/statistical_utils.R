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

# family_mask restricts the multiple-testing family to the transitions
# actually eligible for selection (observed_probability > 0). Never-observed
# transitions can never pass select_transition_matrix's own positivity
# requirement, so including them in the correction inflates the family size
# for no statistical benefit (independent filtering; Bourgon et al. 2010).
# Entries outside family_mask are set to 1 (never significant), matching the
# zero_observed_rule already used for raw empirical p-values. family_mask
# defaults to NULL, which corrects across every entry (prior behavior).
adjust_transition_p_values <- function(empirical_p, correction = c("BH", "none"),
                                       family_mask = NULL) {
  correction <- match.arg(correction)
  empirical_p <- name_transition_matrix(empirical_p)
  if (!is.null(family_mask)) {
    if (!identical(dim(family_mask), dim(empirical_p))) {
      stop("family_mask must have the same dimensions as empirical_p",
           call. = FALSE)
    }
    adjusted <- matrix(
      1, nrow = nrow(empirical_p), ncol = ncol(empirical_p),
      dimnames = dimnames(empirical_p)
    )
    adjusted[family_mask] <- if (correction == "none") {
      empirical_p[family_mask]
    } else {
      stats::p.adjust(empirical_p[family_mask], method = correction)
    }
    return(name_transition_matrix(adjusted))
  }
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

# Minimum permutation count for a multiple-testing family of size m at
# significance level alpha to have any resolving power: the empirical
# p-value floor is 1/(B+1), so the strictest BH threshold (rank 1, alpha/m)
# requires B >= m/alpha - 1. margin scales past that bare floor so p-value
# estimates near the decision boundary aren't dominated by Monte Carlo noise
# (a hypothesis right at the floor has an expected exceedance count of only
# ~1 at the bare minimum; margin = 5 raises that to ~5, and so on).
recommended_iterations <- function(m, alpha, margin = 5) {
  if (length(m) != 1L || !is.finite(m) || m <= 0) {
    stop("m must be one positive finite number", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be one finite number strictly between zero and one",
         call. = FALSE)
  }
  if (length(margin) != 1L || !is.finite(margin) || margin <= 0) {
    stop("margin must be one positive finite number", call. = FALSE)
  }
  floor_iterations <- ceiling(m / alpha) - 1
  as.integer(ceiling(margin * (floor_iterations + 1)))
}

# Frequency-based null for node strength: total non-self transitions in +
# out of each node (raw counts, not the row-normalized probabilities used
# elsewhere), compared against the same quantity computed on each of the
# shuffle-null count matrices already generated for edge-level significance
# testing (see MarkovWhistle's null_model = "shuffle", R/markov_whistle.R).
# That permutation preserves each whistle sub-category's total emission
# count exactly, so a node's null strength distribution reflects only how
# much connectivity its own emission frequency predicts by chance - deviation
# from that (via the returned z_score) is evidence of above/below-chance
# participation in transitions, independent of which individual transitions
# later turn out to be BH-significant.
#
# Uses raw counts rather than the row-normalized transition_probabilities_
# matrix because row-normalization forces every node's out-strength to ~1 by
# construction (each row sums to 1), which would swamp the comparison with a
# constant, uninformative term.
#
# exclude_self = FALSE keeps self-transitions in the sum (diagnostic use:
# comparing the two reveals how much of a node's below/above-null deviation
# is driven by self-repetition rather than connections to other nodes).
strength_null_distribution <- function(
    transition_counts_matrix, null_count_matrices, progress_every = 20000L,
    exclude_self = TRUE
) {
  if (!is.matrix(transition_counts_matrix) ||
      nrow(transition_counts_matrix) != ncol(transition_counts_matrix)) {
    stop("transition_counts_matrix must be a square matrix", call. = FALSE)
  }
  n <- nrow(transition_counts_matrix)
  iterations <- length(null_count_matrices)
  if (iterations < 2L) {
    stop("At least two null matrices are required", call. = FALSE)
  }

  node_strength <- function(m) {
    if (exclude_self) diag(m) <- 0
    rowSums(m) + colSums(m)
  }

  observed_strength <- node_strength(transition_counts_matrix)

  null_sum <- numeric(n)
  null_sumsq <- numeric(n)
  exceedance_count <- integer(n)
  for (iteration in seq_len(iterations)) {
    s <- node_strength(null_count_matrices[[iteration]])
    null_sum <- null_sum + s
    null_sumsq <- null_sumsq + s^2
    exceedance_count <- exceedance_count + (s >= observed_strength)
    if (progress_every > 0L &&
        (iteration %% progress_every == 0L || iteration == iterations)) {
      message("Null node-strength pass: ", iteration, " / ", iterations)
    }
  }

  null_mean <- null_sum / iterations
  null_var <- pmax(0, null_sumsq / iterations - null_mean^2) *
    iterations / (iterations - 1)
  null_sd <- sqrt(null_var)
  z_score <- ifelse(
    null_sd > 0, (observed_strength - null_mean) / null_sd, NA_real_
  )
  # Same empirical-p formula as empirical_enrichment_p_values() (ties
  # included, Laplace-smoothed so p is never exactly 0): rank-based, so it
  # doesn't assume the null is normally distributed the way the z-score
  # implicitly does - matters for low-count nodes whose null distribution
  # is visibly discrete/skewed rather than bell-shaped.
  p_value <- (1 + exceedance_count) / (iterations + 1)

  list(
    observed_strength = observed_strength,
    null_mean = null_mean,
    null_sd = null_sd,
    z_score = z_score,
    exceedance_count = exceedance_count,
    p_value = p_value,
    iterations = iterations
  )
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

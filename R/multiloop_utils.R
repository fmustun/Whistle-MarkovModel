# Groups temporally adjacent whistles into "multi-loop" chains using a tight
# inter-whistle interval (IWI) threshold, independent of the longer-window
# Markov transition analysis in R/markov_whistle.R.

# Assigns each whistle to a chain: within a recording, ordered by start_time,
# a new chain begins whenever the gap to the previous whistle
# (start_time[i] - end_time[i-1]) is >= threshold. Chains are purely
# temporal - no requirement that chained whistles share a type/category.
compute_multiloop_chains <- function(
    whistles_list,
    threshold = MULTILOOP_IWI_THRESHOLD,
    list_names = LIST_NAMES
) {
  seq_category <- build_seq_category(whistles_list, list_names)
  type_to_category <- rep(NA_character_, max(whistles_list$whistle_type))
  for (cat_name in list_names) {
    type_to_category[seq_category[[cat_name]]] <- cat_name
  }

  recordings <- unique(whistles_list$recording)
  chain_frames <- vector("list", length(recordings))

  for (i in seq_along(recordings)) {
    rec <- recordings[i]
    video <- whistles_list[whistles_list$recording == rec, ]
    video <- video[order(video$start_time), ]
    n <- nrow(video)
    if (n == 0) next

    gap <- c(Inf, video$start_time[-1] - video$end_time[-n])
    local_chain_id <- cumsum(gap >= threshold)

    chain_frames[[i]] <- data.frame(
      recording = video$recording,
      chain_id = paste(rec, local_chain_id, sep = "_"),
      whistle_type = video$whistle_type,
      whistle_type_chr = video$whistle_type_chr,
      category = type_to_category[video$whistle_type],
      start_time = video$start_time,
      end_time = video$end_time,
      stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, chain_frames)
  result <- result[order(result$recording, result$start_time), ]

  result$chain_length <- as.integer(
    ave(result$whistle_type, result$chain_id, FUN = length)
  )
  result$position_in_chain <- ave(
    seq_len(nrow(result)), result$chain_id, FUN = seq_along
  )

  rownames(result) <- NULL
  result
}

# Collapses a compute_multiloop_chains() whistle-level table to one row per
# chain: recording, chain length, the category/type of the whistle that
# starts the chain, and the chain's overall start/end time.
summarize_multiloop_chains <- function(whistle_level_df) {
  ordered_df <- whistle_level_df[
    order(whistle_level_df$recording, whistle_level_df$start_time),
  ]
  first_idx <- !duplicated(ordered_df$chain_id)
  last_idx <- !duplicated(ordered_df$chain_id, fromLast = TRUE)

  data.frame(
    chain_id = ordered_df$chain_id[first_idx],
    recording = ordered_df$recording[first_idx],
    chain_length = ordered_df$chain_length[first_idx],
    leading_category = ordered_df$category[first_idx],
    leading_whistle_type_chr = ordered_df$whistle_type_chr[first_idx],
    start_time = ordered_df$start_time[first_idx],
    end_time = ordered_df$end_time[last_idx],
    stringsAsFactors = FALSE
  )
}

# Replicates the transition-pair selection MarkovWhistle() (R/markov_whistle.R)
# uses to build transition_counts_matrix - the empirical pair selection only,
# not the null model - and, for each recorded transition instance, flags
# whether its two whistles belong to the same multi-loop chain (see
# compute_multiloop_chains()). One row per transition instance (i.e. this is
# instance-level, matching how transition_counts_matrix is incremented, not
# deduplicated by whistle-type pair). Also carries from_row_id/to_row_id -
# each instance's row position in the whistles_list passed in (1-based,
# assigned fresh here so it's correct regardless of whistles_list's own row
# names) - so a caller can look a transition instance's whistle_type back up
# under a permutation of whistles_list$whistle_type (see
# compute_permutation_transition_significance(), used by
# scripts/08_standalone_network.R).
compute_markov_transition_multiloop_overlap <- function(
    whistles_list,
    time_window = TIME_WINDOW,
    threshold = MULTILOOP_IWI_THRESHOLD,
    list_names = LIST_NAMES
) {
  whistle_chains <- compute_multiloop_chains(
    whistles_list, threshold = threshold, list_names = list_names
  )
  whistles_list$row_id <- seq_len(nrow(whistles_list))

  recordings <- unique(whistles_list$recording)
  transition_frames <- vector("list", length(recordings))

  for (i in seq_along(recordings)) {
    rec <- recordings[i]
    video <- whistles_list[whistles_list$recording == rec, ]
    video <- video[order(video$start_time), ]
    chain_ids <- whistle_chains$chain_id[whistle_chains$recording == rec]
    n <- nrow(video)
    if (n < 2) next

    rows <- vector("list", n - 1L)
    for (whi in seq_len(n - 1L)) {
      sta_t <- video$start_time[whi]
      following_whistles <- video$whistle_type[
        sta_t + time_window[1] < video$start_time &
          video$start_time < sta_t + time_window[2]
      ]
      if (length(following_whistles) == 0) next

      time_interval <- video$start_time[whi + 1L] - video$end_time[whi]
      to_position <- whi + 1L

      # Mirrors MarkovWhistle(): an overlapping immediate-next whistle
      # (negative gap) is skipped in favor of the whistle after it, when one
      # exists within the window; otherwise no transition is recorded at all
      # for this whi (see R/markov_whistle.R).
      if (time_interval < 0) {
        if (length(following_whistles) > 1 && whi + 2L <= n) {
          to_position <- whi + 2L
          time_interval <- video$start_time[whi + 2L] - video$end_time[whi]
        } else {
          next
        }
      }

      rows[[whi]] <- data.frame(
        recording = rec,
        from_position = whi,
        to_position = to_position,
        from_row_id = video$row_id[whi],
        to_row_id = video$row_id[to_position],
        from_whistle_type = video$whistle_type[whi],
        to_whistle_type = video$whistle_type[to_position],
        from_whistle_type_chr = video$whistle_type_chr[whi],
        to_whistle_type_chr = video$whistle_type_chr[to_position],
        transition_interval = time_interval,
        within_multiloop_chain = chain_ids[whi] == chain_ids[to_position],
        stringsAsFactors = FALSE
      )
    }
    transition_frames[[i]] <- do.call(rbind, rows)
  }

  result <- do.call(rbind, transition_frames)
  rownames(result) <- NULL
  result
}

# Counts adjacent same-whistle_type pairs within a multi-loop chain (e.g. a
# 3-length chain A -> A -> B contributes one A self-repeat pair). Operates
# directly on compute_multiloop_chains()'s whistle-level output; singleton
# chains contribute nothing, so no separate chain_length filtering is
# needed. Returns one row per whistle_type that has at least one self-repeat.
compute_self_repeat_counts <- function(whistle_chains) {
  ordered <- whistle_chains[
    order(whistle_chains$chain_id, whistle_chains$position_in_chain),
  ]
  same_chain_next <- c(ordered$chain_id[-1], NA) == ordered$chain_id
  same_type_next <- c(ordered$whistle_type[-1], NA) == ordered$whistle_type
  is_self_repeat <- !is.na(same_chain_next) & same_chain_next & same_type_next

  counts <- table(ordered$whistle_type[is_self_repeat])
  data.frame(
    whistle_type = as.integer(names(counts)),
    n_self_repeats = as.integer(counts),
    stringsAsFactors = FALSE
  )
}

# One row per (chain_id, whistle_type) actually present in a multi-loop
# chain (chain_length >= 2), collapsing repeated positions of the same type
# within a chain (e.g. A -> A -> B contributes one A row, one B row) and
# dropping position/order entirely: "which types ever share a chain", not
# "which exact/sorted sequence does a chain spell out" (that's
# compute_multiloop_patterns()).
compute_multiloop_type_membership <- function(whistle_chains) {
  multiloop <- whistle_chains[whistle_chains$chain_length >= 2, ]
  deduped <- multiloop[!duplicated(multiloop[c("chain_id", "whistle_type")]), ]
  deduped[
    order(deduped$chain_id, deduped$whistle_type),
    c("chain_id", "whistle_type")
  ]
}

# Assigns each multi-loop chain (chain_length >= 2) its sub-category
# pattern, under three notions of "the same multi-loop": exact sequence
# (whistle_type_chr in order, e.g. "A -> A -> B"), same composition, any
# order (sorted whistle_type_chr *with* repeats, e.g. "A + A + B" - so
# "A -> A -> B" and "A -> B -> A" collapse to the same pattern here, but
# "A -> B" does not), and same composition ignoring repeat count (sorted
# *distinct* whistle_type_chr, e.g. "A + B" - so "A -> A -> B" and "A -> B"
# now collapse together too, since both ever contain exactly {A, B}). One
# row per chain_id; singleton chains are excluded.
compute_multiloop_patterns <- function(whistle_chains) {
  multiloop <- whistle_chains[whistle_chains$chain_length >= 2, ]
  ordered <- multiloop[
    order(multiloop$chain_id, multiloop$position_in_chain),
  ]
  ordered_pattern <- tapply(
    ordered$whistle_type_chr, ordered$chain_id,
    function(x) paste(x, collapse = " -> ")
  )
  unordered_pattern <- tapply(
    ordered$whistle_type_chr, ordered$chain_id,
    function(x) paste(sort(x), collapse = " + ")
  )
  composition_pattern <- tapply(
    ordered$whistle_type_chr, ordered$chain_id,
    function(x) paste(sort(unique(x)), collapse = " + ")
  )
  data.frame(
    chain_id = names(ordered_pattern),
    ordered_pattern = as.character(ordered_pattern),
    unordered_pattern = as.character(unordered_pattern[names(ordered_pattern)]),
    composition_pattern = as.character(composition_pattern[names(ordered_pattern)]),
    stringsAsFactors = FALSE
  )
}

# Ordered whistle_type_chr sequence for every multi-loop chain (chain_length
# >= 2), as a named list (chain_id -> character vector, in position_in_chain
# order) - the raw ordered material behind compute_multiloop_patterns()'s
# collapsed ordered_pattern string, needed wherever a caller wants the
# sequence itself rather than a joined string (e.g. the chain-alignment
# composition clustering in R/chain_clustering_utils.R).
extract_multiloop_chain_sequences <- function(whistle_chains) {
  multiloop <- whistle_chains[whistle_chains$chain_length >= 2, ]
  ordered <- multiloop[
    order(multiloop$chain_id, multiloop$position_in_chain),
  ]
  split(ordered$whistle_type_chr, ordered$chain_id)
}

# Empirical enrichment significance test shared by the stand-alone whistle
# network (scripts/08_standalone_network.R) and the multi-loop composition
# network (scripts/10_multiloop_composition_network.R) - the same "shuffle"
# null model + empirical-p + BH-correction idea script 01 uses for the full
# whistle-level network (MarkovWhistle() in R/markov_whistle.R,
# empirical_enrichment_p_values()/adjust_transition_p_values() in
# R/statistical_utils.R), but computed directly at the handful of cells
# that can ever matter instead of full n_nodes x n_nodes matrices - a dense
# n_nodes^2 matrix per null iteration would be mostly zero and, at the
# iteration counts BH correction needs (tens of thousands), unnecessary
# memory/compute to build and hold per iteration. Since
# adjust_transition_p_values' family_mask already
# restricts the multiple-testing family - and therefore which cells can
# ever end up significant - to cells where transition_counts_matrix > 0,
# only those cells need a p-value at all; every other cell is forced to
# p = 1 regardless.
#
# Null model: NULL_MODEL = "shuffle" in R/config.R - a global permutation of
# which label sits at which position, keeping every position's underlying
# timing (and therefore which position-pairs count as a transition, encoded
# in `from_idx`/`to_idx`) fixed. Because that pairing is fixed and
# independent of label identity, one null iteration is just a permutation
# of `labels` followed by a cheap tabulate()/match() lookup at the fixed set
# of observed edges - no need to rebuild timing structure per iteration the
# way MarkovWhistle() rebuilds whistle-level transitions.
#
# labels: integer node id per permutable unit (one per whistle in script 08;
#   one per stand-alone whistle or multi-loop chain in script 10), in a
#   fixed, arbitrary order - `sample(labels)` defines one null permutation.
# from_idx / to_idx: equal-length integer indices into `labels` - transition
#   instance k pairs labels[from_idx[k]] -> labels[to_idx[k]]. Fixed across
#   iterations; only `labels` is permuted.
#
# Alongside the edge-level test above, each iteration's row_sum (already
# computed for null_probability) plus a matching col_sum give each node's
# total in + out transition count for that iteration for free - the same
# node-strength null strength_null_distribution() (R/statistical_utils.R)
# builds from saved whistle-level null count matrices, but computed here
# from the label permutation directly since script 10 has no
# saved per-iteration count matrices to build a dense n_nodes x n_nodes
# matrix from. Self-transitions are included (row_sum/col_sum both count a
# from==to pair), matching strength_null_distribution(exclude_self = FALSE)
# - the "including self-transitions" strength already used for the
# whistle-level z-score plot (scripts/04_node_network_measures.R).
compute_permutation_transition_significance <- function(
    labels, from_idx, to_idx, n_nodes, transition_counts_matrix,
    iterations, progress_every = 5000L
) {
  observed_edges <- which(transition_counts_matrix > 0, arr.ind = TRUE)
  from_node_id <- observed_edges[, 1]
  to_node_id <- observed_edges[, 2]
  observed_key <- from_node_id + (to_node_id - 1L) * n_nodes
  observed_probability_matrix <- CalcTransitionMatrix(transition_counts_matrix)
  observed_probability <- observed_probability_matrix[observed_edges]
  observed_node_strength <- rowSums(transition_counts_matrix) +
    colSums(transition_counts_matrix)

  exceedance_count <- integer(length(observed_key))
  strength_sum <- numeric(n_nodes)
  strength_sumsq <- numeric(n_nodes)
  strength_exceedance_count <- integer(n_nodes)

  for (iteration in seq_len(iterations)) {
    shuffled_labels <- sample(labels)
    from_node <- shuffled_labels[from_idx]
    to_node <- shuffled_labels[to_idx]

    row_sum <- tabulate(from_node, nbins = n_nodes)
    col_sum <- tabulate(to_node, nbins = n_nodes)
    cell_count <- tabulate(
      match(from_node + (to_node - 1L) * n_nodes, observed_key),
      nbins = length(observed_key)
    )
    null_probability <- ifelse(
      row_sum[from_node_id] == 0, 0, cell_count / row_sum[from_node_id]
    )
    exceedance_count <- exceedance_count +
      (null_probability >= observed_probability)

    node_strength_iteration <- row_sum + col_sum
    strength_sum <- strength_sum + node_strength_iteration
    strength_sumsq <- strength_sumsq + node_strength_iteration^2
    strength_exceedance_count <- strength_exceedance_count +
      (node_strength_iteration >= observed_node_strength)

    if (progress_every > 0L &&
        (iteration %% progress_every == 0L || iteration == iterations)) {
      message("Empirical p-value pass: ", iteration, " / ", iterations)
    }
  }

  p_value <- (1 + exceedance_count) / (iterations + 1)

  p_value_matrix <- matrix(1, n_nodes, n_nodes)
  p_value_matrix[observed_edges] <- p_value
  exceedance_count_matrix <- matrix(0L, n_nodes, n_nodes)
  exceedance_count_matrix[observed_edges] <- exceedance_count

  strength_null_mean <- strength_sum / iterations
  strength_null_var <- pmax(0, strength_sumsq / iterations - strength_null_mean^2) *
    iterations / (iterations - 1)
  strength_null_sd <- sqrt(strength_null_var)
  strength_z_score <- ifelse(
    strength_null_sd > 0,
    (observed_node_strength - strength_null_mean) / strength_null_sd,
    NA_real_
  )
  strength_p_value <- (1 + strength_exceedance_count) / (iterations + 1)

  list(
    p_value_matrix = p_value_matrix,
    exceedance_count_matrix = exceedance_count_matrix,
    observed_probability_matrix = observed_probability_matrix,
    iterations = iterations,
    observed_node_strength = observed_node_strength,
    strength_null_mean = strength_null_mean,
    strength_null_sd = strength_null_sd,
    strength_z_score = strength_z_score,
    strength_exceedance_count = strength_exceedance_count,
    strength_p_value = strength_p_value
  )
}

# Thin wrapper of compute_permutation_transition_significance() for the
# multi-loop composition network (scripts/10_multiloop_composition_network.R):
# a transition is always between a unit (a stand-alone whistle or a whole
# multi-loop chain) and the next unit in start_time order, so
# from_idx/to_idx reduce to `which(within_window)` and `from_idx + 1`.
#
# chain_node_ids: integer node id per unit, in the same recording/start_time
#   order used to build `within_window`.
# within_window: logical, length(chain_node_ids) - 1; TRUE at position i
#   when unit i is recorded as transitioning to unit i + 1.
compute_multiloop_type_significance <- function(
    chain_node_ids, within_window, n_nodes, transition_counts_matrix,
    iterations, progress_every = 5000L
) {
  from_idx <- which(within_window)
  to_idx <- from_idx + 1L
  compute_permutation_transition_significance(
    chain_node_ids, from_idx, to_idx, n_nodes, transition_counts_matrix,
    iterations, progress_every
  )
}

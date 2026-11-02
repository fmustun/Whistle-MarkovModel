# Ordered-categorical distance machinery for multi-loop chain composition
# clustering (scripts/20, scripts/21): a substitution-cost matrix over whistle
# sub-categories (whistle_type_chr) that makes same-main-category swaps
# cheap and cross-category swaps expensive, plus a Needleman-Wunsch global
# alignment distance built on top of it. Kept separate from
# R/multiloop_utils.R, which builds the chains themselves
# (compute_multiloop_chains()) - this file only knows about ordered label
# sequences and a cost matrix, not about whistles/recordings/timing.

# One substitution cost per pair of sub-categories: 0 if identical, epsilon
# if they share a main category (`category_of_sub`, e.g. the LIST_NAMES
# grouping), cross_category_cost otherwise. Symmetric by construction
# (identity and "same category" are both symmetric relations), and satisfies
# the triangle inequality as long as epsilon <= cross_category_cost (checked
# below) - because "same category" is an equivalence relation, no pair of
# epsilon-hops can ever undercut a cross-category hop, since two
# sub-categories that are each within epsilon of a common third must share
# that third's category, and therefore each other's.
#
# cross_category_cost defaults to 1 (paired with gap_penalty = 1 elsewhere -
# see nw_alignment_cost() - a cross-category substitution and a gap cost
# exactly the same, so alignment breaks ties between them arbitrarily).
# Setting cross_category_cost below gap_penalty (e.g. 0.6) makes every
# cross-category substitution strictly cheaper than a gap, so alignment will
# now always prefer substituting across categories over leaving either
# sequence unaligned - a real behavior change from the 0/eps/1 scheme, not
# just a rescaling, worth knowing about when comparing cluster counts
# between the two.
build_chain_substitution_cost <- function(
    sub_levels, category_of_sub, epsilon = 0.3, cross_category_cost = 1
) {
  if (length(cross_category_cost) != 1L || !is.finite(cross_category_cost) ||
      cross_category_cost < 0) {
    stop("cross_category_cost must be one finite non-negative number", call. = FALSE)
  }
  if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon < 0 ||
      epsilon > cross_category_cost) {
    stop("epsilon must be one finite number in [0, cross_category_cost]", call. = FALSE)
  }
  category_of_sub <- category_of_sub[sub_levels]
  if (anyNA(category_of_sub)) {
    stop("category_of_sub is missing an entry for some sub_levels", call. = FALSE)
  }
  same_category <- outer(category_of_sub, category_of_sub, "==")
  cost <- matrix(
    cross_category_cost, length(sub_levels), length(sub_levels),
    dimnames = list(sub_levels, sub_levels)
  )
  cost[same_category] <- epsilon
  diag(cost) <- 0
  cost
}

# Raw Needleman-Wunsch global alignment cost between two label sequences,
# encoded as 1-based row/column indices into `cost_matrix` (so callers doing
# many pairwise alignments over the same label set can encode every sequence
# once up front - see chain_alignment_distance_matrix() below, rather than
# re-matching character labels on every pair). Standard O(m*n) edit-distance
# DP: substitution cost from `cost_matrix`, linear gap cost `gap_penalty` for
# either sequence's unaligned positions.
nw_alignment_cost <- function(seq_a_idx, seq_b_idx, cost_matrix, gap_penalty = 1) {
  m <- length(seq_a_idx)
  n <- length(seq_b_idx)
  if (m == 0L && n == 0L) return(0)
  if (m == 0L) return(n * gap_penalty)
  if (n == 0L) return(m * gap_penalty)

  dp <- matrix(0, m + 1L, n + 1L)
  dp[, 1L] <- (0:m) * gap_penalty
  dp[1L, ] <- (0:n) * gap_penalty
  for (i in seq_len(m)) {
    sub_cost_row <- cost_matrix[seq_a_idx[i], seq_b_idx]
    dp_prev_row <- dp[i, ]
    for (j in seq_len(n)) {
      dp[i + 1L, j + 1L] <- min(
        dp_prev_row[j] + sub_cost_row[j],
        dp[i, j + 1L] + gap_penalty,
        dp[i + 1L, j] + gap_penalty
      )
    }
  }
  dp[m + 1L, n + 1L]
}

# Same DP as nw_alignment_cost(), but returns the *operation counts* (match,
# substitution, gap) along one optimal alignment path instead of just its
# total cost - the input mdl_k_sweep() below needs to build a description
# length from (Rissanen's two-part MDL, generalized here from Gaussian/
# real-valued data to categorical edit operations): "how many bits does
# encoding chain_a as an edit script against reference chain_b actually
# cost", which depends on which TYPE of operation happened at each position,
# not just the cost_matrix's cost for it. Ties are broken toward a diagonal
# move (match/substitution) over an indel whenever they tie on cost - one
# arbitrary but deterministic choice among possibly several optimal
# alignments, same as any traceback-based aligner has to make.
#
# `seq_a` is treated as "the sequence being encoded" and `seq_b` as "the
# reference/medoid" (matching how mdl_k_sweep() calls this): a gap can mean
# two different things for that purpose, so this distinguishes them - an
# INSERTION (seq_a has a symbol the alignment doesn't match to anything in
# seq_b) has to tell the decoder *which* label that extra seq_a symbol is; a
# DELETION (seq_b/the medoid has a symbol the alignment skips over in seq_a)
# doesn't - the decoder already knows the medoid's own label there, it's
# just being told "skip it", no further identity information needed. Same
# distinction for substitution vs. match: a substitution's target label
# (seq_a's own symbol) also has to be told to the decoder; a match doesn't
# (it's already the medoid's label by definition).
#
# Returns a list:
#  - counts: c(match=, substitution=, insertion=, deletion=).
#  - substitution_targets, insertion_targets: the seq_a_idx label at each
#    substitution/insertion respectively (one entry per occurrence, in
#    traceback order) - the label identities a caller needs to build an
#    adaptive "which label" code from, rather than assume a uniform one.
nw_alignment_ops <- function(seq_a_idx, seq_b_idx, cost_matrix, gap_penalty = 1) {
  m <- length(seq_a_idx)
  n <- length(seq_b_idx)
  if (m == 0L && n == 0L) {
    return(list(
      counts = c(match = 0L, substitution = 0L, insertion = 0L, deletion = 0L),
      substitution_targets = integer(0), insertion_targets = integer(0)
    ))
  }

  dp <- matrix(0, m + 1L, n + 1L)
  dp[, 1L] <- (0:m) * gap_penalty
  dp[1L, ] <- (0:n) * gap_penalty
  # move[i+1, j+1]: 1 = diagonal (match/substitution), 2 = up (insertion:
  # seq_a has an unmatched symbol), 3 = left (deletion: seq_b has one).
  move <- matrix(0L, m + 1L, n + 1L)
  if (m > 0L) move[-1L, 1L] <- 2L
  if (n > 0L) move[1L, -1L] <- 3L

  if (m > 0L && n > 0L) {
    for (i in seq_len(m)) {
      sub_cost_row <- cost_matrix[seq_a_idx[i], seq_b_idx]
      dp_prev_row <- dp[i, ]
      for (j in seq_len(n)) {
        diag_cost <- dp_prev_row[j] + sub_cost_row[j]
        up_cost <- dp[i, j + 1L] + gap_penalty
        left_cost <- dp[i + 1L, j] + gap_penalty
        best <- min(diag_cost, up_cost, left_cost)
        dp[i + 1L, j + 1L] <- best
        move[i + 1L, j + 1L] <- if (diag_cost == best) 1L else if (up_cost == best) 2L else 3L
      }
    }
  }

  i <- m; j <- n
  n_match <- 0L; n_sub <- 0L; n_ins <- 0L; n_del <- 0L
  substitution_targets <- integer(0)
  insertion_targets <- integer(0)
  while (i > 0L || j > 0L) {
    mv <- move[i + 1L, j + 1L]
    if (mv == 1L) {
      if (seq_a_idx[i] == seq_b_idx[j]) {
        n_match <- n_match + 1L
      } else {
        n_sub <- n_sub + 1L
        substitution_targets <- c(substitution_targets, seq_a_idx[i])
      }
      i <- i - 1L; j <- j - 1L
    } else if (mv == 2L) {
      n_ins <- n_ins + 1L
      insertion_targets <- c(insertion_targets, seq_a_idx[i])
      i <- i - 1L
    } else {
      n_del <- n_del + 1L
      j <- j - 1L
    }
  }
  list(
    counts = c(match = n_match, substitution = n_sub, insertion = n_ins, deletion = n_del),
    substitution_targets = substitution_targets, insertion_targets = insertion_targets
  )
}

# nw_alignment_cost() normalized by the longer sequence's length, so chain
# pairs of very different lengths aren't automatically far apart purely from
# gap accumulation. normalize = FALSE returns the raw, unnormalized cost
# instead - mdl_k_sweep() below needs this for medoid selection, since
# picking the member that minimizes total *normalized* distance to its
# cluster-mates does not, in general, also minimize total raw edit cost
# (the quantity the description-length accounting actually measures) once
# cluster members span more than one chain length.
nw_alignment_distance <- function(seq_a_idx, seq_b_idx, cost_matrix, gap_penalty = 1, normalize = TRUE) {
  cost <- nw_alignment_cost(seq_a_idx, seq_b_idx, cost_matrix, gap_penalty)
  if (!normalize) return(cost)
  longer <- max(length(seq_a_idx), length(seq_b_idx))
  if (longer == 0L) return(0)
  cost / longer
}

# Pairwise normalized NW-alignment distance across a named list of label
# sequences (character vectors over rownames(cost_matrix)) - e.g. one entry
# per *distinct* multi-loop composition after deduping identical chains (see
# cluster_multiloop_chains() below). Returns a `dist` object, ready for
# hclust()/cluster::silhouette().
# Encodes every sequence to integer indices once up front rather than per
# pair, since character-label cost_matrix lookups are the dominant avoidable
# cost once there are hundreds of sequences and tens of thousands of pairs.
# normalize = FALSE returns raw (unnormalized) pairwise costs instead - see
# nw_alignment_distance()'s normalize argument.
#
# Pairs are generated with utils::combn(n, 2), which enumerates (1,2), (1,3),
# ..., (1,n), (2,3), ..., (n-1,n) - the same order `dist` objects store their
# lower triangle in, so `distances` can be turned directly into a `dist`
# without ever materializing the full n x n matrix.
chain_alignment_distance_matrix <- function(
    chain_sequences, cost_matrix, gap_penalty = 1, normalize = TRUE, progress_every = 20000L
) {
  n <- length(chain_sequences)
  if (n < 2L) {
    stop("Need at least two chain sequences to build a distance matrix", call. = FALSE)
  }
  sub_levels <- rownames(cost_matrix)
  encoded <- lapply(chain_sequences, function(seq) {
    idx <- match(seq, sub_levels)
    if (anyNA(idx)) {
      stop("chain_sequences contains a label not in rownames(cost_matrix)", call. = FALSE)
    }
    idx
  })

  pair_idx <- utils::combn(n, 2)
  n_pairs <- ncol(pair_idx)
  distances <- numeric(n_pairs)
  for (p in seq_len(n_pairs)) {
    i <- pair_idx[1L, p]
    j <- pair_idx[2L, p]
    distances[p] <- nw_alignment_distance(
      encoded[[i]], encoded[[j]], cost_matrix, gap_penalty, normalize = normalize
    )
    if (progress_every > 0L && (p %% progress_every == 0L || p == n_pairs)) {
      message("Chain alignment distance pass: ", p, " / ", n_pairs)
    }
  }

  attributes(distances) <- list(
    Size = n, Labels = names(chain_sequences), Diag = FALSE, Upper = FALSE,
    method = if (normalize) "nw_alignment" else "nw_alignment_raw", class = "dist"
  )
  distances
}

# Mean cluster::silhouette() width for every k in `k_range`, cutting the same
# hclust tree at each k - the model-selection sweep cluster_multiloop_chains()
# below uses to pick a cluster count instead of fixing one by hand. Requires
# the `cluster` package. Singleton clusters (cutree can produce them at
# larger k) get silhouette width 0 by cluster::silhouette()'s own
# convention, which penalizes k values that fragment the tree into trivial
# clusters.
silhouette_k_sweep <- function(hc, dist_obj, k_range) {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("The 'cluster' package is required for silhouette_k_sweep()", call. = FALSE)
  }
  mean_width <- vapply(k_range, function(k) {
    clusters <- stats::cutree(hc, k = k)
    sil <- cluster::silhouette(clusters, dist_obj)
    mean(sil[, "sil_width"])
  }, numeric(1))
  data.frame(k = k_range, mean_silhouette_width = mean_width)
}

# Dedupes identical chain sequences (carrying counts) and builds the
# pairwise NW-alignment distance matrix + hclust tree over the distinct
# compositions - the first half of cluster_multiloop_chains() below, kept as
# its own function for readability.
align_and_cluster_chains <- function(
    chain_sequences, cost_matrix, gap_penalty = 1, linkage_method = "average",
    progress_every = 20000L
) {
  if (length(chain_sequences) < 2L) {
    stop("Need at least two chain sequences to cluster", call. = FALSE)
  }

  pattern_key <- vapply(chain_sequences, paste, collapse = " -> ", FUN.VALUE = character(1))
  chain_ids_by_pattern <- split(names(chain_sequences), pattern_key)
  ordered_patterns <- names(chain_ids_by_pattern)
  n_patterns <- length(ordered_patterns)
  pattern_id <- sprintf("pattern_%0*d", nchar(n_patterns), seq_len(n_patterns))
  pattern_lookup <- setNames(pattern_id, ordered_patterns)

  representative_chain_id <- vapply(chain_ids_by_pattern, `[[`, 1L, FUN.VALUE = character(1))
  pattern_sequences <- setNames(chain_sequences[representative_chain_id], pattern_id)
  n_chains_per_pattern <- setNames(lengths(chain_ids_by_pattern), pattern_id)

  dist_obj <- chain_alignment_distance_matrix(
    pattern_sequences, cost_matrix, gap_penalty = gap_penalty,
    progress_every = progress_every
  )
  hc <- hclust(dist_obj, method = linkage_method)

  list(
    pattern_key = pattern_key,
    pattern_lookup = pattern_lookup,
    pattern_sequences = pattern_sequences,
    n_chains_per_pattern = n_chains_per_pattern,
    dist_obj = dist_obj,
    hc = hc
  )
}

# End-to-end composition clustering for a compute_multiloop_chains() table:
# dedupe identical chains (carrying counts), pairwise NW-alignment distance
# over the distinct compositions, hierarchical clustering, and k chosen by
# max mean silhouette width over k_range. One shared implementation so every
# caller building on these clusters (scripts/20's sub-clustering,
# scripts/21's cluster-level transition network) gets identical,
# reproducible assignments from the same parameters, rather than re-deriving
# the pipeline (and risking drift from R/config.R's own MULTILOOP_CLUSTER_
# EPSILON/GAP_PENALTY/K_RANGE) by hand. Deterministic - no seed needed -
# since hclust/cutree/silhouette have no random component.
#
# Returns a list:
#  - chain_level: one row per multi-loop chain (chain_id, ordered_pattern,
#    chain_length, pattern_id, cluster).
#  - pattern_table: one row per distinct composition (pattern_id,
#    ordered_pattern, pattern_length, n_chains, cluster).
#  - dist_obj, hc: the NW-alignment `dist` object and its hclust tree, over
#    pattern_table's distinct compositions (for further plotting/diagnostics).
#  - silhouette_sweep, best_k, best_width, k_range: the k-selection sweep and
#    its chosen optimum (or, if `force_k` was given, the sweep is still
#    computed for reference but best_k/best_width instead reflect the forced
#    cut).
#  - chain_sequences: the un-deduped, per-chain label sequences
#    (extract_multiloop_chain_sequences()'s own output) - kept on the return
#    value for callers that need the raw sequences alongside the clustering.
#
# force_k: NULL (default) selects k by silhouette as described above. Set to
# an integer to cut the SAME hclust tree at that k instead - e.g. to inspect
# the k a different criterion preferred (mdl_k_sweep()) without re-deriving
# the tree by hand. Must be in [2, n_patterns).
cluster_multiloop_chains <- function(
    whistle_chains, cost_matrix, gap_penalty = 1, linkage_method = "average",
    k_range = 5:40, progress_every = 20000L, force_k = NULL
) {
  chain_sequences <- extract_multiloop_chain_sequences(whistle_chains)
  aligned <- align_and_cluster_chains(
    chain_sequences, cost_matrix, gap_penalty = gap_penalty,
    linkage_method = linkage_method, progress_every = progress_every
  )
  pattern_key <- aligned$pattern_key
  pattern_lookup <- aligned$pattern_lookup
  pattern_sequences <- aligned$pattern_sequences
  n_chains_per_pattern <- aligned$n_chains_per_pattern
  dist_obj <- aligned$dist_obj
  hc <- aligned$hc
  pattern_id <- names(pattern_sequences)
  ordered_patterns <- names(pattern_lookup)
  n_patterns <- length(pattern_id)

  valid_k_range <- k_range[k_range >= 2 & k_range < n_patterns]
  if (length(valid_k_range) == 0L) {
    stop("No valid k in k_range for ", n_patterns, " distinct compositions",
         call. = FALSE)
  }
  silhouette_sweep <- silhouette_k_sweep(hc, dist_obj, k_range = valid_k_range)

  if (is.null(force_k)) {
    best_k <- silhouette_sweep$k[which.max(silhouette_sweep$mean_silhouette_width)]
    best_width <- max(silhouette_sweep$mean_silhouette_width)
  } else {
    if (length(force_k) != 1L || !is.finite(force_k) ||
        force_k < 2L || force_k >= n_patterns) {
      stop("force_k must be one integer in [2, ", n_patterns, ")", call. = FALSE)
    }
    best_k <- as.integer(force_k)
    best_width <- if (best_k %in% valid_k_range) {
      silhouette_sweep$mean_silhouette_width[silhouette_sweep$k == best_k]
    } else {
      NA_real_
    }
  }
  pattern_cluster <- stats::cutree(hc, k = best_k)

  pattern_table <- data.frame(
    pattern_id = pattern_id,
    ordered_pattern = ordered_patterns,
    pattern_length = lengths(pattern_sequences),
    n_chains = as.integer(n_chains_per_pattern[pattern_id]),
    cluster = as.integer(pattern_cluster[pattern_id]),
    stringsAsFactors = FALSE
  )
  pattern_table <- pattern_table[order(-pattern_table$n_chains), ]

  chain_ids_all <- names(chain_sequences)
  chain_level <- data.frame(
    chain_id = chain_ids_all,
    ordered_pattern = unname(pattern_key[chain_ids_all]),
    chain_length = as.integer(lengths(chain_sequences[chain_ids_all])),
    stringsAsFactors = FALSE
  )
  chain_level$pattern_id <- unname(pattern_lookup[chain_level$ordered_pattern])
  chain_level$cluster <- as.integer(pattern_cluster[chain_level$pattern_id])
  chain_level <- chain_level[order(chain_level$cluster, chain_level$chain_id), ]

  list(
    chain_level = chain_level,
    pattern_table = pattern_table,
    dist_obj = dist_obj,
    hc = hc,
    silhouette_sweep = silhouette_sweep,
    best_k = best_k,
    best_width = best_width,
    k_range = valid_k_range,
    linkage_method = linkage_method,
    gap_penalty = gap_penalty,
    chain_sequences = chain_sequences
  )
}

# Two-part minimum description length (Rissanen) for choosing k - how
# many bits does it take to encode all of the (weighted) chain population,
# given a k-cluster partition? Total = L(model) + L(data | model).
#
# The partition itself is re-optimized at every k by PAM/k-medoids
# (cluster::pam() on the raw, unnormalized alignment-cost dissimilarity - see
# nw_alignment_distance()'s normalize argument for why raw, not
# clustering$dist_obj's length-normalized distance), NOT by cutting one
# fixed hclust tree at different heights. An earlier version cut
# clustering$hc instead, which is cheaper but not a fair comparison across k:
# average-linkage merges to minimize *average* pairwise distance, which can
# easily separate a good shared medoid from a chunk of points that were all
# individually close to it - so cutting deeper could raise the raw-cost fit
# at some k even though nothing about the underlying data changed, showing
# up as a spurious jump in the curve. PAM instead re-optimizes both the
# partition AND each cluster's medoid jointly for that specific k (a local
# search, not guaranteed globally optimal, but re-run fresh each time
# instead of inherited from one shared tree), which is what actually makes
# "more clusters can only help" a reasonable expectation for data_bits.
#
# A first version of this function charged every cluster a flat
# BIC_PARAMS_PER_CLUSTER = 2 penalty (for a fitted 3-way match/substitution/
# gap split) and coded *which* label a substitution/gap involved with a
# fixed uniform code (log2(n_sub - 1) / log2(n_sub), never fitted, so never
# charged for). That flat 2-parameter count turned out to never be able to
# create an interior optimum: it's the same tiny handful of bits regardless
# of how much real, recoverable structure a finer k actually explains, so it
# can only ever compete with data_bits' savings at the very smallest k. This
# version instead makes the "which label" code adaptive too - fit from each
# cluster's own empirical target-label frequencies rather than assumed
# uniform - and charges BIC for *that* fitted distribution's own degrees of
# freedom, which scale with how many distinct labels a cluster's
# substitutions/insertions actually involve. A small, homogeneous cluster
# (few distinct target labels) stays cheap either way; a large, heterogeneous
# one (many distinct target labels - which is exactly what makes it
# heterogeneous) now costs real model bits for that diversity, not just an
# undifferentiated per-cluster flat fee. This is a genuine complexity
# accounting fix, not a bigger multiplier on the same flat term.
#
#  - L(model), per cluster:
#     - length_bits + length(medoid) * log2(n_sub): the medoid composition
#       PAM chose, encoded from scratch (length_bits: a fixed-width code
#       wide enough for the longest chain in the data; each of length(medoid)
#       positions drawn from the n_sub possible sub-categories) - the
#       literal payload cost of writing down that reference sequence.
#     - bic_params_r * log2(n_bic): a BIC-style penalty for this cluster's
#       own fitted distributions (n_bic = the total number of chains, the
#       actual independent sampling unit here). bic_params_r = 3 (a 4-way
#       match/substitution/insertion/deletion split has 3 degrees of
#       freedom - see nw_alignment_ops()'s insertion/deletion distinction)
#       + (number of distinct substitution-target labels observed in this
#       cluster, minus 1, floored at 0) + (same for insertion-target
#       labels). The entropy codes below get to assume the decoder already
#       knows all of these fitted frequencies; BIC's classical p * log(n)
#       term is the standard asymptotic charge for having to estimate (and
#       transmit) them in the first place.
#  - L(data | model): within each cluster separately, every member pattern
#    is aligned to that cluster's own medoid (nw_alignment_ops()) to get
#    match/substitution/insertion/deletion counts plus which label was
#    involved at each substitution/insertion, weighted by how many actual
#    chains share that pattern (encoding the same chain 73 times costs 73x,
#    not 1x). Three separate optimal codes, each built from that cluster's
#    own empirical frequencies (Shannon's source coding theorem - using a
#    type's own observed frequency as its code length is provably the
#    cheapest possible categorical code): N_r * H(p_match, p_substitution,
#    p_insertion, p_deletion) for the operation type, plus the entropy of
#    the substitution-target label distribution (only over labels that
#    actually occur - deletion needs no target code at all, per
#    nw_alignment_ops()'s own distinction). Coding each cluster separately
#    (rather than pooling frequencies globally across every cluster first)
#    matters here: entropy of a pooled mixture is always >= the
#    size-weighted average of each part's own entropy - i.e. a shared
#    global code can only ever be equally or less efficient than letting
#    each cluster's code fit its own local mix, never more.
#
# Small k -> few, cheap models, but each medoid is a poor stand-in for its
# (heterogeneous) cluster, so match probability is low and data cost is
# high. Large k -> many, more expensive models, but each medoid is an
# almost-exact stand-in for its (near-homogeneous) cluster, so data cost
# falls toward 0. best_k minimizes the sum - unlike silhouette (undefined at
# k = 1) or the Gap statistic (needs a permutation null to define a k = 1
# baseline), k = 1 (one global medoid, no split at all) is a perfectly
# well-defined point on this same curve (cluster::pam() accepts k = 1
# directly), so k_range may start there.
#
# `clustering` is a cluster_multiloop_chains() return value (only its
# pattern_table/gap_penalty are used - its hc/dist_obj are NOT, since the
# partition is re-fit by PAM here rather than inherited from that tree);
# `cost_matrix` must be the same matrix it was built with (it isn't itself
# stored on `clustering`, so it's passed again here).
mdl_k_sweep <- function(clustering, cost_matrix, k_range, progress_every = 0L) {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("The 'cluster' package is required for mdl_k_sweep()", call. = FALSE)
  }
  # Degrees of freedom for the operation-type split alone (match/
  # substitution/insertion/deletion, a 4-way categorical) - the
  # substitution-/insertion-target label distributions add their own
  # (n_distinct - 1) each, computed per cluster inside total_bits_for_k().
  BIC_OP_TYPE_DOF <- 3

  pattern_table <- clustering$pattern_table
  pattern_id <- pattern_table$pattern_id
  n_patterns <- length(pattern_id)
  n_chains_per_pattern <- pattern_table$n_chains
  n_bic <- sum(n_chains_per_pattern)

  # cluster::pam() itself requires k <= n - 1 (rejects k == n, the trivial
  # every-point-its-own-cluster case) - same upper bound
  # cluster_multiloop_chains() already enforces for its own k_range.
  k_range <- sort(unique(k_range[k_range >= 1L & k_range < n_patterns]))
  if (length(k_range) == 0L) {
    stop("No valid k in k_range for ", n_patterns, " distinct compositions", call. = FALSE)
  }

  sub_levels <- rownames(cost_matrix)
  n_sub <- length(sub_levels)
  chain_sequences_list <- setNames(
    strsplit(pattern_table$ordered_pattern, " -> ", fixed = TRUE), pattern_id
  )
  encoded <- lapply(chain_sequences_list, function(seq) match(seq, sub_levels))
  length_bits <- ceiling(log2(max(lengths(encoded)) + 1))
  gap_penalty <- clustering$gap_penalty

  raw_dist_obj <- chain_alignment_distance_matrix(
    chain_sequences_list, cost_matrix, gap_penalty = gap_penalty,
    normalize = FALSE, progress_every = progress_every
  )

  # Shannon entropy (in bits) of a non-negative count vector, weighted by
  # its own total - shared by the op-type code and both label codes below.
  # 0 for an empty/all-zero vector (nothing to encode).
  weighted_entropy_bits <- function(counts) {
    observed <- counts[counts > 0]
    total <- sum(observed)
    if (total <= 0) return(0)
    sum(observed * (-log2(observed / total)))
  }

  total_bits_for_k <- function(k) {
    pam_fit <- cluster::pam(raw_dist_obj, k = k, diss = TRUE)
    members_by_cluster <- split(seq_len(n_patterns), pam_fit$clustering)

    model_bits <- 0
    data_bits <- 0
    for (cl in seq_along(members_by_cluster)) {
      members <- members_by_cluster[[cl]]
      medoid_seq <- encoded[[pam_fit$id.med[cl]]]

      op_counts <- c(match = 0, substitution = 0, insertion = 0, deletion = 0)
      sub_target_counts <- numeric(n_sub)
      ins_target_counts <- numeric(n_sub)
      for (i in members) {
        ops <- nw_alignment_ops(encoded[[i]], medoid_seq, cost_matrix, gap_penalty)
        weight <- n_chains_per_pattern[[i]]
        op_counts <- op_counts + ops$counts * weight
        if (length(ops$substitution_targets) > 0L) {
          sub_target_counts <- sub_target_counts +
            tabulate(ops$substitution_targets, nbins = n_sub) * weight
        }
        if (length(ops$insertion_targets) > 0L) {
          ins_target_counts <- ins_target_counts +
            tabulate(ops$insertion_targets, nbins = n_sub) * weight
        }
      }

      sub_dof <- max(0L, sum(sub_target_counts > 0) - 1L)
      ins_dof <- max(0L, sum(ins_target_counts > 0) - 1L)
      bic_params_r <- BIC_OP_TYPE_DOF + sub_dof + ins_dof

      model_bits <- model_bits + length_bits + length(medoid_seq) * log2(n_sub) +
        bic_params_r * log2(n_bic)
      data_bits <- data_bits + weighted_entropy_bits(op_counts) +
        weighted_entropy_bits(sub_target_counts) + weighted_entropy_bits(ins_target_counts)
    }

    c(model_bits = model_bits, data_bits = data_bits, total_bits = model_bits + data_bits)
  }

  rows <- lapply(k_range, function(k) {
    result <- total_bits_for_k(k)
    if (progress_every > 0L) message("MDL pass: k = ", k)
    result
  })
  sweep <- data.frame(k = k_range, do.call(rbind, rows))
  rownames(sweep) <- NULL

  list(sweep = sweep, best_k = sweep$k[which.min(sweep$total_bits)])
}

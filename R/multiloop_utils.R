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
# deduplicated by whistle-type pair).
compute_markov_transition_multiloop_overlap <- function(
    whistles_list,
    time_window = TIME_WINDOW,
    threshold = MULTILOOP_IWI_THRESHOLD,
    list_names = LIST_NAMES
) {
  whistle_chains <- compute_multiloop_chains(
    whistles_list, threshold = threshold, list_names = list_names
  )

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

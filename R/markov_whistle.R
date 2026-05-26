MarkovWhistle <- function(
    data,
    time_window = c(0, 10),
    null_model = "shift",
    iterations = 1000,
    threshold = 2,
    max_shift = 1000,
    min_shift = 0,
    list_names = LIST_NAMES
) {
  recordings <- unique(data$recording)
  recording_number <- length(recordings)
  whistle_number <- max(unique(data$whistle_type))
  recording_duration <- data$recording_duration

  seq_category <- list()
  for (wh_name in list_names) {
    seq_category[[wh_name]] <- unique(data$whistle_type)[
      grepl(wh_name, unique(data$whistle_type_chr))
    ]
  }

  whistle_occurrence <- table(data$whistle_type)
  transition_counts_matrix <- matrix(0, ncol = whistle_number, nrow = whistle_number)
  random_whistle_probability_mean_matrix <- matrix(0, ncol = whistle_number, nrow = whistle_number)
  random_whistle_probability_SD_matrix <- matrix(0, ncol = whistle_number, nrow = whistle_number)
  transition_probabilities_matrix_all <- matrix(0, ncol = whistle_number, nrow = whistle_number)

  for (rec in 1:recording_number) {
    video <- data[data$recording == recordings[rec], ]
    video <- video[order(video$start_time), ]
    for (whi in 1:(dim(video)[1] - 1)) {
      sta_t <- video$start_time[whi]
      following_whistles <- video$whistle_type[
        sta_t + time_window[1] < video$start_time &
          video$start_time < sta_t + time_window[2]
      ]
      if (length(following_whistles) != 0) {
        time_interval <- video$start_time[whi + 1] - video$end_time[whi]
        next_wh <- following_whistles[1]
        current_wh_type <- video[whi, ]$whistle_type

        for (cat_name in list_names) {
          if (current_wh_type %in% seq_category[[cat_name]]) {
            current_wh_cat <- cat_name
          }
          if (next_wh %in% seq_category[[cat_name]]) {
            following_wh_cat <- cat_name
          }
        }

        if (time_interval < 0) {
          if (length(following_whistles) > 1) {
            next_wh <- following_whistles[2]
            time_interval <- video$start_time[whi + 2] - video$end_time[whi]
            for (cat_name in list_names) {
              if (next_wh %in% seq_category[[cat_name]]) {
                following_wh_cat <- cat_name
              }
            }
            transition_counts_matrix[
              video$whistle_type[whi], next_wh
            ] <- transition_counts_matrix[
              video$whistle_type[whi], next_wh
            ] + 1
          }
        } else {
          transition_counts_matrix[
            video$whistle_type[whi], next_wh
          ] <- transition_counts_matrix[
            video$whistle_type[whi], next_wh
          ] + 1
        }
      }
    }
  }

  message("Compute null model ")
  pb <- txtProgressBar(min = 1, max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  cores <- parallel::detectCores()
  cl <- makeSOCKcluster(max(1, cores - 2))
  registerDoSNOW(cl)

  Compute_random_network <- function(iteration) {
    rdm_whistle_occ_mat_current <- matrix(
      data = 0, nrow = whistle_number, ncol = whistle_number
    )
    data_copy <- data

    if (null_model == "shuffle") {
      data_copy$whistle_type <- sample(data_copy$whistle_type, replace = FALSE)
    }

    for (rec in 1:recording_number) {
      video <- data_copy[data_copy$recording == recordings[rec], ]
      if (null_model == "shift") {
        for (whi in rownames(video)) {
          shi <- runif(1, min_shift, max_shift)
          video[whi, ]$start_time <- video[whi, ]$start_time + shi -
            floor((video[whi, ]$start_time + shi) / recording_duration[rec]) *
              recording_duration[rec]
          video[whi, ]$end_time <- video[whi, ]$end_time + shi -
            floor((video[whi, ]$end_time + shi) / recording_duration[rec]) *
              recording_duration[rec]
        }
      }
      video <- video[order(video$start_time), ]
      for (whi in 1:(dim(video)[1] - 1)) {
        sta_t <- video$start_time[whi]
        following_whistles_rd <- video$whistle_type[
          sta_t + time_window[1] < video$start_time &
            video$start_time < sta_t + time_window[2]
        ]
        if (length(following_whistles_rd) != 0) {
          time_interval <- video$start_time[whi + 1] - video$end_time[whi]
          next_wh <- following_whistles_rd[1]
          current_wh_type <- video[whi, ]$whistle_type

          for (cat_name in list_names) {
            if (current_wh_type %in% seq_category[[cat_name]]) {
              current_wh_cat <- cat_name
            }
            if (next_wh %in% seq_category[[cat_name]]) {
              following_wh_cat <- cat_name
            }
          }

          if (time_interval < 0) {
            if (length(following_whistles_rd) > 1) {
              next_wh <- following_whistles_rd[2]
              time_interval <- video$start_time[whi + 2] - video$end_time[whi]
              for (cat_name in list_names) {
                if (next_wh %in% seq_category[[cat_name]]) {
                  following_wh_cat <- cat_name
                }
              }
              rdm_whistle_occ_mat_current[
                video$whistle_type[whi], next_wh
              ] <- rdm_whistle_occ_mat_current[
                video$whistle_type[whi], next_wh
              ] + 1
            }
          } else {
            rdm_whistle_occ_mat_current[
              video$whistle_type[whi], next_wh
            ] <- rdm_whistle_occ_mat_current[
              video$whistle_type[whi], next_wh
            ] + 1
          }
        }
      }
    }

    rdm_whistle_occ_mat_current
  }

  random_whistle_occurrences_matrices <- foreach(
    ite = 1:iterations, .options.snow = opts
  ) %dopar% {
    Compute_random_network(ite)
  }

  close(pb)
  stopCluster(cl)

  transition_probabilities_matrix_all <- CalcTransitionMatrix(transition_counts_matrix)
  random_whistle_transition_probabilities_matrices <- lapply(
    random_whistle_occurrences_matrices, CalcTransitionMatrix
  )

  a_mat <- matrix(0, ncol = whistle_number, nrow = whistle_number)
  b_mat <- matrix(0, ncol = whistle_number, nrow = whistle_number)

  for (i in 1:iterations) {
    a_mat <- a_mat + (
      transition_probabilities_matrix_all >
        random_whistle_transition_probabilities_matrices[[i]]
    )
    b_mat <- b_mat + (
      transition_probabilities_matrix_all <
        random_whistle_transition_probabilities_matrices[[i]]
    )
  }

  p_values <- b_mat / iterations
  inv_p_values <- (1 - p_values) * (transition_probabilities_matrix_all > 0)

  random_whistle_probability_mean_matrix <- apply(
    simplify2array(random_whistle_transition_probabilities_matrices), 1:2, mean
  )
  random_whistle_probability_SD_matrix <- apply(
    simplify2array(random_whistle_transition_probabilities_matrices), 1:2, sd
  )
  transition_probabilities_matrix_significant <- transition_probabilities_matrix_all
  transition_probabilities_matrix_significant[
    random_whistle_probability_mean_matrix +
      threshold * random_whistle_probability_SD_matrix >
      transition_probabilities_matrix_all
  ] <- 0

  message("Done")
  list(
    transition_probabilities_matrix_significant = transition_probabilities_matrix_significant,
    transition_probabilities_matrix_all = transition_probabilities_matrix_all,
    transition_counts_matrix = transition_counts_matrix,
    p_value_matrix = p_values,
    inv_p_value_matrix = inv_p_values,
    random_whistle_occurrences_matrices = random_whistle_occurrences_matrices,
    random_whistle_probability_mean_matrix = random_whistle_probability_mean_matrix,
    random_whistle_probability_SD_matrix = random_whistle_probability_SD_matrix,
    whistle_occurrence = whistle_occurrence,
    time_window = time_window,
    null_model = null_model,
    iterations_number = iterations,
    threshold = threshold,
    max_shifth_value = max_shift,
    min_shift_value = min_shift
  )
}

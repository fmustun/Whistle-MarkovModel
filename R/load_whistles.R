load_whistles <- function(path = DATA_PATH) {
  if (!file.exists(path)) {
    stop(
      "Whistle data not found at: ", path, "\n",
      "See data/README.md for setup instructions.",
      call. = FALSE
    )
  }
  whistles_list <- read.table(path, sep = ",", header = TRUE)
  whistles_list <- as.data.frame(whistles_list)
  whistles_list$recording <- as.factor(whistles_list$recording)
  whistles_list$whistle_type_chr <- as.character(whistles_list$whistle_type_chr)
  subset(whistles_list, !grepl("Unknown", whistle_type_chr))
}

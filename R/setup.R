# Resolve repository root and source R modules.
get_repo_root <- function() {
  if (exists("REPO_ROOT", envir = .GlobalEnv)) {
    return(get("REPO_ROOT", envir = .GlobalEnv))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

source_repo <- function() {
  repo_root <- get_repo_root()
  assign("REPO_ROOT", repo_root, envir = .GlobalEnv)
  source(file.path(repo_root, "R", "config.R"), local = FALSE)
  source(file.path(repo_root, "R", "load_whistles.R"), local = FALSE)
  source(file.path(repo_root, "R", "markov_whistle.R"), local = FALSE)
  source(file.path(repo_root, "R", "graph_utils.R"), local = FALSE)
  source(file.path(repo_root, "R", "plot_io.R"), local = FALSE)
  invisible(repo_root)
}

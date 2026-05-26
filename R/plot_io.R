# Save figures under outputs/plots (see PLOTS_DIR in config.R).

plots_dir <- function(repo_root) {
  file.path(repo_root, PLOTS_DIR)
}

ensure_plots_dir <- function(repo_root) {
  out <- plots_dir(repo_root)
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  invisible(out)
}

plot_path <- function(name, repo_root, ext = "pdf") {
  file.path(plots_dir(repo_root), paste0(name, ".", ext))
}

with_pdf_plot <- function(path, expr, width = 7, height = 7, ...) {
  pdf(path, width = width, height = height, ...)
  on.exit(dev.off(), add = TRUE)
  force(expr)
  message("Saved plot ", path)
  invisible(path)
}

save_ggplot <- function(plot, name, repo_root, width = 8, height = 6, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to save plots", call. = FALSE)
  }
  path <- plot_path(name, repo_root)
  ggplot2::ggsave(filename = path, plot = plot, width = width, height = height, ...)
  message("Saved plot ", path)
  invisible(path)
}

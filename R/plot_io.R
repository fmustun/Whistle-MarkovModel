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

# Draws expr once (into the pdf device), then replays the recorded plot
# into a png device so base-graphics plots (igraph, etc.) don't need to be
# evaluated twice.
with_pdf_plot <- function(
    path, expr, width = 7, height = 7, png_res = 300, bg = "white", ...) {
  # pdf()'s own default bg is "transparent"; recordPlot() bakes whatever
  # bg was active into the display list, so without this the replayed png
  # below would inherit a transparent background instead of a white one.
  pdf(path, width = width, height = height, bg = bg, ...)
  on.exit(dev.off(), add = TRUE)
  # Non-interactive devices default to displaylist = "inhibit", so
  # recordPlot() below would otherwise capture nothing to replay.
  dev.control(displaylist = "enable")
  force(expr)
  recorded <- recordPlot()
  dev.off()
  on.exit()

  png_path <- sub("\\.pdf$", ".png", path)
  png(png_path, width = width, height = height, units = "in", res = png_res)
  on.exit(dev.off(), add = TRUE)
  replayPlot(recorded)

  message("Saved plot ", path, " and ", png_path)
  invisible(c(path, png_path))
}

save_ggplot <- function(plot, name, repo_root, width = 8, height = 6, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to save plots", call. = FALSE)
  }
  paths <- vapply(c("pdf", "png"), function(ext) {
    path <- plot_path(name, repo_root, ext = ext)
    ggplot2::ggsave(
      filename = path, plot = plot, width = width, height = height, ...
    )
    path
  }, character(1))
  message("Saved plot ", paths[["pdf"]], " and ", paths[["png"]])
  invisible(paths)
}

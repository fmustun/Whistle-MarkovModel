#!/usr/bin/env Rscript
#
# Composites three already-generated figures into one "network analysis"
# summary figure (A4, PNG + PDF): the community-detection network (top,
# with the category legend - outputs/plots/MM_legend.png - composited into
# its blank canvas space), the betweenness-vs-strength z-score scatter
# (middle, outputs/plots/markov_betweenness_vs_strength_zscore.png), and
# the turn-taking panel (bottom - the first of turn_taking_and_centrality
# .png's three stacked ggplot panels).
#
# Pure raster compositing via the system ImageMagick `convert`/`identify`
# binaries (the R `magick` package isn't installed here and needs system
# headers this environment doesn't have). PDF export does NOT go through
# ImageMagick: this system's policy.xml blocks the PDF/PS/EPS/XPS coders
# outright (a common Ghostscript-CVE hardening default - see
# /etc/ImageMagick-6/policy.xml), so the final PNG is instead wrapped into a
# single-page PDF via Python's Pillow (already installed), which writes PDF
# natively without invoking Ghostscript at all.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  REPO_ROOT <- normalizePath(
    file.path(dirname(sub("^--file=", "", file_arg[1])), ".."),
    mustWork = TRUE
  )
} else {
  REPO_ROOT <- normalizePath(getwd(), mustWork = TRUE)
}
setwd(REPO_ROOT)

source(file.path(REPO_ROOT, "R", "config.R"))

if (!nzchar(Sys.which("convert")) || !nzchar(Sys.which("identify"))) {
  stop("ImageMagick's `convert`/`identify` binaries are required but not found on PATH.", call. = FALSE)
}
if (!nzchar(Sys.which("python3"))) {
  stop("python3 (with Pillow) is required for PDF export but not found on PATH.", call. = FALSE)
}

plots_root <- file.path(REPO_ROOT, PLOTS_DIR)
network_measures_dir <- file.path(plots_root, "network_measures")

community_detection_path <- file.path(plots_root, "community_detection.png")
legend_path <- file.path(plots_root, "MM_legend.png")
zscore_path <- file.path(plots_root, "markov_betweenness_vs_strength_zscore.png")
turn_taking_path <- file.path(network_measures_dir, "turn_taking_and_centrality.png")

required <- c(
  community_detection_path = community_detection_path,
  legend_path = legend_path,
  zscore_path = zscore_path,
  turn_taking_path = turn_taking_path
)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop(
    "Missing required input plot(s):\n",
    paste0("  - ", missing, collapse = "\n"),
    "\nRun scripts/03_global_network_measures.R and scripts/02_node_network_measures.R first.",
    call. = FALSE
  )
}

scratch_dir <- tempfile("network_analysis_figure_")
dir.create(scratch_dir)
on.exit(unlink(scratch_dir, recursive = TRUE), add = TRUE)

im_convert <- function(...) {
  convert_args <- c(...)
  status <- system2("convert", convert_args)
  if (!identical(status, 0L)) {
    stop("ImageMagick `convert` failed: convert ", paste(convert_args, collapse = " "), call. = FALSE)
  }
  invisible(NULL)
}

im_size <- function(path) {
  # "%wx%h" (no literal space in the format string), not "%w %h": passing a
  # format string containing a space through system2() somehow gets it
  # split into two argv elements by the time ImageMagick sees them (visible
  # even though the same call works fine typed directly in a shell), so
  # identify receives a bogus second "file" argument. Avoiding the embedded
  # space sidesteps it entirely.
  out <- system2("identify", c("-format", "%wx%h", path), stdout = TRUE)
  as.integer(strsplit(out, "x", fixed = TRUE)[[1]])
}

im_minima <- function(path) {
  out <- system2(
    "convert",
    c(path, "-colorspace", "Gray", "-format", "%[fx:minima]", "info:"),
    stdout = TRUE
  )
  as.numeric(out)
}

trim_image <- function(src, dst, fuzz = "2%") {
  im_convert(src, "-fuzz", fuzz, "-trim", "+repage", dst)
}

# --- Step 1: composite the legend into community_detection's blank space ---
#
# The top-right corner of community_detection.png (x: 1890-2700, y: 0-810 in
# its native 2700x2700 canvas) is reliably blank there - found by a fine
# grid scan of pixel minima before hardcoding this box. Re-verified below
# (rather than trusted blindly) since a materially different Leiden
# partition could shift the layout enough to encroach on it.
legend_box <- list(x = 1890L, y = 0L, w = 810L, h = 810L)

box_crop <- file.path(scratch_dir, "legend_box_check.png")
im_convert(
  community_detection_path,
  "-crop", sprintf("%dx%d+%d+%d", legend_box$w, legend_box$h, legend_box$x, legend_box$y),
  "+repage", box_crop
)
box_minima <- im_minima(box_crop)
if (box_minima < 0.97) {
  stop(
    "The hardcoded legend placement box in community_detection.png is no ",
    "longer blank (min pixel intensity ", round(box_minima, 3), " < 0.97) - ",
    "the community layout has likely changed enough to overlap it. Re-pick ",
    "a blank region (see legend_box in this script) before re-running.",
    call. = FALSE
  )
}

legend_native <- im_size(legend_path)
padding <- 70L
usable <- list(
  x = legend_box$x + padding, y = legend_box$y + padding,
  w = legend_box$w - 2L * padding, h = legend_box$h - 2L * padding
)
legend_scale <- min(usable$w / legend_native[1], usable$h / legend_native[2])
legend_w <- as.integer(round(legend_native[1] * legend_scale))
legend_h <- as.integer(round(legend_native[2] * legend_scale))

legend_resized <- file.path(scratch_dir, "legend_resized.png")
im_convert(legend_path, "-resize", sprintf("%dx%d", legend_w, legend_h), legend_resized)

legend_x <- as.integer(round(usable$x + (usable$w - legend_w) / 2))
legend_y <- as.integer(round(usable$y + (usable$h - legend_h) / 2))

top_panel_raw <- file.path(scratch_dir, "top_panel_raw.png")
im_convert(
  community_detection_path, legend_resized,
  "-geometry", sprintf("+%d+%d", legend_x, legend_y), "-composite",
  top_panel_raw
)
top_panel <- file.path(scratch_dir, "top_panel.png")
trim_image(top_panel_raw, top_panel)
message("Top panel: community detection + legend")

# --- Step 2: middle panel (used as-is, just trimmed of its thin margin) ---
middle_panel <- file.path(scratch_dir, "middle_panel.png")
trim_image(zscore_path, middle_panel)
message("Middle panel: betweenness vs. strength z-score")

# --- Step 3: bottom panel - first (top) of turn_taking_and_centrality's ---
# three equal-height stacked ggplot panels (g1/g3/g4 in
# scripts/02_node_network_measures.R, arranged via gridExtra::arrangeGrob()
# with no explicit `heights`, so each gets exactly 1/3 of the image height).
turn_taking_size <- im_size(turn_taking_path)
panel_h <- as.integer(round(turn_taking_size[2] / 3))
bottom_panel_raw <- file.path(scratch_dir, "bottom_panel_raw.png")
im_convert(
  turn_taking_path,
  "-crop", sprintf("%dx%d+0+0", turn_taking_size[1], panel_h),
  "+repage", bottom_panel_raw
)
bottom_panel <- file.path(scratch_dir, "bottom_panel.png")
trim_image(bottom_panel_raw, bottom_panel)
message("Bottom panel: turn-taking (P(next node in a different category))")

# --- Step 4: lay all three out on an A4 page, common width, stacked -------
DPI <- 300L
PAGE_W <- 2480L # 210mm at 300dpi
PAGE_H <- 3508L # 297mm at 300dpi
margin <- 100L
gap <- 50L
content_w <- PAGE_W - 2L * margin
content_h <- PAGE_H - 2L * margin - 2L * gap

panel_paths <- c(top_panel, middle_panel, bottom_panel)
sizes <- lapply(panel_paths, im_size)
aspect_h_over_w <- vapply(sizes, function(s) s[2] / s[1], numeric(1))

# A common width shared by all three panels (rather than equal height
# bands) keeps their left/right edges aligned - solved so the stacked
# heights, at that common width, exactly fill the page's content height.
common_width <- min(content_h / sum(aspect_h_over_w), content_w)

resize_to_width <- function(src, dst, width) {
  im_convert(src, "-resize", sprintf("%dx", as.integer(round(width))), dst)
}

final_paths <- file.path(scratch_dir, c("panel1_final.png", "panel2_final.png", "panel3_final.png"))
for (i in seq_along(panel_paths)) {
  resize_to_width(panel_paths[i], final_paths[i], common_width)
}
final_sizes <- lapply(final_paths, im_size)
total_h <- sum(vapply(final_sizes, function(s) s[2], numeric(1))) + 2L * gap

# Guard against rounding pushing the total just past the content budget:
# shrink once more rather than letting a panel spill past the page margin.
if (total_h > content_h) {
  common_width <- common_width * (content_h / total_h)
  for (i in seq_along(panel_paths)) {
    resize_to_width(panel_paths[i], final_paths[i], common_width)
  }
  final_sizes <- lapply(final_paths, im_size)
  total_h <- sum(vapply(final_sizes, function(s) s[2], numeric(1))) + 2L * gap
}

y_start <- margin + max(0L, as.integer(round((content_h - total_h) / 2)))
positions <- vector("list", length(final_paths))
y_cursor <- y_start
for (i in seq_along(final_paths)) {
  w_i <- final_sizes[[i]][1]
  h_i <- final_sizes[[i]][2]
  x_i <- margin + as.integer(round((content_w - w_i) / 2))
  positions[[i]] <- list(x = x_i, y = y_cursor, w = w_i, h = h_i)
  y_cursor <- y_cursor + h_i + gap
}

canvas <- file.path(scratch_dir, "canvas.png")
im_convert(
  "-size", sprintf("%dx%d", PAGE_W, PAGE_H),
  "-units", "PixelsPerInch", "-density", as.character(DPI),
  "xc:white", canvas
)

composite_args <- canvas
for (i in seq_along(final_paths)) {
  composite_args <- c(
    composite_args, final_paths[i],
    "-geometry", sprintf("+%d+%d", positions[[i]]$x, positions[[i]]$y),
    "-composite"
  )
}
composited <- file.path(scratch_dir, "composited.png")
im_convert(composite_args, composited)

ensure_plots_dir_exists <- dir.create(plots_root, showWarnings = FALSE, recursive = TRUE)

out_png <- file.path(plots_root, "network_analysis.png")
im_convert(composited, "-units", "PixelsPerInch", "-density", as.character(DPI), out_png)
message("Saved ", out_png)

out_pdf <- file.path(plots_root, "network_analysis.pdf")
# Written to a script file rather than passed inline via `python3 -c <code>`:
# system2() on this system doesn't reliably preserve embedded spaces/
# newlines within a single argv element (confirmed above for identify's
# -format string too) - a multi-line -c argument gets mangled into several
# shell words. A file path has no such whitespace, so passing it as the
# script argument sidesteps the issue entirely.
pil_script <- file.path(scratch_dir, "png_to_pdf.py")
writeLines(
  c(
    "import sys",
    "from PIL import Image",
    "src, dst, dpi = sys.argv[1], sys.argv[2], float(sys.argv[3])",
    "Image.open(src).convert('RGB').save(dst, resolution=dpi)"
  ),
  pil_script
)
pdf_status <- system2("python3", c(pil_script, out_png, out_pdf, as.character(DPI)))
if (!identical(pdf_status, 0L)) {
  stop("Failed to write PDF via Pillow (python3 exit status ", pdf_status, ")", call. = FALSE)
}
message("Saved ", out_pdf)

#!/usr/bin/env Rscript

# Overlays the longest observed multi-loop chains directly on script 01's
# significant Markov transition network (same nodes/layout/edges as
# markov_significant_network.pdf, unchanged - same background scripts
# 06/12 already reuse). Each qualifying chain (chain_length >
# LONG_CHAIN_MIN_LENGTH) is drawn as an open path through its whistles'
# node positions, in temporal order, one color per chain:
#  - a step between two DIFFERENT consecutive sub-categories is an arrow
#    from the earlier node to the later one.
#  - a run of consecutive whistles of the SAME sub-category (a self-repeat,
#    same notion as compute_self_repeat_counts()) is drawn as a small open
#    loop at that node instead of a zero-length segment, labeled "xN" for
#    an N-whistle run - several of the longest chains here are dominated by
#    exactly this (e.g. one 13-whistle chain is 9 consecutive repeats of a
#    single sub-category).
# Paths are open lines/loops (no polygon fill), unlike script 12's filled
# triplet triangles, since a temporal path has no interior to fill.
# A chain with a member sub-category outside the significant network can't
# be placed and is skipped entirely - same handling as scripts 06/12.

suppressPackageStartupMessages({
  library(igraph)
})

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
source(file.path(REPO_ROOT, "R", "load_whistles.R"))
source(file.path(REPO_ROOT, "R", "graph_utils.R"))
source(file.path(REPO_ROOT, "R", "multiloop_utils.R"))
source(file.path(REPO_ROOT, "R", "plot_io.R"))

ensure_plots_dir(REPO_ROOT, "multiloops")
multiloop_dir <- file.path(REPO_ROOT, OUTPUT_DIR, "multiloop_analysis")
dir.create(multiloop_dir, recursive = TRUE, showWarnings = FALSE)

# User-specified: "the longest loops, maybe the ones with more than 7" - not
# a derived statistical threshold, just how many whistles counts as "long"
# here.
LONG_CHAIN_MIN_LENGTH <- 7
LOOP_RADIUS <- 0.045

model_rds <- file.path(REPO_ROOT, OUTPUT_DIR, "markov_model.rds")
if (!file.exists(model_rds)) {
  stop("Run scripts/01_run_markov_model.R first. Missing: ", model_rds, call. = FALSE)
}

saved <- readRDS(model_rds)
gra1 <- saved$gra1

register_fcircle_shape()

# Same function/seed/graph as script 01 uses for markov_significant_network.pdf
# (and scripts 06/12 reuse the same way) - guarantees identical node positions.
coordinates <- compute_markov_layout(gra1, PLOT_SEED)
valid_node_names <- V(gra1)$name

message("Loading whistles from ", DATA_PATH)
whistles_list <- load_whistles(DATA_PATH)

message(
  "Grouping whistles into multi-loop chains (IWI < ",
  MULTILOOP_IWI_THRESHOLD * 1000, "ms)..."
)
whistle_chains <- compute_multiloop_chains(
  whistles_list, threshold = MULTILOOP_IWI_THRESHOLD, list_names = LIST_NAMES
)
long_chains <- whistle_chains[whistle_chains$chain_length > LONG_CHAIN_MIN_LENGTH, ]
long_chains <- long_chains[order(long_chains$chain_id, long_chains$position_in_chain), ]
chain_ids <- unique(long_chains$chain_id)

message(
  length(chain_ids), " multi-loop chain(s) with chain_length > ",
  LONG_CHAIN_MIN_LENGTH
)

# A chain can only be traced on this layout if every whistle type it visits
# is itself a node in the significant network (same requirement scripts
# 06/12 apply to multi-loop edges/pairs/triplets).
in_network_by_chain <- vapply(chain_ids, function(cid) {
  types <- long_chains$whistle_type[long_chains$chain_id == cid]
  all(as.character(types) %in% valid_node_names)
}, logical(1))
n_dropped <- sum(!in_network_by_chain)
chain_ids <- chain_ids[in_network_by_chain]

message(
  length(chain_ids), " / ", length(chain_ids) + n_dropped,
  " qualifying chains plotted (", n_dropped,
  " dropped: a visited sub-category has no significant Markov transition)"
)

chain_colors <- setNames(
  grDevices::hcl.colors(length(chain_ids), palette = "Dark 3"), chain_ids
)

# Run-length-encodes one chain's ordered whistle_type sequence into steps
# (one row per maximal run of the same sub-category), so a stutter of N
# consecutive same-type whistles becomes one "self-repeat of length N" step
# rather than N-1 invisible zero-length segments.
build_chain_steps <- function(cid) {
  sub <- long_chains[long_chains$chain_id == cid, ]
  run <- rle(as.character(sub$whistle_type))
  data.frame(
    whistle_type = run$values,
    run_length = run$lengths,
    stringsAsFactors = FALSE
  )
}
chain_steps <- setNames(lapply(chain_ids, build_chain_steps), chain_ids)

path_table <- do.call(rbind, lapply(chain_ids, function(cid) {
  sub <- long_chains[long_chains$chain_id == cid, ]
  data.frame(
    chain_id = cid,
    recording = sub$recording[1],
    chain_length = sub$chain_length[1],
    position_in_chain = sub$position_in_chain,
    whistle_type = sub$whistle_type,
    whistle_type_chr = sub$whistle_type_chr,
    stringsAsFactors = FALSE
  )
}))
path_csv <- file.path(multiloop_dir, "longest_multiloop_paths.csv")
write.csv(path_table, path_csv, row.names = FALSE)
message("Saved ", path_csv)

# Self-repeat loops are drawn just outside the node, in a direction that
# fans out per chain (so multiple long chains that all stutter on the same
# popular sub-category - several of these do, on SW_Neo_Category_20 - don't
# draw their loops directly on top of one another).
node_outward_angle <- function(pos) {
  if (all(pos == 0)) 0 else atan2(pos["y"], pos["x"])
}

message("Plotting longest multi-loop paths on the significant network...")
with_pdf_plot(
  plot_path("longest_multiloop_paths_overlay", REPO_ROOT, subdir = "multiloops"),
  width = 10,
  height = 10,
  {
    plot_markov_graph(
      gra1,
      list_names = LIST_NAMES,
      seed = PLOT_SEED,
      coordinates = coordinates,
      show_legend = FALSE,
      caption = c(
        "background: script 01's significant Markov transition network (same nodes/layout/edges, unchanged)",
        paste0(
          "colored path = one observed multi-loop chain with chain_length > ",
          LONG_CHAIN_MIN_LENGTH, " (", length(chain_ids), " chains, ", n_dropped,
          " skipped: a visited sub-category has no significant Markov transition)"
        ),
        "arrow = transition to a different sub-category; open loop = N consecutive whistles of the same sub-category (labeled xN)"
      )
    )

    for (cid in chain_ids) {
      steps <- chain_steps[[cid]]
      pts <- coordinates[steps$whistle_type, , drop = FALSE]
      color <- chain_colors[[cid]]
      chain_angle_offset <- (match(cid, chain_ids) - 1) * (2 * pi / 7)

      n_steps <- nrow(steps)
      if (n_steps > 1) {
        arrows(
          pts[-n_steps, "x"], pts[-n_steps, "y"],
          pts[-1, "x"], pts[-1, "y"],
          col = adjustcolor(color, alpha.f = 0.85), lwd = 2.2,
          length = 0.08, angle = 25, code = 2
        )
      }

      repeated <- which(steps$run_length > 1)
      for (i in repeated) {
        base_angle <- node_outward_angle(pts[i, ]) + chain_angle_offset
        loop_center <- pts[i, ] + LOOP_RADIUS * c(cos(base_angle), sin(base_angle))
        symbols(
          loop_center["x"], loop_center["y"], circles = LOOP_RADIUS * 0.6,
          fg = color, bg = NA, lwd = 2, add = TRUE, inches = FALSE
        )
        text(
          loop_center["x"], loop_center["y"] + LOOP_RADIUS * 1.1,
          labels = paste0("x", steps$run_length[i]),
          col = color, cex = 0.6, xpd = NA
        )
      }
    }

    legend(
      "topleft", inset = 0.01, xpd = NA, bty = "n", cex = 0.55,
      legend = paste0(chain_ids, " (n=", sapply(chain_ids, function(cid) {
        long_chains$chain_length[long_chains$chain_id == cid][1]
      }), ")"),
      col = chain_colors, lty = 1, lwd = 2.5
    )
  }
)

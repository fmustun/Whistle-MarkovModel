#!/usr/bin/env Rscript

# Manuscript Fig. S5: a 2x2 panel of script 01's significant Markov
# transition network (same nodes/layout as markov_significant_network.pdf,
# unchanged), drawn four times. Each panel
# highlights several example multi-loop chains (IWI < MULTILOOP_IWI_
# THRESHOLD) at once, each in its own color, as a path through its
# whistles' node positions:
#  - background (non-highlighted) edges are all drawn in a uniform light
#    gray with a small arrowhead, rather than the usual per-edge
#    weight-based shading, so they read as unobtrusive context.
#  - a highlighted chain's transitions are drawn through that same
#    plot.igraph() call - reusing the actual significant-network edge (just
#    recolored) where one already connects the same two nodes in the same
#    direction, or added as a new edge styled the same way (curvature,
#    arrowhead shape) otherwise - so a highlighted step looks exactly like
#    a background one except for color, width, and arrowhead size.
#  - a run of consecutive whistles of the SAME sub-category (a self-repeat,
#    same notion as compute_self_repeat_counts()) has no edge to reuse
#    (gra1 has none - self-transitions are BH-significant-border-only
#    there), so it stays a small open loop drawn next to that node instead,
#    labeled "xN" for an N-whistle run, in the chain's own color.
#
# Within one panel, the example chains (3 or 4 of them, MIN/MAX_EXAMPLES_
# PER_PANEL below) are chosen, in order of preference, to:
#  1. never share a transition (an unordered consecutive-sub-category pair,
#     self-repeats included) with each other, so every color in a panel
#     marks genuinely distinct structure - a hard requirement.
#  2. never have a segment of one chain's path geometrically cross a
#     segment of another chain's path (paths_cross()) - a soft preference,
#     relaxed (chains allowed to cross) only if honoring it would leave a
#     panel with fewer than MIN_EXAMPLES_PER_PANEL chains.
# Chain identities are not reused across panels either (so all four panels
# show different examples), but the same *transition* can appear
# highlighted in more than one panel - only the within-panel constraints
# above are enforced, per request. Picked longest-chain-first (greedy)
# among chains that qualify, so examples stay as visually substantial as
# those constraints allow.
# Each highlighted transition is labeled with its position in that chain's
# temporal order (1, 2, 3, ...), in the chain's own color, so the drawn
# path's direction of travel is legible even where paths bunch up.
# A chain with a member sub-category outside the significant network can't
# be placed and is skipped entirely.

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

N_PANELS <- 4L
# "3 or 4" non-overlapping examples per panel, per request - not a derived
# quantity, just how many distinct colors stay legible against one
# background (hcl.colors' "Dark 3" qualitative palette). The
# selection below always tries for MAX first and only settles for fewer
# (down to MIN) if the non-crossing preference or the candidate pool
# forces it.
MIN_EXAMPLES_PER_PANEL <- 3L
MAX_EXAMPLES_PER_PANEL <- 4L
LOOP_RADIUS <- 0.045
LABEL_OFFSET <- 0.03
BACKGROUND_EDGE_COLOR <- "gray85"
BACKGROUND_EDGE_WIDTH <- 1
HIGHLIGHT_EDGE_WIDTH <- 2.6
# Bigger than plot_markov_graph()'s normal default (0.95) for highlighted
# edges, so a chain's transitions stand out clearly; smaller for the
# light-gray background, so it reads as unobtrusive context rather than
# competing with the highlighted examples.
GREY_ARROW_SIZE <- 0.6
HIGHLIGHT_ARROW_SIZE <- 1.4

model_rds <- file.path(REPO_ROOT, OUTPUT_DIR, "markov_model.rds")
if (!file.exists(model_rds)) {
  stop("Run scripts/01_run_markov_model.R first. Missing: ", model_rds, call. = FALSE)
}

saved <- readRDS(model_rds)
gra1 <- saved$gra1

register_fcircle_shape()

# Same function/seed/graph as script 01 uses for markov_significant_network.pdf
# - guarantees identical node positions across all four panels and Fig. 4.
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
loop_whistles <- whistle_chains[whistle_chains$chain_length >= 2, ]
loop_whistles <- loop_whistles[
  order(loop_whistles$chain_id, loop_whistles$position_in_chain),
]
chain_ids_all <- unique(loop_whistles$chain_id)
chain_length_by_id <- tapply(
  loop_whistles$chain_length, loop_whistles$chain_id, function(x) x[1]
)

# A chain can only be traced on this layout if every whistle type it visits
# is itself a node in the significant network (same requirement scripts
# 06/12/13 apply to multi-loop edges/pairs/triplets/paths).
in_network_by_chain <- vapply(chain_ids_all, function(cid) {
  types <- loop_whistles$whistle_type[loop_whistles$chain_id == cid]
  all(as.character(types) %in% valid_node_names)
}, logical(1))
n_dropped <- sum(!in_network_by_chain)
chain_ids_valid <- chain_ids_all[in_network_by_chain]

message(
  length(chain_ids_valid), " / ", length(chain_ids_all),
  " multi-loop chains qualify (", n_dropped,
  " dropped: a visited sub-category has no significant Markov transition)"
)

# Run-length-encodes one chain's ordered whistle_type sequence into steps
# (one row per maximal run of the same sub-category), so a stutter of N
# consecutive same-type whistles becomes one "self-repeat of length N" step
# rather than N-1 invisible zero-length segments.
build_chain_steps <- function(cid) {
  sub <- loop_whistles[loop_whistles$chain_id == cid, ]
  run <- rle(as.character(sub$whistle_type))
  data.frame(
    whistle_type = run$values, run_length = run$lengths,
    stringsAsFactors = FALSE
  )
}

# The set of transitions a chain's path would draw: an unordered
# consecutive-sub-category pair for every step-to-step arrow, plus a
# same-type self pair for every self-repeat run. Unordered (pmin/pmax)
# because two chains that traverse the same node pair in opposite
# directions would still draw over the same screen segment.
chain_edge_keys <- function(steps) {
  keys <- character(0)
  self_repeat <- steps$run_length > 1
  if (any(self_repeat)) {
    keys <- c(keys, paste0(
      steps$whistle_type[self_repeat], "|", steps$whistle_type[self_repeat]
    ))
  }
  n_steps <- nrow(steps)
  if (n_steps > 1) {
    a <- steps$whistle_type[-n_steps]
    b <- steps$whistle_type[-1]
    keys <- c(keys, paste0(pmin(a, b), "|", pmax(a, b)))
  }
  unique(keys)
}

# Unique directed (from, to) whistle_type pairs a chain's arrows would
# connect, in first-encountered order - the actual edges added/recolored
# on the plotting graph (unlike chain_edge_keys()'s unordered pairs, which
# only exist to check for overlap).
chain_directed_edges <- function(steps) {
  n_steps <- nrow(steps)
  if (n_steps <= 1) {
    return(data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE))
  }
  pairs <- data.frame(
    from = steps$whistle_type[-n_steps], to = steps$whistle_type[-1],
    stringsAsFactors = FALSE
  )
  pairs[!duplicated(pairs), ]
}

# Every step-to-step transition a chain visits, in temporal order,
# *without* deduplication (unlike chain_directed_edges()) - used only to
# number each drawn arrow with its position in the chain (1, 2, 3, ...),
# since a chain that revisits the same directed pair should still get a
# fresh number at each occurrence.
chain_ordered_transitions <- function(steps) {
  n_steps <- nrow(steps)
  if (n_steps <= 1) {
    return(data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE))
  }
  data.frame(
    from = steps$whistle_type[-n_steps], to = steps$whistle_type[-1],
    stringsAsFactors = FALSE
  )
}

# Proper segment intersection test (do open segments p1-p2 and p3-p4 cross
# in their interior?), via the standard orientation/sign test. Two segments
# that only touch at a shared endpoint are not treated as crossing by
# paths_cross() below (that's normal - two chains passing through the same
# hub node - not a crossing worth avoiding).
segments_intersect <- function(p1, p2, p3, p4) {
  orient <- function(a, b, c) {
    sign((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]))
  }
  d1 <- orient(p3, p4, p1)
  d2 <- orient(p3, p4, p2)
  d3 <- orient(p1, p2, p3)
  d4 <- orient(p1, p2, p4)
  (d1 != d2) && (d3 != d4) && d1 != 0 && d2 != 0 && d3 != 0 && d4 != 0
}

# Does any segment of chain A's path cross any segment of chain B's path?
# Segment pairs that share a node (same from/to id on both sides) are
# skipped - meeting at a common node is not a "crossing" to avoid, only
# genuinely separate segments passing through each other in open space are.
paths_cross <- function(edges_a, edges_b) {
  if (nrow(edges_a) == 0 || nrow(edges_b) == 0) return(FALSE)
  for (i in seq_len(nrow(edges_a))) {
    a_from <- edges_a$from[i]
    a_to <- edges_a$to[i]
    p1 <- coordinates[a_from, ]
    p2 <- coordinates[a_to, ]
    for (j in seq_len(nrow(edges_b))) {
      b_from <- edges_b$from[j]
      b_to <- edges_b$to[j]
      shared_node <- a_from == b_from || a_from == b_to ||
        a_to == b_from || a_to == b_to
      if (shared_node) next
      if (segments_intersect(p1, p2, coordinates[b_from, ], coordinates[b_to, ])) {
        return(TRUE)
      }
    }
  }
  FALSE
}

# Greedy, longest-chain-first: at each step, take the next-longest
# remaining candidate (excluding chain ids already used in an earlier
# panel) whose transitions don't intersect any already-chosen chain's
# transitions *within this panel*. Deterministic (no randomness/seed
# needed) - ties in chain_length broken by chain_id for a stable order.
candidate_order <- chain_ids_valid[
  order(-chain_length_by_id[chain_ids_valid], chain_ids_valid)
]

select_panel_chains <- function(exclude_chain_ids) {
  used_edges <- character(0)
  chosen <- character(0)
  chosen_edges <- list()

  # Pass 1: require both no shared transition AND no path crossing.
  for (cid in candidate_order) {
    if (cid %in% exclude_chain_ids || cid %in% chosen) next
    if (length(chosen) >= MAX_EXAMPLES_PER_PANEL) break
    steps <- build_chain_steps(cid)
    edges <- chain_edge_keys(steps)
    if (length(intersect(edges, used_edges)) > 0L) next
    directed <- chain_directed_edges(steps)
    crosses <- any(vapply(
      chosen_edges, function(e) paths_cross(directed, e), logical(1)
    ))
    if (crosses) next
    chosen <- c(chosen, cid)
    chosen_edges[[cid]] <- directed
    used_edges <- union(used_edges, edges)
  }

  # Pass 2 (only if pass 1 fell short of the minimum): fill remaining slots
  # requiring only no shared transition, crossings now allowed.
  if (length(chosen) < MIN_EXAMPLES_PER_PANEL) {
    for (cid in candidate_order) {
      if (cid %in% exclude_chain_ids || cid %in% chosen) next
      if (length(chosen) >= MAX_EXAMPLES_PER_PANEL) break
      edges <- chain_edge_keys(build_chain_steps(cid))
      if (length(intersect(edges, used_edges)) == 0L) {
        chosen <- c(chosen, cid)
        used_edges <- union(used_edges, edges)
      }
    }
  }
  chosen
}

global_used_chain_ids <- character(0)
panel_chain_ids <- vector("list", N_PANELS)
for (p in seq_len(N_PANELS)) {
  chosen <- select_panel_chains(global_used_chain_ids)
  panel_chain_ids[[p]] <- chosen
  global_used_chain_ids <- union(global_used_chain_ids, chosen)
  message(
    "Panel ", p, ": ", length(chosen),
    " / ", MIN_EXAMPLES_PER_PANEL, "-", MAX_EXAMPLES_PER_PANEL,
    " non-overlapping example chains selected"
  )
}

path_table <- do.call(rbind, lapply(seq_len(N_PANELS), function(p) {
  chain_ids <- panel_chain_ids[[p]]
  if (length(chain_ids) == 0) return(NULL)
  do.call(rbind, lapply(seq_along(chain_ids), function(k) {
    cid <- chain_ids[k]
    sub <- loop_whistles[loop_whistles$chain_id == cid, ]
    data.frame(
      panel = p,
      example_in_panel = k,
      chain_id = cid,
      recording = sub$recording[1],
      chain_length = sub$chain_length[1],
      position_in_chain = sub$position_in_chain,
      whistle_type = sub$whistle_type,
      whistle_type_chr = sub$whistle_type_chr,
      stringsAsFactors = FALSE
    )
  }))
}))
path_csv <- file.path(multiloop_dir, "loop_examples_panel_paths.csv")
write.csv(path_table, path_csv, row.names = FALSE)
message("Saved ", path_csv)

# Self-repeat loops are drawn just outside the node, in the direction
# pointing away from the plot origin. No per-chain fan-out offset is needed:
# two chains in the same panel can never share a self-repeat at the same
# node, since that would be a shared "type|type" edge key and the selection
# above already forbids it.
node_outward_angle <- function(pos) {
  if (all(pos == 0)) 0 else atan2(pos["y"], pos["x"])
}

draw_self_repeat_markers <- function(steps, color) {
  pts <- coordinates[steps$whistle_type, , drop = FALSE]
  repeated <- which(steps$run_length > 1)
  for (i in repeated) {
    base_angle <- node_outward_angle(pts[i, ])
    loop_center <- pts[i, ] + LOOP_RADIUS * c(cos(base_angle), sin(base_angle))
    symbols(
      loop_center["x"], loop_center["y"], circles = LOOP_RADIUS * 0.6,
      fg = color, bg = NA, lwd = 2.2, add = TRUE, inches = FALSE
    )
    text(
      loop_center["x"], loop_center["y"] + LOOP_RADIUS * 1.1,
      labels = paste0("x", steps$run_length[i]),
      col = color, cex = 0.55, xpd = NA
    )
  }
}

# Labels each transition arrow with its position in the chain's temporal
# order (1, 2, 3, ...), offset perpendicular to the segment (rather than
# sitting on top of the line) so the number stays legible - offset grows a
# little further out each time the same directed pair recurs later in the
# same chain, so repeated back-and-forth steps don't stack their numbers
# exactly on top of one another.
draw_transition_order_labels <- function(steps, color) {
  n_steps <- nrow(steps)
  if (n_steps <= 1) return(invisible(NULL))
  pts <- coordinates[steps$whistle_type, , drop = FALSE]
  occurrence_count <- list()
  for (i in seq_len(n_steps - 1)) {
    p1 <- pts[i, ]
    p2 <- pts[i + 1, ]
    seg <- p2 - p1
    seg_len <- sqrt(sum(seg^2))
    if (seg_len == 0) next
    perp <- c(-seg["y"], seg["x"]) / seg_len
    pair_key <- paste(steps$whistle_type[i], steps$whistle_type[i + 1])
    prior <- occurrence_count[[pair_key]]
    if (is.null(prior)) prior <- 0L
    occurrence_count[[pair_key]] <- prior + 1L
    offset <- LABEL_OFFSET + prior * LABEL_OFFSET * 0.9
    label_pos <- (p1 + p2) / 2 + offset * perp
    text(
      label_pos["x"], label_pos["y"], labels = as.character(i),
      col = color, cex = 0.95, font = 2, xpd = NA
    )
  }
}

# Builds one panel's plotting graph: gra1 with every edge defaulted to the
# light-gray background style, then each chosen chain's directed
# transitions layered on top in its own color - reusing the matching gra1
# edge (recolored) when one already connects those two nodes in that
# direction, or adding a new edge styled the same way otherwise, so
# highlighted and background edges are rendered through the exact same
# plot.igraph() call (identical curvature/arrowhead shape, only color/
# width/arrow-size differ).
build_panel_graph <- function(chain_ids, chain_colors) {
  gra_panel <- gra1
  E(gra_panel)$plot_color <- BACKGROUND_EDGE_COLOR
  E(gra_panel)$plot_width <- BACKGROUND_EDGE_WIDTH
  E(gra_panel)$plot_arrow_size <- GREY_ARROW_SIZE

  for (cid in chain_ids) {
    color <- chain_colors[[cid]]
    directed_edges <- chain_directed_edges(build_chain_steps(cid))
    for (r in seq_len(nrow(directed_edges))) {
      from_id <- directed_edges$from[r]
      to_id <- directed_edges$to[r]
      eid <- get.edge.ids(gra_panel, c(from_id, to_id), directed = TRUE)
      if (eid == 0) {
        gra_panel <- add_edges(gra_panel, c(from_id, to_id))
        eid <- ecount(gra_panel)
      }
      E(gra_panel)$plot_color[eid] <- color
      E(gra_panel)$plot_width[eid] <- HIGHLIGHT_EDGE_WIDTH
      E(gra_panel)$plot_arrow_size[eid] <- HIGHLIGHT_ARROW_SIZE
    }
  }
  gra_panel
}

message("Plotting 2x2 panel of loop examples...")
with_pdf_plot(
  plot_path("loop_examples_panel", REPO_ROOT, subdir = "multiloops"),
  width = 12,
  height = 12,
  {
    layout(matrix(seq_len(N_PANELS), nrow = 2, byrow = TRUE))

    for (p in seq_len(N_PANELS)) {
      chain_ids <- panel_chain_ids[[p]]
      if (length(chain_ids) == 0) {
        par(mar = c(0.5, 0.5, 0.5, 0.5))
        plot.new()
        text(0.5, 0.5, "no non-overlapping\nexample chain available", cex = 0.9)
        next
      }

      chain_colors <- setNames(
        grDevices::hcl.colors(length(chain_ids), palette = "Dark 3"), chain_ids
      )
      gra_panel <- build_panel_graph(chain_ids, chain_colors)

      plot_markov_graph(
        gra_panel,
        list_names = LIST_NAMES,
        seed = PLOT_SEED,
        coordinates = coordinates,
        show_legend = FALSE,
        edge_color = E(gra_panel)$plot_color,
        edge_width = E(gra_panel)$plot_width,
        edge_arrow_size = E(gra_panel)$plot_arrow_size,
        edge_arrow_width = E(gra_panel)$plot_arrow_size
      )

      for (cid in chain_ids) {
        steps <- build_chain_steps(cid)
        color <- chain_colors[[cid]]
        draw_self_repeat_markers(steps, color)
        draw_transition_order_labels(steps, color)
      }

      legend(
        "topleft", inset = 0.01, xpd = NA, bty = "n", cex = 0.58,
        legend = paste0(
          "multi-loop ", seq_along(chain_ids), " (n=",
          vapply(chain_ids, function(cid) chain_length_by_id[[cid]], numeric(1)), ")"
        ),
        col = chain_colors, lty = 1, lwd = 2.5
      )
    }
  }
)

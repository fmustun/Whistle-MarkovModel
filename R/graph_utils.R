build_seq_category <- function(whistles_list, list_names = LIST_NAMES) {
  seq_category <- list()
  for (wh_name in list_names) {
    seq_category[[wh_name]] <- unique(whistles_list$whistle_type)[
      grepl(wh_name, unique(whistles_list$whistle_type_chr))
    ]
  }
  seq_category
}

assign_vertex_category_color <- function(Gra, seq_category, list_names = LIST_NAMES, list_colors = LIST_COLORS) {
  V(Gra)$color <- ifelse(
    V(Gra)$name %in% seq_category[["SW_Neo"]], list_colors[1],
    ifelse(
      V(Gra)$name %in% seq_category[["SW_Luna"]], list_colors[2],
      ifelse(
        V(Gra)$name %in% seq_category[["SW_Yosefa"]], list_colors[3],
        ifelse(
          V(Gra)$name %in% seq_category[["SW_Nikita"]], list_colors[4],
          ifelse(
            V(Gra)$name %in% seq_category[["SW_Nana"]], list_colors[5],
            ifelse(
              V(Gra)$name %in% seq_category[["SW_Dana"]], list_colors[6],
              ifelse(
                V(Gra)$name %in% seq_category[["SW_Shy"]], list_colors[7],
                ifelse(
                  V(Gra)$name %in% seq_category[["NSW_3"]], list_colors[8],
                  ifelse(
                    V(Gra)$name %in% seq_category[["NSW_6"]], list_colors[9],
                    ifelse(
                      V(Gra)$name %in% seq_category[["NSW_9"]], list_colors[10],
                      list_colors[11]
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )

  V(Gra)$category <- ifelse(
    V(Gra)$name %in% seq_category[["SW_Neo"]], list_names[1],
    ifelse(
      V(Gra)$name %in% seq_category[["SW_Luna"]], list_names[2],
      ifelse(
        V(Gra)$name %in% seq_category[["SW_Yosefa"]], list_names[3],
        ifelse(
          V(Gra)$name %in% seq_category[["SW_Nikita"]], list_names[4],
          ifelse(
            V(Gra)$name %in% seq_category[["SW_Nana"]], list_names[5],
            ifelse(
              V(Gra)$name %in% seq_category[["SW_Dana"]], list_names[6],
              ifelse(
                V(Gra)$name %in% seq_category[["SW_Shy"]], list_names[7],
                ifelse(
                  V(Gra)$name %in% seq_category[["NSW_3"]], list_names[8],
                  ifelse(
                    V(Gra)$name %in% seq_category[["NSW_6"]], list_names[9],
                    ifelse(
                      V(Gra)$name %in% seq_category[["NSW_9"]], list_names[10],
                      list_names[11]
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )
  Gra
}

compute_graph <- function(
    transition_matrix,
    Markov_Model,
    Whistles_List,
    sub_division = TRUE,
    list_names = LIST_NAMES,
    list_colors = LIST_COLORS
) {
  Gra <- graph_from_adjacency_matrix(
    transition_matrix, mode = "directed", weighted = TRUE
  )
  vertex_ids <- as.character(seq_len(vcount(Gra)))
  V(Gra)$name <- vertex_ids
  whistle_chr_map <- tapply(
    Whistles_List$whistle_type_chr,
    Whistles_List$whistle_type,
    function(x) x[1]
  )
  V(Gra)$sub_category <- whistle_chr_map[vertex_ids]
  V(Gra)$occurrences <- as.numeric(Markov_Model$whistle_occurrence[vertex_ids])
  Gra <- delete.vertices(Gra, degree(Gra) == 0)

  if (sub_division) {
    seq_category <- build_seq_category(Whistles_List, list_names)
    Gra <- assign_vertex_category_color(Gra, seq_category, list_names, list_colors)
  } else {
    V(Gra)$color <- list_colors
    V(Gra)$category <- list_names
  }

  Gra
}

compute_graph_from_selection <- function(
    transition_matrix,
    selection_values,
    Markov_Model,
    Whistles_List,
    list_names = LIST_NAMES,
    list_colors = LIST_COLORS,
    edge_attribute_name = "adjusted_p",
    raw_p_values = NULL
) {
  graph_with_loops <- compute_graph(
    transition_matrix,
    Markov_Model = Markov_Model,
    Whistles_List = Whistles_List,
    sub_division = TRUE,
    list_names = list_names,
    list_colors = list_colors
  )
  if (ecount(graph_with_loops) > 0L) {
    endpoints <- ends(graph_with_loops, E(graph_with_loops), names = TRUE)
    index <- cbind(as.integer(endpoints[, 1L]), as.integer(endpoints[, 2L]))
    graph_with_loops <- set_edge_attr(
      graph_with_loops, edge_attribute_name, value = selection_values[index]
    )
    if (!is.null(raw_p_values)) {
      graph_with_loops <- set_edge_attr(
        graph_with_loops, "empirical_p", value = raw_p_values[index]
      )
    }
  }
  loop_ids <- which(diag(transition_matrix) > 0)
  frame_color <- ifelse(
    as.integer(V(graph_with_loops)$name) %in% loop_ids, "black", "gray"
  )
  graph_with_loops <- set_vertex_attr(
    graph_with_loops, "vertex.frame.col", value = frame_color
  )
  graph_no_loops <- delete_edges(
    graph_with_loops, which(which_loop(graph_with_loops))
  )
  list(
    with_loops = graph_with_loops,
    no_loops = graph_no_loops,
    significant_loop_ids = loop_ids
  )
}

# Runs cluster_leiden(gra, objective_function = "modularity", ...) once per
# value in `resolutions`, returning one row per resolution with the
# resulting community count and modularity. cluster_leiden() - unlike
# cluster_louvain() - does not store modularity on its returned communities
# object (modularity(cl) errors with "Modularity was not calculated"), so it
# is recomputed explicitly from the membership vector at that same
# resolution_parameter. A fixed seed is set before every call (rather than
# once before the loop) so results are reproducible regardless of how many
# resolutions are swept or in what order, and so real-vs-null comparisons at
# the same resolution aren't confounded by different RNG draws.
leiden_resolution_sweep <- function(
    gra, resolutions, weights = NULL, n_iterations = -1, seed = NULL
) {
  if (is.null(weights) && "weight" %in% edge_attr_names(gra)) {
    weights <- E(gra)$weight
  }
  n_communities <- rep(NA_integer_, length(resolutions))
  modularity_values <- rep(NA_real_, length(resolutions))
  for (i in seq_along(resolutions)) {
    if (!is.null(seed)) set.seed(seed)
    cl <- cluster_leiden(
      gra,
      objective_function = "modularity",
      weights = weights,
      resolution_parameter = resolutions[i],
      n_iterations = n_iterations
    )
    n_communities[i] <- length(cl)
    modularity_values[i] <- modularity(
      gra, membership(cl), weights = weights, resolution = resolutions[i]
    )
  }
  data.frame(
    resolution = resolutions,
    n_communities = n_communities,
    modularity = modularity_values
  )
}

# Rebuilds `gra_template`'s edge weights from `null_counts_matrix` instead of
# the observed data, holding the edge set (topology) fixed - only the edges
# already in gra_template are reweighted, using that same (i, j) pair's
# transition probability under the null. Used to benchmark community
# detection against the shuffle-null count matrices already generated for
# the edge-significance test (MarkovWhistle's null_model = "shuffle",
# R/markov_whistle.R): holding topology fixed isolates whether the observed
# *weight pattern* on the significant edges - not just which edges are
# significant - drives the community structure, and keeps null realizations
# comparable to gra_template in size/density instead of confounding the
# comparison with a much denser unfiltered null network. Requires CalcTransi
# tionMatrix (ccber) to already be on the search path, same convention as
# R/multiloop_utils.R's own use of it.
reweight_graph_from_null_counts <- function(gra_template, null_counts_matrix) {
  null_probability <- CalcTransitionMatrix(null_counts_matrix)
  endpoints <- ends(gra_template, E(gra_template), names = TRUE)
  index <- cbind(as.integer(endpoints[, 1L]), as.integer(endpoints[, 2L]))
  set_edge_attr(gra_template, "weight", value = null_probability[index])
}

# Per-group "compactness": each group's mean within-group transition
# probability divided by its mean between-group transition probability,
# where "between" is strictly outgoing - edges from a node in this group to
# a node in any other group. Edges coming IN from other groups are excluded
# from the denominator on purpose: this is meant to measure the group's own
# tendency to stay inside itself once a whistle is in it, not how attractive
# a destination it is for other groups (those two can differ sharply - e.g.
# a group whose own out-transitions are highly self-loyal but which other
# groups frequently transition into would be pulled toward "not compact" if
# incoming edges counted toward its between-mean, even though nothing about
# its own behavior is diffuse). Above 1 means a node in the group is more
# likely, on average, to transition to another node in the same group than
# to leave it; below 1 means the reverse. `group` is
# deliberately generic - the grouping doesn't have to be a detected
# community: scripts/03_community_detection.R uses this both for Leiden's
# communities and for the main whistle category (V(gra1)$category), and
# those two groupings need not agree on how compact any given node's group
# is.
#
# Takes `gra`'s edges as-is (each weight used directly, no aggregation
# across direction) - pass the *directed* significant-transition graph
# (e.g. gra1), not the undirected/summed graph a caller may have clustered
# on, since collapsing a reciprocal pair into one summed weight (as
# as.undirected()'s default "collapse" does) would no longer be a
# probability in [0, 1]. `membership` is still whatever grouping was used
# on that undirected collapse (or any other grouping entirely) - it's just
# a per-vertex label, independent of which graph's edges are then
# classified with it.
#
# Whether self-loops count toward a group's "within" side depends entirely
# on whether `gra` has them: a self-loop is always within its own node's
# group by construction (it can never be a between-group edge), so passing
# a self-loops-included graph (e.g. bh_analysis/bh_corrected_network.rds's
# graph_with_loops, rather than gra1, which has them stripped out) raises
# every group's within side, never its between side.
#
# edge_counts optionally reweights the mean by each edge's underlying
# observation count (E(gra)$weight's numerator, e.g. MarkovModel$
# transition_counts_matrix indexed the same way reweight_graph_from_null_
# counts() is) - same length/order as E(gra). NULL (default) is a plain
# unweighted mean, treating every edge as one observation regardless of how
# many times that specific transition actually occurred; supplying counts
# instead weights a transition's contribution by how often it was actually
# observed, so a probability backed by many observations moves the mean
# more than one backed by few.
#
# A group with no internal edges (e.g. a singleton) gets
# mean_within_probability = NA and compactness = NA; one with no outgoing
# edges to any other group (everything reaching it is incoming-only, or it's
# a fully isolated component) gets mean_between_probability = NA and
# compactness = NA - both left as NA (not silently dropped or zeroed) since
# "compactness relative to the rest of the network" is genuinely undefined
# without something to compare against.
community_compactness <- function(gra, membership, edge_counts = NULL) {
  if (length(membership) != vcount(gra)) {
    stop("membership must have one entry per vertex of gra", call. = FALSE)
  }
  if (!is.null(edge_counts) && length(edge_counts) != ecount(gra)) {
    stop("edge_counts must have one entry per edge of gra", call. = FALSE)
  }
  endpoints <- ends(gra, E(gra), names = FALSE)
  from_group <- membership[endpoints[, 1L]]
  to_group <- membership[endpoints[, 2L]]
  within_edge <- from_group == to_group
  weight <- E(gra)$weight

  weighted_mean <- function(w, n) {
    if (length(w) == 0L) return(NA_real_)
    if (is.null(edge_counts)) mean(w) else sum(w * n) / sum(n)
  }

  groups <- sort(unique(membership))
  rows <- lapply(groups, function(group_id) {
    within_idx <- within_edge & from_group == group_id
    between_idx <- !within_edge & from_group == group_id
    mean_within <- weighted_mean(weight[within_idx], edge_counts[within_idx])
    mean_between <- weighted_mean(weight[between_idx], edge_counts[between_idx])
    data.frame(
      group = group_id,
      n_nodes = sum(membership == group_id),
      n_within_edges = sum(within_idx),
      mean_within_probability = mean_within,
      n_between_edges = sum(between_idx),
      mean_between_probability = mean_between,
      compactness = mean_within / mean_between
    )
  })
  do.call(rbind, rows)
}

mycircle <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1 / 140 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }

  mapply(
    coords[, 1], coords[, 2], vertex.color, vertex.frame.color,
    vertex.size, vertex.frame.width,
    FUN = function(x, y, bg, fg, size, lwd) {
      symbols(
        x = x, y = y, bg = bg, fg = fg, lwd = lwd,
        circles = size, add = TRUE, inches = FALSE
      )
    }
  )
}

register_fcircle_shape <- function() {
  add.vertex.shape(
    "fcircle",
    clip = igraph.shape.noclip,
    plot = mycircle,
    parameters = list(vertex.frame.color = 1, vertex.frame.width = 1)
  )
}

compute_markov_layout <- function(gra, seed = PLOT_SEED, half_width = 1.3) {
  set.seed(seed)
  # layout_with_fr uses the edge `weight` attribute (raw transition
  # probability, ~0.008-1 here) as attraction strength by default, so
  # near-1 edges pull their nodes almost on top of each other. Raising it
  # to a small power flattens that range so strongly-linked nodes still
  # end up closer together without collapsing into overlapping points.
  # (A near-zero exponent, i.e. ~ignoring weight entirely, was tried too:
  # combined with the wider box below it reads as a tangled hairball,
  # since node position then mostly reflects raw topology rather than
  # transition strength - 0.25 keeps the clusters visually coherent.)
  layout_weights <- E(gra)$weight^0.25
  coordinates <- layout_with_fr(
    gra, niter = 2000, grid = "nogrid", weights = layout_weights
  )
  #coordinates <- layout_with_sugiyama(gra)
  # igraph 2.x dropped layout_with_fr's area/repulserad knobs (they're
  # silently no-ops now), so there's no direct repulsion dial anymore.
  # vertex circles are drawn at a fixed absolute size (see mycircle()),
  # so normalizing into a wider box - while plot_markov_graph derives
  # xlim/ylim from this same range - grows every node's gap relative to
  # that fixed circle size. That's the practical equivalent of turning up
  # repulsion: raising half_width from the original 0.95 pushes every
  # node further from its neighbors relative to its (fixed-size) circle.
  coordinates <- norm_coords(
    coordinates,
    xmin = -half_width, xmax = half_width,
    ymin = -half_width, ymax = half_width
  )
  rownames(coordinates) <- V(gra)$name
  colnames(coordinates) <- c("x", "y")
  coordinates
}

layout_for_graph <- function(coordinates, gra) {
  if (is.null(rownames(coordinates))) {
    stop("Layout coordinates must have vertex IDs as row names", call. = FALSE)
  }
  missing <- setdiff(V(gra)$name, rownames(coordinates))
  if (length(missing)) {
    stop("Layout lacks graph vertices: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  coordinates[V(gra)$name, , drop = FALSE]
}

rescale_from_reference <- function(values, reference_range, output_range) {
  if (length(reference_range) != 2L || !all(is.finite(reference_range)) ||
      diff(reference_range) <= 0) {
    stop("reference_range must contain two finite increasing values",
         call. = FALSE)
  }
  output_range[1L] +
    (values - reference_range[1L]) / diff(reference_range) * diff(output_range)
}

plot_markov_graph <- function(
    gra,
    list_names = LIST_NAMES,
    list_colors = LIST_COLORS,
    seed = PLOT_SEED,
    coordinates = NULL,
    weight_reference_range = NULL,
    show_legend = TRUE,
    main = NULL,
    # Per-vertex label (the number drawn inside each node) color/font
    # override, e.g. so a caller mixing more than one kind of node (script
    # 17's whistle vs. multi-loop-cluster nodes) can make one kind visually
    # distinct from the other. NULL (default) omits both from the
    # plot.igraph() call entirely, leaving every existing caller's look
    # (igraph's own default label color) exactly as it was.
    vertex_label_color = NULL,
    vertex_label_font = NULL,
    # Per-vertex label TEXT override (distinct from vertex_label_color/_font
    # above, which only restyle whatever plot.igraph's own default label -
    # V(gra)$name - already is). NULL (default) reproduces that default
    # exactly. A caller mixing more than one node universe with different
    # numbering schemes (e.g. scripts/10's whistle sub-categories, which
    # keep their whistle_type id from the original per-whistle-type network,
    # alongside multi-loop cluster nodes numbered to continue right after)
    # can pass the display numbers here without renumbering V(gra)$name
    # itself, which every other part of that caller's script still relies on
    # as a dense 1..n_nodes index into its own matrices.
    vertex_label = NULL,
    # Per-vertex circle size / border width override (same NULL-preserves-
    # prior-behavior convention as vertex_label_color/_font above) - e.g. so
    # scripts/10 can draw its multi-loop-cluster nodes larger and
    # thicker-bordered than the standalone-whistle nodes they share a plot
    # with, without changing vertex.frame.color's existing meaning
    # (significant self-transition). NULL reproduces the original constants
    # (vertex.size = 5, vertex.frame.width = 2.5) for every other caller.
    vertex_size = NULL,
    vertex_frame_width = NULL,
    # Edge color/width override (same NULL-preserves-prior-behavior
    # convention as the vertex_* overrides above) - e.g. so a caller can
    # render every edge as a uniform light gray background (a scalar
    # recycles across all edges) instead of the default per-edge alpha/width
    # derived from transition weight, to make a separately-drawn highlighted
    # subset of edges stand out against it.
    edge_color = NULL,
    edge_width = NULL,
    # Arrow head size/width override (same NULL-preserves-prior-behavior
    # convention as edge_color/edge_width above) - a scalar or per-edge
    # vector, e.g. so a caller mixing a light-gray background edge set with
    # a separately colored highlighted subset (both drawn through this same
    # plot.igraph() call, so they share identical curvature/arrow styling)
    # can shrink just the background arrowheads. NULL reproduces the
    # original constants (0.95/0.95) for every other caller.
    edge_arrow_size = NULL,
    edge_arrow_width = NULL,
    # Extra caption lines drawn under the network (e.g. explaining a
    # non-default vertex/edge encoding) via mtext(). Rendered inside this
    # function - not after it returns - so it shares the same par(mar=...)
    # the network itself was drawn with; doing this from the caller after
    # plot_markov_graph() returns would use the restored (larger) default
    # margins instead and misalign the text against the already-rendered plot.
    caption = NULL
) {
  if (is.null(coordinates)) coordinates <- compute_markov_layout(gra, seed)
  coordinates <- layout_for_graph(coordinates, gra)
  # plot.igraph does not touch par("mar"), so the default axis/title
  # margins (5.1/4.1/4.1/2.1 lines) otherwise eat a large, unused border
  # around the network. Shrink to near-zero (more at the bottom when the
  # legend or a caption needs room) so the layout fills the figure.
  # Every existing caller uses either show_legend or caption, never both, so
  # the two single-element branches below are unchanged from before; the
  # combined branch is additive room for the legend row on top of the
  # caption's own lines, for a caller that wants both at once.
  bottom_margin <- if (show_legend && !is.null(caption)) {
    3 + length(caption) * 1.3
  } else if (show_legend) {
    3
  } else if (!is.null(caption)) {
    length(caption) * 1.3 + 0.5
  } else {
    0.5
  }
  old_par <- par(mar = c(bottom_margin, 0.5, 0.5, 0.5))
  on.exit(par(old_par), add = TRUE)
  if (is.null(edge_color)) {
    if (is.null(weight_reference_range)) {
      weight_reference_range <- range(E(gra)$weight)
    }
    edge_alpha <- rescale_from_reference(
      E(gra)$weight, weight_reference_range, c(0.1, 1)
    )
    edge_alpha <- pmin(1, pmax(0.01, edge_alpha))
    edge_color <- rgb(140 / 255, 140 / 255, 140 / 255, edge_alpha)
  }
  if (is.null(edge_width)) {
    if (is.null(weight_reference_range)) {
      weight_reference_range <- range(E(gra)$weight)
    }
    edge_width <- rescale_from_reference(
      E(gra)$weight, weight_reference_range, c(3, 7)
    )
    edge_width <- pmax(0, edge_width)
  }

  # Derived from the coordinates actually passed in (rather than a
  # hardcoded range) so this can never drift out of sync with whatever
  # box compute_markov_layout() (or a caller-supplied layout) used - a
  # mismatch here would clip nodes near the edge out of the plot.
  plot_lim <- range(coordinates) * 1.05

  plot.igraph(
    gra,
    layout = coordinates,
    rescale = FALSE,
    edge.curved = 0.2,
    vertex.size = if (is.null(vertex_size)) 5 else vertex_size,
    edge.color = edge_color,
    edge.width = edge_width,
    edge.arrow.size = if (is.null(edge_arrow_size)) 0.95 else edge_arrow_size,
    edge.arrow.width = if (is.null(edge_arrow_width)) 0.95 else edge_arrow_width,
    vertex.shape = "fcircle",
    vertex.frame.color = V(gra)$vertex.frame.col,
    vertex.frame.width = if (is.null(vertex_frame_width)) 2.5 else vertex_frame_width,
    vertex.label = if (is.null(vertex_label)) V(gra)$name else vertex_label,
    vertex.label.color = vertex_label_color,
    vertex.label.font = vertex_label_font,
    ylim = plot_lim,
    xlim = plot_lim,
    asp = 1,
    main = main
  )
  if (show_legend) {
    legend(
      "bottom", inset = c(0, -0.08), xpd = NA, ncol = 5,
      legend = list_names, pch = 21, pt.bg = list_colors,
      col = "gray", pt.cex = 1.05, cex = 0.62, bty = "n",
      x.intersp = 0.45, y.intersp = 0.8
    )
  }
  if (!is.null(caption)) {
    # Push the caption below the legend row (rather than under the same
    # near-axis lines it would otherwise share) when both are drawn -
    # matches the extra room the combined bottom_margin branch above adds.
    caption_line_offset <- if (show_legend) 2.5 else 0
    mtext(caption, side = 1, line = seq_along(caption) * 1.3 - 0.7 + caption_line_offset,
          cex = 0.6, adj = 0)
  }
  invisible(coordinates)
}

# Plots `gra` split into its connected components (weak/undirected sense)
# instead of one shared layout, and returns a ready-to-write node table.
#
# Filtering a Markov network down to a subset of edges (multi-loop-only,
# inter-event-only, ...) frequently leaves it not fully connected - e.g. a
# pair of sub-categories that only ever transition to each other, nowhere
# near the rest of the network in the data. Forcing
# compute_markov_layout()'s force-directed layout to place every component
# in one shared coordinate space stretches the whole plot to make room for
# these outliers, at the cost of legibility for the (much larger) main
# component. Instead, the largest component gets a full-size main panel
# with its own independent layout, and every other component gets its own
# compact panel stacked below it (smallest scaled down further via
# half_width, since two points normalized to fill a wide box otherwise
# read as needlessly spread out) - so nothing is dropped from the plot,
# but nothing distorts the main network's layout either.
#
# Returns a data frame (id, sub_category, category, occurrences,
# has_significant_self_transition, component, component_size, x, y), one
# row per plotted node, sorted by component (largest first) then id -
# callers write this straight to a nodes CSV rather than repeating the
# per-component layout/attribute bookkeeping themselves.
plot_network_by_component <- function(
    gra, path, width = 9, height = 9, list_names = LIST_NAMES,
    seed = PLOT_SEED, caption = NULL
) {
  components_info <- components(gra, mode = "weak")
  component_order <- order(-components_info$csize)
  n_components <- components_info$no
  message(
    n_components, " connected component(s) in the plotted network (sizes: ",
    paste(components_info$csize[component_order], collapse = ", "), ")"
  )

  component_graphs <- lapply(component_order, function(component_id) {
    induced_subgraph(gra, V(gra)[components_info$membership == component_id])
  })
  main_size <- components_info$csize[component_order[1]]
  component_coordinates <- lapply(component_graphs, function(g_i) {
    half_width <- max(0.3, 1.3 * sqrt(vcount(g_i) / main_size))
    compute_markov_layout(g_i, seed = seed, half_width = half_width)
  })
  # plot_markov_graph() defaults weight_reference_range to that call's own
  # range(E(gra)$weight), which is degenerate (min == max, so
  # rescale_from_reference() errors) whenever a component has only a
  # single edge - true of every small disconnected component typically
  # seen here. Pass one shared range, computed across all edges before the
  # split, to every panel so edge widths stay on the same, non-degenerate
  # scale.
  weight_reference_range <- range(E(gra)$weight)

  main_caption <- caption
  if (n_components > 1) {
    main_caption <- c(main_caption, paste0(
      "largest connected component shown above (", vcount(component_graphs[[1]]),
      " / ", vcount(gra), " plotted nodes); ",
      n_components - 1L, " smaller disconnected component(s) shown below"
    ))
  }

  with_pdf_plot(
    path,
    width = width,
    height = height,
    {
      if (n_components == 1) {
        plot_markov_graph(
          component_graphs[[1]],
          list_names = list_names,
          seed = seed,
          coordinates = component_coordinates[[1]],
          weight_reference_range = weight_reference_range,
          show_legend = FALSE,
          caption = main_caption
        )
      } else {
        n_other <- n_components - 1L
        # Row 1 = the main component, spanning every column; row 2 = one
        # column per smaller component, sized well below the main row.
        layout(
          matrix(
            c(rep(1L, n_other), seq_len(n_other) + 1L),
            nrow = 2, byrow = TRUE
          ),
          heights = c(5, 1.4)
        )
        plot_markov_graph(
          component_graphs[[1]],
          list_names = list_names,
          seed = seed,
          coordinates = component_coordinates[[1]],
          weight_reference_range = weight_reference_range,
          show_legend = FALSE,
          caption = main_caption
        )
        for (i in seq_len(n_other)) {
          g_i <- component_graphs[[i + 1L]]
          plot_markov_graph(
            g_i,
            list_names = list_names,
            seed = seed,
            coordinates = component_coordinates[[i + 1L]],
            weight_reference_range = weight_reference_range,
            show_legend = FALSE,
            # A `main=` title uses plot.igraph()'s default title cex, sized
            # for the large main panel - comically oversized on these
            # small ones. caption's smaller, fixed cex reads correctly
            # regardless of panel size.
            caption = paste0(
              "component ", i + 1L, " (", vcount(g_i), " node",
              if (vcount(g_i) > 1) "s" else "", ")"
            )
          )
        }
      }
    }
  )

  node_table <- do.call(rbind, lapply(seq_along(component_graphs), function(i) {
    g_i <- component_graphs[[i]]
    coords_i <- component_coordinates[[i]]
    data.frame(
      node_id = as.integer(V(g_i)$name),
      sub_category = V(g_i)$sub_category,
      category = V(g_i)$category,
      occurrences = V(g_i)$occurrences,
      has_significant_self_transition = V(g_i)$vertex.frame.col == "black",
      component = i,
      component_size = vcount(g_i),
      x = coords_i[, "x"],
      y = coords_i[, "y"],
      stringsAsFactors = FALSE
    )
  }))
  node_table[order(node_table$component, node_table$node_id), ]
}

order_nodes_by_category <- function(df, list_names, order_col) {
  nodes_ordered <- c()
  for (i in seq_along(list_names)) {
    sub_df <- df[df$category == list_names[i], , drop = FALSE]
    if (nrow(sub_df) == 0) next
    ord <- order(sub_df[[order_col]], decreasing = TRUE, na.last = TRUE)
    nodes_ordered <- c(nodes_ordered, as.character(sub_df$nodes)[ord])
  }
  df$nodes <- factor(as.character(df$nodes), levels = unique(nodes_ordered))
  df
}

# Scatter of node betweenness vs. total strength, each normalized to [0, 1]
# by dividing by its own max within `gra`, with a marginal histogram (plus
# a fitted-exponential density curve, rate = 1 / mean) along each axis -
# the two centrality measures are typically right-skewed/heavy-tailed, and
# an exponential is the simplest single-parameter model for that shape.
# Red dashed lines mark each measure's empirical 90th percentile; a node is
# labeled with its id when it exceeds *either* line (an outlier on
# betweenness, strength, or both) - unremarkable nodes clustered near the
# origin are left unlabeled to keep that region legible.
plot_betweenness_strength_scatter <- function(
    gra, path, width = 8, height = 8, percentile = 0.9, n_bins = 10,
    # Per-vertex outlier-label color/font override (see plot_markov_graph()'s
    # same-purpose parameters) - NULL (default) reproduces the prior
    # behavior (plain black labels) exactly for every existing caller.
    vertex_label_color = NULL, vertex_label_font = NULL,
    # Per-vertex outlier-label TEXT override (see plot_markov_graph()'s
    # same-purpose parameter for the full rationale) - NULL (default)
    # reproduces the prior behavior (V(gra)$name) exactly.
    vertex_label = NULL,
    # Per-vertex point size / border width override, same NULL-preserves-
    # prior-behavior convention (see plot_markov_graph()'s vertex_size/
    # vertex_frame_width) - lets a caller mixing more than one kind of node
    # (scripts/10) make one kind's points visibly bigger/heavier-bordered
    # than the other. NULL reproduces the original constants (cex = 1.6,
    # point border lwd = 1) for every other caller.
    vertex_cex = NULL, vertex_lwd = NULL,
    # Extra caption lines appended below plot_markov_graph()-style caption's
    # usual two (percentile-threshold, labeled-node) lines - e.g. explaining
    # a non-default vertex_cex/vertex_lwd encoding. NULL omits them, leaving
    # every other caller's caption exactly as it was.
    extra_caption = NULL,
    # Precomputed, unnormalized per-vertex strength, in the same vertex
    # order as V(gra) (i.e. betweenness(gra)'s own order) - lets a caller
    # plot a different strength definition (e.g. raw transition counts with
    # self-transitions included, divided by occurrence count) on the same
    # normalized/marginal-histogram/outlier-labeling layout instead of this
    # function's own strength(gra, mode = "all"). NULL (default) reproduces
    # every existing caller's behavior exactly. strength_label customizes
    # the y-axis text to match.
    strength_values = NULL, strength_label = "strength (normalized)",
    # FALSE drops the two-line explanatory caption below the plot (and
    # shrinks the bottom margin that made room for it) - TRUE (default)
    # reproduces every existing caller's behavior exactly.
    show_caption = TRUE,
    # Multiplier applied to every text element's cex (axis titles, axis tick
    # labels, outlier point labels, caption) - 1 (default) reproduces every
    # existing caller's sizing exactly; a caller preparing a figure for a
    # smaller reproduction (e.g. Fig. 5c) can bump this up
    # without touching point/marker size (vertex_cex/vertex_lwd above are
    # separate).
    text_cex = 1,
    # Fraction of each axis's data range added as blank margin on both ends
    # of that axis - 0 (default) reproduces every existing caller's exact-
    # fit xaxs/yaxs = "i" framing, which puts a point exactly at the data
    # min/max (e.g. betweenness = 0 or 1) right on the panel border, half
    # clipped by it. A small positive value pushes the frame out so those
    # points sit fully inside it instead. Only the axis/window limits move -
    # histogram bin edges (breaks) stay anchored to the true data range, so
    # bin content is unaffected.
    axis_padding = 0
) {
  if (!is.null(strength_values) && length(strength_values) != vcount(gra)) {
    stop("strength_values must have one entry per vertex of gra", call. = FALSE)
  }
  betweenness_norm <- betweenness(gra)
  betweenness_norm <- betweenness_norm / max(betweenness_norm)
  strength_raw <- if (is.null(strength_values)) strength(gra, mode = "all") else strength_values

  # The default strength(gra, mode = "all") always normalizes to [0, 1] by
  # dividing by its max ("default" mode below), exactly reproducing every
  # existing caller's prior behavior - valid because it's always
  # non-negative with a meaningful floor at 0. A caller-supplied
  # strength_values instead gets one of two treatments depending on
  # whether it can go negative:
  #  - non-negative (e.g. a per-occurrence rate): "raw" mode, left in its
  #    own units with the axis spanning its own observed [min, max], since
  #    those units (e.g. "transitions per occurrence") are meaningful to
  #    read directly and a floor of 0 may sit well below the lowest point
  #    actually observed.
  #  - can be negative (e.g. a null-model z-score): "minmax" mode,
  #    min-max normalized to [0, 1] over its own observed range - dividing
  #    by max isn't valid here (it wouldn't map negative values into
  #    [0, 1] at all), but min-max rescaling still maps the minimum
  #    actually observed to axis position 0, keeping both axes on a common
  #    [0, 1] scale.
  use_default_strength <- is.null(strength_values)
  strength_nonneg <- min(strength_raw) >= 0
  strength_scale_mode <- if (use_default_strength) {
    "default"
  } else if (strength_nonneg) {
    "raw"
  } else {
    "minmax"
  }
  if (strength_scale_mode == "raw") {
    strength_norm <- strength_raw
    strength_ylim <- range(strength_raw)
  } else {
    strength_norm <- if (strength_scale_mode == "default") {
      strength_raw / max(strength_raw)
    } else {
      strength_range <- range(strength_raw)
      (strength_raw - strength_range[1]) / diff(strength_range)
    }
    strength_ylim <- c(0, 1)
  }

  betweenness_threshold <- quantile(betweenness_norm, percentile, names = FALSE)
  strength_threshold <- quantile(strength_norm, percentile, names = FALSE)
  labeled <- betweenness_norm > betweenness_threshold |
    strength_norm > strength_threshold

  # Exponential MLE: rate = 1 / mean. Guarded against a degenerate
  # all-zero measure (e.g. a graph with no betweenness-carrying paths),
  # which would otherwise divide by zero.
  betweenness_rate <- if (mean(betweenness_norm) > 0) {
    1 / mean(betweenness_norm)
  } else {
    NA
  }

  breaks <- seq(0, 1, length.out = n_bins + 1)
  curve_x <- seq(0, 1, length.out = 200)
  strength_breaks <- if (strength_scale_mode == "raw") {
    seq(strength_ylim[1], strength_ylim[2], length.out = n_bins + 1)
  } else {
    breaks
  }
  strength_curve_x <- if (strength_scale_mode == "raw") {
    seq(strength_ylim[1], strength_ylim[2], length.out = 200)
  } else {
    curve_x
  }
  # Exponential fit whenever strength is non-negative - its support is
  # [0, Inf), the standard model for a right-skewed centrality-like
  # measure, and applies equally to the default strength(gra, "all") and a
  # non-negative custom strength_values (e.g. a per-occurrence rate).
  # Falls back to a normal fit (mean/sd) once negative values are possible
  # (e.g. a null-model z-score), since an exponential can't represent
  # negative mass at all.
  strength_fit_y <- if (strength_nonneg) {
    strength_rate <- if (mean(strength_norm) > 0) 1 / mean(strength_norm) else NA
    if (is.finite(strength_rate)) dexp(strength_curve_x, rate = strength_rate) else NULL
  } else {
    strength_sd <- sd(strength_norm)
    if (strength_sd > 0) {
      dnorm(strength_curve_x, mean = mean(strength_norm), sd = strength_sd)
    } else {
      NULL
    }
  }

  # Padded axis/window limits (see axis_padding above) - betweenness's data
  # range is always [0, 1]; strength's is strength_ylim, whichever mode
  # produced it.
  x_lim <- c(0, 1) + c(-1, 1) * axis_padding * 1
  y_lim <- strength_ylim + c(-1, 1) * axis_padding * diff(strength_ylim)

  # Shared margins keep the three panels aligned: top-hist and main share
  # left/right (both in layout's left column); right-hist and main share
  # bottom/top (both in layout's bottom row).
  left_margin <- 4.5
  right_margin <- 0.5
  bottom_margin <- if (show_caption) 7.2 + length(extra_caption) * 1.3 else 4.5
  top_margin <- 1
  point_cex <- if (is.null(vertex_cex)) 1.6 else vertex_cex
  point_lwd <- if (is.null(vertex_lwd)) 1 else vertex_lwd

  with_pdf_plot(
    path,
    width = width,
    height = height,
    {
      layout(
        matrix(c(2, 0, 1, 3), nrow = 2, byrow = TRUE),
        widths = c(4, 1), heights = c(1, 4)
      )
      # Applies to every subsequent base-graphics call in this device (plot,
      # hist, axis) that doesn't set its own cex.axis/cex.lab explicitly.
      par(cex.axis = text_cex, cex.lab = text_cex)

      par(mar = c(bottom_margin, left_margin, top_margin - 0.5, right_margin))
      plot(
        betweenness_norm, strength_norm,
        xlim = x_lim, ylim = y_lim, xaxs = "i", yaxs = "i",
        pch = 21, bg = V(gra)$color, col = "black", cex = point_cex, lwd = point_lwd,
        xlab = "betweenness (normalized)", ylab = strength_label,
        las = 1, yaxt = if (strength_scale_mode == "raw") "n" else "s"
      )
      if (strength_scale_mode == "raw") {
        # yaxs = "i" already makes the panel span exactly strength_ylim
        # (the true observed min/max) with no padding, but axis()'s default
        # tick placement only labels "pretty" round numbers inside that
        # range - since the true min/max are rarely round, neither gets its
        # own label, which visually reads as if the axis starts/stops short
        # of where the lowest/highest point actually is. Forcing the exact
        # (rounded for display) min and max in as their own ticks, on top
        # of the usual pretty in-between ones, makes the axis's true extent
        # unambiguous.
        strength_axis_at <- sort(unique(c(
          round(strength_ylim, 2), pretty(strength_ylim)
        )))
        strength_axis_at <- strength_axis_at[
          strength_axis_at >= strength_ylim[1] & strength_axis_at <= strength_ylim[2]
        ]
        axis(2, at = strength_axis_at, las = 1)
      }
      abline(v = betweenness_threshold, col = "red", lty = 3, lwd = 1.5)
      abline(h = strength_threshold, col = "red", lty = 3, lwd = 1.5)
      if (any(labeled)) {
        label_color <- if (is.null(vertex_label_color)) "black" else vertex_label_color[labeled]
        label_font <- if (is.null(vertex_label_font)) 1 else vertex_label_font[labeled]
        label_text <- if (is.null(vertex_label)) V(gra)$name else vertex_label
        text(
          betweenness_norm[labeled], strength_norm[labeled],
          labels = label_text[labeled], pos = 3, cex = 0.8 * text_cex, xpd = NA,
          col = label_color, font = label_font
        )
      }
      if (show_caption) {
        strength_norm_caption <- switch(
          strength_scale_mode,
          default = "betweenness/strength each normalized by their own max in this network",
          raw = "betweenness normalized by its own max; strength (y-axis) left in its own raw units, axis spans its observed [min, max]",
          minmax = "betweenness normalized by its own max; strength (y-axis) min-max normalized to [0, 1] over its own observed range (can be negative, so not divided by max)"
        )
        strength_fit_caption <- if (strength_nonneg) {
          "blue curve on marginal histograms = fitted exponential density (rate = 1 / mean)"
        } else {
          "blue curve on betweenness histogram = fitted exponential density; on strength histogram = fitted normal density (mean, sd)"
        }
        mtext(
          c(
            paste0(
              strength_norm_caption, "; red dashed lines = ",
              percentile * 100, "th percentile of each"
            ),
            paste0("labeled nodes exceed either threshold; ", strength_fit_caption),
            extra_caption
          ),
          side = 1, line = c(4.2, 5.5, 5.5 + 1.3 * seq_along(extra_caption)),
          cex = 0.55 * text_cex, adj = 0
        )
      }

      par(mar = c(top_margin - 0.5, left_margin, top_margin, right_margin))
      hist(
        betweenness_norm, breaks = breaks, freq = FALSE,
        col = "gray70", border = "white", main = "", xlab = "", ylab = "",
        xlim = x_lim, xaxs = "i", xaxt = "n", las = 1
      )
      if (is.finite(betweenness_rate)) {
        lines(curve_x, dexp(curve_x, rate = betweenness_rate),
              col = "steelblue", lwd = 2)
      }

      par(mar = c(bottom_margin, right_margin, top_margin - 0.5, right_margin + 0.5))
      h_right <- hist(strength_norm, breaks = strength_breaks, plot = FALSE)
      right_xmax <- max(h_right$density, if (!is.null(strength_fit_y)) max(strength_fit_y) else 0)
      plot.new()
      plot.window(xlim = c(0, right_xmax * 1.05), ylim = y_lim, yaxs = "i")
      rect(
        xleft = 0, ybottom = h_right$breaks[-length(h_right$breaks)],
        xright = h_right$density, ytop = h_right$breaks[-1],
        col = "gray70", border = "white"
      )
      if (!is.null(strength_fit_y)) {
        lines(strength_fit_y, strength_curve_x, col = "steelblue", lwd = 2)
      }
      axis(1)
      box()
    }
  )
}

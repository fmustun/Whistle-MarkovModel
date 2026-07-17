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

compute_graph_p_value_significant <- function(
    p_value,
    Markov_Model,
    Whistles_List,
    sub_division = TRUE,
    sub_network = NULL,
    list_names = LIST_NAMES
) {
  seq_category <- build_seq_category(Whistles_List, list_names)

  if (!is.null(sub_network)) {
    subv <- seq_category[[sub_network]]
    num_cat <- max(unique(Whistles_List$whistle_type))
    transition_probabilities_matrix_p_value <- matrix(0, num_cat, num_cat)
    transition_probabilities_matrix_p_value[subv, ] <-
      Markov_Model$transition_probabilities_matrix_all[subv, ]
    transition_probabilities_matrix_p_value[, subv] <-
      Markov_Model$transition_probabilities_matrix_all[, subv]
    transition_probabilities_matrix_p_value[
      Markov_Model$p_value_matrix >= p_value
    ] <- 0
  } else {
    transition_probabilities_matrix_p_value <-
      Markov_Model$transition_probabilities_matrix_all
    transition_probabilities_matrix_p_value[
      Markov_Model$p_value_matrix >= p_value
    ] <- 0
  }

  edges_p_value <- c()
  for (i in seq_len(dim(transition_probabilities_matrix_p_value)[1])) {
    v_tmp <- Markov_Model$inv_p_value_matrix[i, ][
      transition_probabilities_matrix_p_value[i, ] > 0
    ]
    edges_p_value <- c(edges_p_value, v_tmp)
  }

  Gra <- compute_graph(
    transition_probabilities_matrix_p_value,
    Markov_Model = Markov_Model,
    Whistles_List = Whistles_List,
    sub_division = sub_division,
    list_names = list_names
  )
  Gra <- set_edge_attr(Gra, "inv_p_value", index = E(Gra), edges_p_value)

  v.frame <- ifelse(
    V(Gra)$name %in% which(diag(transition_probabilities_matrix_p_value) > 0),
    "black", "gray"
  )
  Gra <- set_vertex_attr(Gra, "vertex.frame.col", index = V(Gra), v.frame)
  Gra <- delete.edges(Gra, which(which_loop(Gra)))

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
  bottom_margin <- if (show_legend) {
    3
  } else if (!is.null(caption)) {
    length(caption) * 1.3 + 0.5
  } else {
    0.5
  }
  old_par <- par(mar = c(bottom_margin, 0.5, 0.5, 0.5))
  on.exit(par(old_par), add = TRUE)
  if (is.null(weight_reference_range)) {
    weight_reference_range <- range(E(gra)$weight)
  }
  edge_alpha <- rescale_from_reference(
    E(gra)$weight, weight_reference_range, c(0.1, 1)
  )
  edge_width <- rescale_from_reference(
    E(gra)$weight, weight_reference_range, c(3, 7)
  )
  edge_alpha <- pmin(1, pmax(0.01, edge_alpha))
  edge_width <- pmax(0, edge_width)

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
    vertex.size = 5,
    edge.color = rgb(
      140 / 255, 140 / 255, 140 / 255,
      edge_alpha
    ),
    edge.width = edge_width,
    edge.arrow.size = 0.95,
    edge.arrow.width = 0.95,
    vertex.shape = "fcircle",
    vertex.frame.color = V(gra)$vertex.frame.col,
    vertex.frame.width = 2.5,
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
    mtext(caption, side = 1, line = seq_along(caption) * 1.3 - 0.7,
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
    gra, path, width = 8, height = 8, percentile = 0.9, n_bins = 10
) {
  betweenness_norm <- betweenness(gra)
  betweenness_norm <- betweenness_norm / max(betweenness_norm)
  strength_norm <- strength(gra, mode = "all")
  strength_norm <- strength_norm / max(strength_norm)

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
  strength_rate <- if (mean(strength_norm) > 0) 1 / mean(strength_norm) else NA

  breaks <- seq(0, 1, length.out = n_bins + 1)
  curve_x <- seq(0, 1, length.out = 200)

  # Shared margins keep the three panels aligned: top-hist and main share
  # left/right (both in layout's left column); right-hist and main share
  # bottom/top (both in layout's bottom row).
  left_margin <- 4.5
  right_margin <- 0.5
  bottom_margin <- 7.2
  top_margin <- 1

  with_pdf_plot(
    path,
    width = width,
    height = height,
    {
      layout(
        matrix(c(2, 0, 1, 3), nrow = 2, byrow = TRUE),
        widths = c(4, 1), heights = c(1, 4)
      )

      par(mar = c(bottom_margin, left_margin, top_margin - 0.5, right_margin))
      plot(
        betweenness_norm, strength_norm,
        xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
        pch = 21, bg = V(gra)$color, col = "black", cex = 1.6,
        xlab = "betweenness (normalized)", ylab = "strength (normalized)",
        las = 1
      )
      abline(v = betweenness_threshold, col = "red", lty = 3, lwd = 1.5)
      abline(h = strength_threshold, col = "red", lty = 3, lwd = 1.5)
      if (any(labeled)) {
        text(
          betweenness_norm[labeled], strength_norm[labeled],
          labels = V(gra)$name[labeled], pos = 3, cex = 0.8, xpd = NA
        )
      }
      mtext(
        c(
          paste0(
            "betweenness/strength each normalized by their own max in this network; red dashed lines = ",
            percentile * 100, "th percentile of each"
          ),
          "labeled nodes exceed either threshold; blue curve on marginal histograms = fitted exponential density (rate = 1 / mean)"
        ),
        side = 1, line = c(4.2, 5.5), cex = 0.55, adj = 0
      )

      par(mar = c(top_margin - 0.5, left_margin, top_margin, right_margin))
      hist(
        betweenness_norm, breaks = breaks, freq = FALSE,
        col = "gray70", border = "white", main = "", xlab = "", ylab = "",
        xlim = c(0, 1), xaxs = "i", xaxt = "n", las = 1
      )
      if (is.finite(betweenness_rate)) {
        lines(curve_x, dexp(curve_x, rate = betweenness_rate),
              col = "steelblue", lwd = 2)
      }

      par(mar = c(bottom_margin, right_margin, top_margin - 0.5, right_margin + 0.5))
      h_right <- hist(strength_norm, breaks = breaks, plot = FALSE)
      right_xmax <- max(h_right$density, if (is.finite(strength_rate)) {
        dexp(0, rate = strength_rate)
      } else {
        0
      })
      plot.new()
      plot.window(xlim = c(0, right_xmax * 1.05), ylim = c(0, 1), yaxs = "i")
      rect(
        xleft = 0, ybottom = h_right$breaks[-length(h_right$breaks)],
        xright = h_right$density, ytop = h_right$breaks[-1],
        col = "gray70", border = "white"
      )
      if (is.finite(strength_rate)) {
        lines(dexp(curve_x, rate = strength_rate), curve_x,
              col = "steelblue", lwd = 2)
      }
      axis(1)
      box()
    }
  )
}

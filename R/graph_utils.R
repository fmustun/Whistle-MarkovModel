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
  vertex.size <- 1 / 200 * params("vertex", "size")
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

compute_markov_layout <- function(gra, seed = PLOT_SEED) {
  set.seed(seed)
  coordinates <- layout_with_fr(gra, niter = 500, grid = "nogrid")
  coordinates <- norm_coords(
    coordinates, xmin = -0.85, xmax = 0.85, ymin = -0.85, ymax = 0.85
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
    main = NULL
) {
  if (is.null(coordinates)) coordinates <- compute_markov_layout(gra, seed)
  coordinates <- layout_for_graph(coordinates, gra)
  if (is.null(weight_reference_range)) {
    weight_reference_range <- range(E(gra)$weight)
  }
  edge_alpha <- rescale_from_reference(
    E(gra)$weight, weight_reference_range, c(0.01, 1)
  )
  edge_width <- rescale_from_reference(
    E(gra)$weight, weight_reference_range, c(3, 7)
  )
  edge_alpha <- pmin(1, pmax(0.01, edge_alpha))
  edge_width <- pmax(0, edge_width)

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
    edge.arrow.size = 0.5,
    edge.arrow.width = 0.75,
    vertex.shape = "fcircle",
    vertex.frame.color = V(gra)$vertex.frame.col,
    vertex.frame.width = 2.5,
    ylim = c(-0.85, 0.85),
    xlim = c(-0.85, 0.85),
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
  invisible(coordinates)
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

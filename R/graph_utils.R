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

plot_markov_graph <- function(gra, list_names = LIST_NAMES, seed = PLOT_SEED) {
  set.seed(seed)
  plot.igraph(
    gra,
    layout = layout_with_fr(gra, niter = 500, grid = "nogrid"),
    edge.curved = 0.2,
    vertex.size = 5,
    edge.color = rgb(
      140 / 255, 140 / 255, 140 / 255,
      rescale(E(gra)$weight, c(0.01, 1))
    ),
    edge.width = rescale(E(gra)$weight, c(3, 7)),
    edge.arrow.size = 0.5,
    edge.arrow.width = 0.75,
    vertex.shape = "fcircle",
    vertex.frame.color = V(gra)$vertex.frame.col,
    vertex.frame.width = 2.5,
    ylim = c(-0.85, 0.85),
    xlim = c(-0.85, 0.85),
    asp = 1
  )
  a <- legend(-2, 1, legend = list_names)
  x <- (a$text$x + a$rect$left) / 2
  y <- a$text$y
  symbols(
    x, y,
    circles = rep(1 / 30, length(list_names)),
    inches = FALSE, add = TRUE,
    bg = unique(V(gra)$color), col = "gray"
  )
  invisible(gra)
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

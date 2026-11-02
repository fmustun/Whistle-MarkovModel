#!/usr/bin/env Rscript

# Which whistle sub-categories tend to occur together in the same
# multi-loop chain (IWI < MULTILOOP_IWI_THRESHOLD), ignoring order/position
# entirely? This is a different question from the two multi-loop-pattern
# analyses already in the repo:
#  - compute_multiloop_edge_counts() (script 06) only counts *adjacent*
#    different-type pairs (A -> B right next to each other).
#  - compute_multiloop_patterns() (script 04/07) matches a chain's whole
#    composition as one unit (e.g. "A + A + B"), so "A + B + C" and
#    "A + B" never contribute to the same pattern.
# Here, any two (or three) distinct sub-categories that both appear
# somewhere in the same chain count as "co-occurring", regardless of how
# many whistles are in the chain, their order, or repeats. Built from
# compute_multiloop_type_membership()/compute_type_cooccurrence_pairs()/
# _triplets() (R/multiloop_utils.R).

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
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

register_fcircle_shape()

message("Loading whistles from ", DATA_PATH)
whistles_list <- load_whistles(DATA_PATH)

message(
  "Grouping whistles into multi-loop chains (IWI < ",
  MULTILOOP_IWI_THRESHOLD * 1000, "ms)..."
)
whistle_chains <- compute_multiloop_chains(
  whistles_list, threshold = MULTILOOP_IWI_THRESHOLD, list_names = LIST_NAMES
)

membership <- compute_multiloop_type_membership(whistle_chains)
n_multiloop_chains <- length(unique(membership$chain_id))
message(
  n_multiloop_chains, " multi-loop chains (chain_length >= 2) contribute ",
  nrow(membership), " (chain, sub-category) memberships"
)

# whistle_type (int) <-> whistle_type_chr / category lookups, built once and
# reused for every table/plot below.
type_lookup <- unique(whistle_chains[, c("whistle_type", "whistle_type_chr", "category")])
type_chr_of <- setNames(type_lookup$whistle_type_chr, type_lookup$whistle_type)
category_of <- setNames(type_lookup$category, type_lookup$whistle_type)
category_color_of <- setNames(LIST_COLORS[seq_along(LIST_NAMES)], LIST_NAMES)

# How often each sub-category appears in a multi-loop chain at all
# (marginal frequency) - used for node sizing below.
marginal_counts <- table(membership$whistle_type)

# -- pairwise / triplet co-occurrence tables --------------------------------

pair_counts <- compute_type_cooccurrence_pairs(membership)
pair_counts$type_a_chr <- type_chr_of[as.character(pair_counts$type_a)]
pair_counts$type_b_chr <- type_chr_of[as.character(pair_counts$type_b)]
pair_counts$category_a <- category_of[as.character(pair_counts$type_a)]
pair_counts$category_b <- category_of[as.character(pair_counts$type_b)]
pair_counts <- pair_counts[
  , c(
    "type_a", "type_a_chr", "category_a",
    "type_b", "type_b_chr", "category_b", "n_chains"
  )
]

triplet_counts <- compute_type_cooccurrence_triplets(membership)
triplet_counts$type_a_chr <- type_chr_of[as.character(triplet_counts$type_a)]
triplet_counts$type_b_chr <- type_chr_of[as.character(triplet_counts$type_b)]
triplet_counts$type_c_chr <- type_chr_of[as.character(triplet_counts$type_c)]
triplet_counts$category_a <- category_of[as.character(triplet_counts$type_a)]
triplet_counts$category_b <- category_of[as.character(triplet_counts$type_b)]
triplet_counts$category_c <- category_of[as.character(triplet_counts$type_c)]
triplet_counts <- triplet_counts[
  , c(
    "type_a", "type_a_chr", "category_a",
    "type_b", "type_b_chr", "category_b",
    "type_c", "type_c_chr", "category_c", "n_chains"
  )
]

message(
  nrow(pair_counts), " distinct co-occurring pairs, ",
  nrow(triplet_counts), " distinct co-occurring triplets observed ",
  "(order/position ignored)"
)

pair_csv <- file.path(multiloop_dir, "multiloop_type_pair_cooccurrence.csv")
write.csv(pair_counts, pair_csv, row.names = FALSE)
message("Saved ", pair_csv)

triplet_csv <- file.path(multiloop_dir, "multiloop_type_triplet_cooccurrence.csv")
write.csv(triplet_counts, triplet_csv, row.names = FALSE)
message("Saved ", triplet_csv)

# -- distribution plot: top co-occurring pairs ------------------------------

# "Which types most commonly co-occur" is exactly the top of this ranked
# list; capped at 25 so the bar labels stay legible (pair_counts itself, the
# full ranked table, is already saved above).
top_n_pairs <- min(25, nrow(pair_counts))
top_pairs <- head(pair_counts, top_n_pairs)
top_pairs$pair_label <- paste(top_pairs$type_a_chr, top_pairs$type_b_chr, sep = " + ")
top_pairs$pair_label <- factor(top_pairs$pair_label, levels = rev(top_pairs$pair_label))
# A pair's two sub-categories either share one main category or span two;
# colored the same way script 07 colors mixed-category multi-loop-type
# nodes (single category -> that category's color, mixed -> a fixed gray).
MIXED_CATEGORY_COLOR <- "gray50"
top_pairs$pair_category <- ifelse(
  top_pairs$category_a == top_pairs$category_b, top_pairs$category_a, "mixed"
)

pair_distribution_plot <- ggplot(
  top_pairs, aes(x = pair_label, y = n_chains, fill = pair_category)
) +
  geom_col(color = "black", width = 0.75) +
  geom_text(aes(label = n_chains), hjust = -0.2, size = 3) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    breaks = c(LIST_NAMES, "mixed"),
    values = c(category_color_of, mixed = MIXED_CATEGORY_COLOR)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_minimal() +
  theme(plot.margin = margin(5.5, 20, 5.5, 5.5), plot.title.position = "plot") +
  labs(
    title = paste0(
      "Top ", top_n_pairs, " co-occurring sub-category pairs within a multi-loop chain"
    ),
    subtitle = paste0(
      "order/position ignored; out of ", nrow(pair_counts),
      " distinct pairs observed across ", n_multiloop_chains, " multi-loop chains"
    ),
    x = NULL,
    y = "number of multi-loop chains containing both sub-categories",
    fill = "category"
  )

save_ggplot(
  pair_distribution_plot, "multiloop_type_cooccurrence_top_pairs", REPO_ROOT,
  width = 9, height = 8, subdir = "multiloops", bg = "white"
)

# -- distribution plot: top co-occurring triplets ---------------------------

top_n_triplets <- min(25, nrow(triplet_counts))
top_triplets <- head(triplet_counts, top_n_triplets)
top_triplets$triplet_label <- paste(
  top_triplets$type_a_chr, top_triplets$type_b_chr, top_triplets$type_c_chr,
  sep = " + "
)
top_triplets$triplet_label <- factor(
  top_triplets$triplet_label, levels = rev(top_triplets$triplet_label)
)
# Same single-category-vs-mixed coloring as the pair distribution above, just
# extended to all three sub-categories agreeing rather than two.
top_triplets$triplet_category <- ifelse(
  top_triplets$category_a == top_triplets$category_b &
    top_triplets$category_b == top_triplets$category_c,
  top_triplets$category_a, "mixed"
)

triplet_distribution_plot <- ggplot(
  top_triplets, aes(x = triplet_label, y = n_chains, fill = triplet_category)
) +
  geom_col(color = "black", width = 0.75) +
  geom_text(aes(label = n_chains), hjust = -0.2, size = 3) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    breaks = c(LIST_NAMES, "mixed"),
    values = c(category_color_of, mixed = MIXED_CATEGORY_COLOR)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_minimal() +
  theme(plot.margin = margin(5.5, 20, 5.5, 5.5), plot.title.position = "plot") +
  labs(
    title = paste0(
      "Top ", top_n_triplets, " co-occurring sub-category triplets within a multi-loop chain"
    ),
    subtitle = paste0(
      "order/position ignored; out of ", nrow(triplet_counts),
      " distinct triplets observed across ", n_multiloop_chains, " multi-loop chains"
    ),
    x = NULL,
    y = "number of multi-loop chains containing all three sub-categories",
    fill = "category"
  )

save_ggplot(
  triplet_distribution_plot, "multiloop_type_cooccurrence_top_triplets", REPO_ROOT,
  width = 9, height = 8, subdir = "multiloops", bg = "white"
)

# -- hypergraph of pairs (2-uniform: an ordinary co-occurrence network) -----

# A single shared co-occurrence is easy to hit by chance in a chain of 2-6
# whistles; requiring it in >= 2 chains before drawing an edge is the same
# "recurs" bar script 04's summary already uses for repeated sub-category
# patterns (n_occurrences > 1) - the minimum needed to call something a
# recurring pairing rather than a one-off.
pairs_recurring <- pair_counts[pair_counts$n_chains >= 2, ]
message(
  nrow(pairs_recurring), " / ", nrow(pair_counts),
  " pairs recur (co-occur in >= 2 multi-loop chains); plotting those as the pair network"
)

gra_pairs <- graph_from_data_frame(
  data.frame(
    from = as.character(pairs_recurring$type_a),
    to = as.character(pairs_recurring$type_b)
  ),
  directed = FALSE
)
E(gra_pairs)$weight <- pairs_recurring$n_chains
V(gra_pairs)$sub_category <- type_chr_of[V(gra_pairs)$name]
V(gra_pairs)$occurrences <- as.numeric(marginal_counts[V(gra_pairs)$name])
V(gra_pairs)$vertex.frame.col <- "gray"

seq_category <- build_seq_category(whistles_list, LIST_NAMES)
gra_pairs <- assign_vertex_category_color(gra_pairs, seq_category, LIST_NAMES, LIST_COLORS)

message("Plotting pair co-occurrence network...")
with_pdf_plot(
  plot_path("multiloop_type_pair_cooccurrence_network", REPO_ROOT, subdir = "multiloops"),
  width = 10,
  height = 10,
  plot_markov_graph(
    gra_pairs,
    seed = PLOT_SEED,
    show_legend = TRUE,
    caption = c(
      "each node = one whistle sub-category seen in a multi-loop chain; edge = both sub-categories co-occurred in the same chain, any order/position",
      paste0(
        "edges require co-occurrence in >= 2 distinct multi-loop chains (",
        nrow(pairs_recurring), " / ", nrow(pair_counts), " observed pairs); edge width = number of chains sharing that pair"
      ),
      "undirected: edge direction/arrowheads carry no meaning here"
    )
  )
)

# -- hypergraph of triplets (3-uniform, drawn as a bipartite incidence graph)

# igraph has no native hyperedge (a single edge spanning 3+ nodes at once),
# so each recurring triplet is instead drawn as its own small connector
# node wired to its 3 member sub-categories - the standard bipartite/
# incidence representation of a hypergraph. Same >= 2 chains "recurs" bar as
# the pair network above.
triplets_recurring <- triplet_counts[triplet_counts$n_chains >= 2, ]
message(
  nrow(triplets_recurring), " / ", nrow(triplet_counts),
  " triplets recur (co-occur in >= 2 multi-loop chains); plotting those as the triplet hypergraph"
)

if (nrow(triplets_recurring) == 0) {
  message("No recurring triplets to plot; skipping the triplet hypergraph.")
} else {
  hyperedge_id <- paste0("H", seq_len(nrow(triplets_recurring)))
  incidence_edges <- do.call(rbind, lapply(seq_len(nrow(triplets_recurring)), function(i) {
    data.frame(
      from = hyperedge_id[i],
      to = as.character(unlist(
        triplets_recurring[i, c("type_a", "type_b", "type_c")]
      ))
    )
  }))
  gra_triplets <- graph_from_data_frame(incidence_edges, directed = FALSE)

  is_hyperedge <- V(gra_triplets)$name %in% hyperedge_id
  node_color <- character(vcount(gra_triplets))
  node_color[!is_hyperedge] <- category_color_of[
    category_of[V(gra_triplets)$name[!is_hyperedge]]
  ]
  node_color[is_hyperedge] <- "gray20"

  node_shape <- ifelse(is_hyperedge, "square", "fcircle")
  node_frame <- ifelse(is_hyperedge, "black", "gray40")
  node_label <- ifelse(is_hyperedge, NA, V(gra_triplets)$name)

  node_size <- numeric(vcount(gra_triplets))
  node_size[!is_hyperedge] <- rescale_from_reference(
    as.numeric(marginal_counts[V(gra_triplets)$name[!is_hyperedge]]),
    range(marginal_counts), c(5, 11)
  )
  hyperedge_row <- match(V(gra_triplets)$name[is_hyperedge], hyperedge_id)
  node_size[is_hyperedge] <- rescale_from_reference(
    triplets_recurring$n_chains[hyperedge_row],
    range(triplets_recurring$n_chains), c(3, 7)
  )

  set.seed(PLOT_SEED)
  coordinates <- layout_with_fr(gra_triplets, niter = 2000, grid = "nogrid")
  coordinates <- norm_coords(coordinates, xmin = -1.3, xmax = 1.3, ymin = -1.3, ymax = 1.3)
  plot_lim <- range(coordinates) * 1.05

  message("Plotting triplet co-occurrence hypergraph...")
  with_pdf_plot(
    plot_path("multiloop_type_triplet_cooccurrence_hypergraph", REPO_ROOT, subdir = "multiloops"),
    width = 10,
    height = 10,
    {
      old_par <- par(mar = c(5.5, 0.5, 0.5, 0.5))
      on.exit(par(old_par), add = TRUE)
      plot.igraph(
        gra_triplets,
        layout = coordinates,
        rescale = FALSE,
        vertex.shape = node_shape,
        vertex.color = node_color,
        vertex.frame.color = node_frame,
        vertex.frame.width = 2,
        vertex.size = node_size,
        vertex.label = node_label,
        vertex.label.cex = 0.55,
        vertex.label.color = "black",
        edge.color = "gray60",
        edge.width = 1,
        edge.curved = 0,
        ylim = plot_lim,
        xlim = plot_lim,
        asp = 1
      )
      legend(
        "bottomleft", inset = c(0, -0.02), xpd = NA, ncol = 2,
        legend = LIST_NAMES, pch = 21, pt.bg = LIST_COLORS[seq_along(LIST_NAMES)],
        col = "gray", pt.cex = 1.05, cex = 0.55, bty = "n",
        x.intersp = 0.45, y.intersp = 0.8
      )
      legend(
        "bottomright", inset = c(0, -0.02), xpd = NA,
        legend = "triplet (size ~ n chains)", pch = 22, pt.bg = "gray20",
        col = "black", pt.cex = 1.3, cex = 0.55, bty = "n"
      )
      mtext(
        c(
          paste0(
            "circle = whistle sub-category (size ~ how many multi-loop chains it appears in); square = one recurring 3-way co-occurrence, wired to its 3 sub-categories"
          ),
          paste0(
            "edges require co-occurrence in >= 2 distinct multi-loop chains (",
            nrow(triplets_recurring), " / ", nrow(triplet_counts), " observed triplets); order/position ignored"
          )
        ),
        side = 1, line = c(1.3, 2.6), cex = 0.55, adj = 0
      )
    }
  )
}

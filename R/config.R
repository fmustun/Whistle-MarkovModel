# Analysis defaults (paths relative to repository root).

DATA_PATH <- "data/AllWhistlesSubClustering_final.csv"
OUTPUT_DIR <- "outputs"
PLOTS_DIR <- "outputs/plots"

LIST_NAMES <- c(
  "SW_Neo", "SW_Luna", "SW_Yosefa", "SW_Nikita", "SW_Nana",
  "SW_Dana", "SW_Shy", "NSW_3", "NSW_6", "NSW_9"
)

LIST_COLORS <- c(
  "#61D04F", "#2297E6", "#F5C710", "#ee0000", "#A6E1DF",
  "#ae3450", "#ffc2f5", "#a06f05", "#feaf92", "#e6e6fa"
)

TIME_WINDOW <- c(0, 5.94)
NULL_MODEL <- "shuffle"
GRAPH_P_VALUE <- 0.05
# Edge selection for a newly generated model. "none" reproduces the original
# raw empirical-p network; "BH" adjusts the family of transitions actually
# observed at least once (never-observed transitions are excluded from the
# correction family - see family_mask in adjust_transition_p_values,
# R/statistical_utils.R).
MULTIPLE_TESTING_CORRECTION <- "BH"
GRAPH_ALPHA <- 0.05
# Null-model permutations. The empirical p-value floor is 1/(ITERATIONS + 1),
# so BH correction across m observed transitions at GRAPH_ALPHA needs
# ITERATIONS >= m / GRAPH_ALPHA - 1 just to have any resolving power, plus a
# margin above that floor so p-value estimates near the decision boundary
# aren't dominated by Monte Carlo noise. recommended_iterations(m, alpha,
# margin) in R/statistical_utils.R computes this once you know m (the count
# of transitions observed at least once - see a prior run's
# outputs/bh_analysis/bh_analysis_summary.txt, or
# sum(MarkovModel$transition_probabilities_matrix_all > 0)). For this
# dataset m = 1214, so recommended_iterations(1214, 0.05, margin = 5) ~=
# 121,400; the manuscript run uses 150,000. Lower it (e.g. 1000) only for
# quick development runs - those are not statistically powered.
ITERATIONS <- 150000L
SEED <- 0L
MIN_SHIFT <- 0
MAX_SHIFT <- 1500
THRESHOLD <- 2
PLOT_SEED <- 42L

# Inter-whistle interval threshold (seconds) for grouping consecutive
# whistles into a multi-loop chain (start_time[i+1] - end_time[i] < this).
MULTILOOP_IWI_THRESHOLD <- 0.25

# Ordered-categorical clustering of multi-loop chain compositions
# (cluster_multiloop_chains(), R/chain_clustering_utils.R) - shared here
# rather than pinned locally in any one script, since both scripts/09's
# sub-clustering and scripts/10's cluster-level transition network need
# identical parameters to stay consistent with each other.
#
# EPSILON/CROSS_CATEGORY_COST are user-specified starting points, not
# derived: substitution cost between two sub-categories is 0 if identical,
# EPSILON if they share a main category, CROSS_CATEGORY_COST otherwise.
# Tying CROSS_CATEGORY_COST to GAP_PENALTY matters - leaving them unequal
# (e.g. CROSS_CATEGORY_COST = 0.6 with GAP_PENALTY = 1) was tried and
# fragmented clusters by chain length on top of main category, since a
# cross-category substitution became strictly cheaper than a gap; keeping
# them equal (as below) avoids that.
MULTILOOP_CLUSTER_EPSILON <- 0.3
MULTILOOP_CLUSTER_CROSS_CATEGORY_COST <- 1
MULTILOOP_CLUSTER_GAP_PENALTY <- 1
MULTILOOP_CLUSTER_LINKAGE <- "average"
MULTILOOP_CLUSTER_K_RANGE <- 5:40

# Leiden community detection (scripts/03, scripts/03b). COMMUNITY_RESOLUTION
# is the single resolution_parameter used for the community detection
# (Fig. 5a) - shared here (rather than pinned locally in scripts/03) so
# scripts/03b's resolution sweep can mark it on its plots and stay in sync.
# 0.5 was chosen from scripts/03b's sweep output
# (outputs/leiden_resolution_sweep.csv), not picked a priori: the
# number-of-communities curve is flat there (a stable 5-community plateau
# from resolution 0.25-0.55, before it starts climbing at 0.6), and the
# two-tailed shuffle-null benchmark at that resolution is significant
# (empirical_p_modularity_two_tailed < 0.05) - the real network is *less*
# modular than the null there, i.e. the real network supports a coarser,
# more cohesive partition than randomly reweighting the same significant
# edges would give.
#
# COMMUNITY_RESOLUTION_RANGE is the sweep range scripts/03b scans looking
# for that plateau. COMMUNITY_NULL_REALIZATIONS is how many of the
# ITERATIONS shuffle-null count matrices (already generated for the
# edge-significance test - MarkovModel$random_whistle_occurrences_matrices,
# R/markov_whistle.R) scripts/03b subsamples to benchmark the sweep against:
# running the full sweep across all ITERATIONS (150k) null realizations
# would be ~1500x more Leiden calls for no real gain in the resolution of
# the resulting null envelope, so a much smaller subsample (still large
# enough for a stable 95% band) is used instead.
COMMUNITY_RESOLUTION <- 0.5
COMMUNITY_RESOLUTION_RANGE <- seq(0.25, 2.5, by = 0.05)
COMMUNITY_NULL_REALIZATIONS <- 100L

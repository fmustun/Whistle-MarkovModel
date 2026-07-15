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
# dataset m ~= 1214, so recommended_iterations(1214, 0.05, margin = 5) ~=
# 121,400. 1000 below is a fast default for development iteration, not a
# statistically powered run.
ITERATIONS <- 150000L
SEED <- 0L
MIN_SHIFT <- 0
MAX_SHIFT <- 1500
THRESHOLD <- 2
NULL_NETWORK_ITERATIONS <- 1000L
PLOT_SEED <- 42L

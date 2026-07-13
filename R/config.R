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
  "#ae3450", "#ffc2f5", "#fefeb1", "#feaf92", "#e6e6fa"
)

TIME_WINDOW <- c(0, 5.94)
NULL_MODEL <- "shuffle"
GRAPH_P_VALUE <- 0.05
# Edge selection for a newly generated model. "none" reproduces the original
# raw empirical-p network; "BH" adjusts the complete ordered state-pair family.
MULTIPLE_TESTING_CORRECTION <- "BH"
GRAPH_ALPHA <- 0.05
ITERATIONS <- 1000L
SEED <- 0L
MIN_SHIFT <- 0
MAX_SHIFT <- 1500
THRESHOLD <- 2
NULL_NETWORK_ITERATIONS <- 1000L
PLOT_SEED <- 22L

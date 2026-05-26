# MarkovModel

Markov transition network analysis for dolphin whistle sub-categories: null-model significance, network visualization, node-level metrics, and global network statistics.

## Requirements

- R ≥ 4.0
- Packages: `igraph`, `ggplot2`, `scales`, `foreach`, `doSNOW`, `tnet`, `dplyr`, `gridExtra`, `ggrepel`, `ccber`

Install `ccber` (provides `CalcTransitionMatrix`):

```r
install.packages("remotes")
remotes::install_github("hadley/ccber")
```

Install other dependencies:

```r
install.packages(c("igraph", "ggplot2", "scales", "foreach", "doSNOW", "tnet",
                   "dplyr", "gridExtra", "ggrepel"))
```

## Data

See [data/README.md](data/README.md). Place `AllWhistlesSubClustering_final.csv` in `data/` (copy or symlink).

## Usage

From the repository root:

```bash
Rscript scripts/01_run_markov_model.R
Rscript scripts/02_node_network_measures.R
Rscript scripts/03_global_network_measures.R
```

Or in R:

```r
setwd("/path/to/MarkovModel")
source("scripts/01_run_markov_model.R")
source("scripts/02_node_network_measures.R")
source("scripts/03_global_network_measures.R")
```

### Script overview

| Script | Purpose | Output |
|--------|---------|--------|
| `01_run_markov_model.R` | Build Markov model (null model + p-values), plot significant-transition network | `outputs/markov_model.rds` |
| `02_node_network_measures.R` | Node metrics, turn-taking probabilities, bar/scatter plots | `outputs/node_measures.csv` |
| `03_global_network_measures.R` | Degree-preserving random networks, clustering p-value, small-world, communities | `outputs/global_null_metrics.rds` |

Script 02 and 03 require `outputs/markov_model.rds` from script 01.

## Parameters

Defaults are in `R/config.R`:

- `ITERATIONS = 1000` — null-model permutations (use `10000` for publication parity with the original monolithic script)
- `TIME_WINDOW = c(0, 5.94)` — seconds after each whistle for counting transitions
- `NULL_MODEL = "shuffle"` — category shuffle null model
- `GRAPH_P_VALUE = 0.05` — edge significance threshold

## Runtime

Script 01 is CPU-intensive (parallel null model). With 1000 iterations on the full dataset (~8500 whistles), expect several minutes depending on core count.

## Original analysis

Refactored from `Markov_Model/Batch_2/MM_batch_2_null_model_not_forced.R`.

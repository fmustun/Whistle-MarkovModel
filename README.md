# Whistle-MarkovModel

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
Rscript scripts/04_bh_corrected_network.R --input=/path/to/markov_model.rds
```

Or in R:

```r
setwd("/path/to/Whistle-MarkovModel")
source("scripts/01_run_markov_model.R")
source("scripts/02_node_network_measures.R")
source("scripts/03_global_network_measures.R")
```

### Script overview

| Script | Purpose | Output |
|--------|---------|--------|
| `01_run_markov_model.R` | Build Markov model (null model + p-values), plot significant-transition network | `outputs/markov_model.rds`, `outputs/plots/markov_significant_network.pdf` |
| `02_node_network_measures.R` | Node metrics, turn-taking probabilities, bar/scatter plots | `outputs/node_measures.csv`, `outputs/plots/*.pdf` |
| `03_global_network_measures.R` | Degree-preserving random networks, clustering p-value, small-world, communities | `outputs/global_null_metrics.rds`, `outputs/plots/*.pdf` |
| `04_bh_corrected_network.R` | Apply valid empirical p-values and BH correction to an existing saved null model | `outputs/bh_analysis/` with p/q matrices, retained edges, fixed coordinates, corrected network RDS, and directly comparable figures |

Figures are written as PDFs under `outputs/plots/` (created automatically). Data files remain in `outputs/`.

Scripts 02 and 03 require `outputs/markov_model.rds` from script 01. Script 04
reuses the null matrices in an existing model; it does not regenerate them.

## Parameters

Defaults are in `R/config.R`:

- `ITERATIONS = 1000` — null-model permutations (use `10000` for publication parity with the original monolithic script)
- `TIME_WINDOW = c(0, 5.94)` — seconds after each whistle for counting transitions
- `NULL_MODEL = "shuffle"` — category shuffle null model
- `GRAPH_P_VALUE = 0.05` — edge significance threshold
- `MULTIPLE_TESTING_CORRECTION = "BH"` — `"BH"` for the revised network or `"none"` for valid uncorrected empirical p-values
- `GRAPH_ALPHA = 0.05` — adjusted-value threshold for script 04

## Benjamini-Hochberg network

Script 04 never generates null models. It reads a saved `markov_model.rds` and
tests the complete family of all ordered state pairs. For the publication model
this is 86 x 86 = 7,396 hypotheses, including auto-loops and zero-observed
transitions. The one-sided enrichment p-value is
`(1 + count(null >= observed)) / (B + 1)`, where ties are exceedances;
zero-observed transitions are assigned p = 1 but remain in the BH family.
An edge is retained only when its observed transition probability is positive
and its adjusted value is below `GRAPH_ALPHA`. Transition probabilities remain
the edge weights.

Run against an existing 150,000-realization model without regenerating it:

```bash
Rscript scripts/04_bh_corrected_network.R \
  --input=/path/to/markov_model.rds \
  --output-dir=/path/to/new/bh_outputs
```

For a checked, previously computed raw empirical-p matrix, add
`--precomputed-p=/path/to/empirical_p_values.rds`. This skips only the
exceedance pass; matrix dimensions, values, zero-observed assignments, BH
counts, and output structure are still validated. Use `--correction=none` to
produce the valid uncorrected network.

The plotting stage computes and exports the original seed-22 85-node layout,
subsets those exact coordinates for BH-retained nodes, uses the original
graph's absolute weight range for edge width and transparency, and represents
significant loops as black node borders. All configured category colours are
preserved.

## Lightweight validation

```bash
Rscript tests/test_statistical_utils.R
Rscript tests/test_bh_correction.R
```

## Runtime

Script 01 is CPU-intensive (parallel null model). With 1000 iterations on the full dataset (~8500 whistles), expect several minutes depending on core count.

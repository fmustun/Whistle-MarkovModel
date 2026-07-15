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
Rscript scripts/04_multiloop_analysis.R
```

Or in R:

```r
setwd("/path/to/Whistle-MarkovModel")
source("scripts/01_run_markov_model.R")
source("scripts/02_node_network_measures.R")
source("scripts/03_global_network_measures.R")
source("scripts/04_multiloop_analysis.R")
```

### Script overview

| Script | Purpose | Output |
|--------|---------|--------|
| `01_run_markov_model.R` | Build Markov model (null model + empirical p-values + BH correction), plot significant-transition network | `outputs/markov_model.rds`, `outputs/plots/markov_significant_network.pdf`, `outputs/bh_analysis/` |
| `02_node_network_measures.R` | Node metrics, turn-taking probabilities, bar/scatter plots | `outputs/node_measures.csv`, `outputs/plots/*.pdf` |
| `03_global_network_measures.R` | Degree-preserving random networks, clustering p-value, small-world, communities | `outputs/global_null_metrics.rds`, `outputs/plots/*.pdf` |
| `04_multiloop_analysis.R` | Group temporally adjacent whistles into multi-loop chains (250ms IWI threshold), summarize the chain-length population and its overlap with the Markov transitions | `outputs/multiloop_analysis/`, `outputs/plots/multiloop_combined_summary.pdf` |

Figures are written as PDFs under `outputs/plots/` (created automatically). Data files remain in `outputs/`.

Scripts 02 and 03 require `outputs/markov_model.rds` from script 01. Script 04 is standalone (only needs the raw whistle CSV).

## Parameters

Defaults are in `R/config.R`:

- `ITERATIONS = 1000` — null-model permutations; a fast development default, not a statistically powered run (see "Benjamini-Hochberg network" below for how to derive a properly powered value)
- `TIME_WINDOW = c(0, 5.94)` — seconds after each whistle for counting transitions
- `NULL_MODEL = "shuffle"` — category shuffle null model
- `GRAPH_P_VALUE = 0.05` — edge significance threshold
- `MULTIPLE_TESTING_CORRECTION = "BH"` — `"BH"` for the revised network or `"none"` for valid uncorrected empirical p-values
- `GRAPH_ALPHA = 0.05` — adjusted-value threshold for edge selection in script 01
- `MULTILOOP_IWI_THRESHOLD = 0.25` — inter-whistle interval (seconds) below which two consecutive whistles are grouped into the same multi-loop chain in script 04

## Benjamini-Hochberg network

`MarkovWhistle()` (called from script 01) computes null-model permutations and,
from them, empirical enrichment p-values directly as part of building the
model. The one-sided enrichment p-value is
`(1 + count(null >= observed)) / (B + 1)`, where ties are exceedances and B is
`ITERATIONS`.

BH correction is applied only to the transitions **observed at least once**
in the data (auto-loops included) — never-observed transitions are assigned
p = 1 and structurally excluded from the correction family, since they can
never pass script 01's `observed_probability > 0` selection requirement
regardless of their adjusted value. Including them would only inflate the
family size and make the correction more conservative for no statistical
benefit ("independent filtering"; Bourgon et al. 2010). This restriction is
implemented via the `family_mask` argument to `adjust_transition_p_values()`
in `R/statistical_utils.R`. Script 01 retains an edge only when its observed
transition probability is positive and its adjusted value is below
`GRAPH_ALPHA`. Transition probabilities remain the edge weights.

**Choosing `ITERATIONS`:** the empirical p-value floor is `1/(ITERATIONS+1)`,
so BH correction across `m` observed transitions at `GRAPH_ALPHA` needs
`ITERATIONS >= m/GRAPH_ALPHA - 1` just to have any resolving power at all —
below that bound nothing can ever be significant, regardless of effect size.
A safety margin above that bare floor keeps p-value estimates near the
decision boundary from being dominated by Monte Carlo noise.
`recommended_iterations(m, alpha, margin = 5)` in `R/statistical_utils.R`
computes this; find `m` (the number of transitions observed at least once)
from a prior run's `outputs/bh_analysis/bh_analysis_summary.txt` or via
`sum(MarkovModel$transition_probabilities_matrix_all > 0)`. Do not choose
`ITERATIONS` by trying values until the retained network "looks right" —
that reintroduces exactly the bias this formula is meant to avoid.

To change the correction method or alpha, edit `MULTIPLE_TESTING_CORRECTION`/
`GRAPH_ALPHA` in `R/config.R` and rerun script 01 (it regenerates the null
model each time).

Alongside `outputs/markov_model.rds` and
`outputs/plots/markov_significant_network.pdf`, script 01 writes diagnostic
outputs to `outputs/bh_analysis/`:

- `empirical_p_values.csv`/`.rds`, `bh_q_values.csv`/`.rds` (or
  `uncorrected_p_values.*` when `MULTIPLE_TESTING_CORRECTION = "none"`),
  `empirical_exceedance_counts.rds` — the full p/q matrices
- `node_layout_coordinates.csv` — per-node id/label/category and the
  seed-22 layout coordinates for the retained network
- `bh_retained_edges.csv` — per-edge transition probability, empirical p,
  adjusted value, and auto-loop flag for every retained transition
- `bh_validation_counts.csv` — retained-node and significant-transition counts
- `bh_corrected_network.rds` — bundles the matrices, selection, graph, and
  layout coordinates together
- `bh_analysis_summary.txt` — plain-text run summary

Script 01 does not pin the retained counts to a specific expected value —
the previous pin (`ITERATIONS = 150000`, expecting `retained_nodes = 66`
etc.) was tied to the old unfiltered hypothesis family and to an iteration
count that had been tuned to reproduce a target result rather than derived
from the resolution requirement above. Once you've run a derived
`ITERATIONS` value once and trust the result, call
`assert_transition_selection_counts(counts, c(...))` (see
`R/statistical_utils.R`) to pin a new regression baseline.

## Multi-loop whistle analysis

Script 04 groups temporally adjacent whistles into "multi-loop" chains,
independent of the Markov transition analysis above: within each recording,
ordered by `start_time`, a new chain starts whenever the gap to the previous
whistle (`start_time[i] - end_time[i-1]`) is `>= MULTILOOP_IWI_THRESHOLD`
(250ms by default). Chaining is purely temporal — chained whistles are not
required to share a type or category. A chain's "leading category" is the
category of the whistle that starts it.

`R/multiloop_utils.R` provides `compute_multiloop_chains()` (whistle-level
chain assignment) and `summarize_multiloop_chains()` (one row per chain).
Script 04 writes `outputs/multiloop_analysis/multiloop_whistles.csv` and
`multiloop_chains.csv`, plus a combined three-panel figure at
`outputs/plots/multiloop_combined_summary.pdf`:

1. Chain-length distribution broken down by leading category (the category
   of the whistle that starts each chain).
2. Multi-loop chains split by whether all whistles in the chain share a
   single main category (e.g. all `SW_Neo`) or mix more than one — "main
   category" is the `LIST_NAMES`-level grouping, not the finer
   `whistle_type_chr` sub-category (e.g. `SW_Neo_Category_20`).
3. The Markov-transition/multi-loop overlap ratios (see below).

A chain of length 1 has no neighbor within the IWI threshold, so it isn't a
multi-loop by definition. The CSVs retain these singleton rows (with
`chain_length = 1`) for completeness, and the summary reports their count/
percentage for context, but they are excluded from every multi-loop measure
(mean/median/max chain length, chain-length counts, category composition)
and from the figure.

### Overlap with the Markov transition pairs

Script 04 also checks how much of the Markov model's transition pairs
(the ones `MarkovWhistle()`/script 01 actually counts into
`transition_counts_matrix`, using the much wider `TIME_WINDOW`) connect two
whistles that happen to be in the same tight-IWI multi-loop chain.
`compute_markov_transition_multiloop_overlap()` (`R/multiloop_utils.R`)
replicates script 01's empirical pair-selection logic (not the null model)
instance-by-instance and flags each transition as `within_multiloop_chain`
or not; its instance count matches `sum(transition_counts_matrix)` from a
prior script 01 run exactly, confirming it selects the same pairs. Note that
an overlapping pair (negative IWI gap) is itself never recorded as a Markov
transition — script 01 skips straight to the whistle after it when one is
available in the window — so some multi-loop-adjacent whistle pairs never
appear as a transition at all. Instance-level results are written to
`outputs/multiloop_analysis/markov_transition_multiloop_overlap.csv`; the
overall percentage plus a same-type vs. different-type breakdown appear in
`multiloop_summary.txt` and in the third panel of the combined figure.

## Lightweight validation

```bash
Rscript tests/test_statistical_utils.R
Rscript tests/test_bh_correction.R
Rscript tests/test_multiloop_utils.R
```

## Runtime

Script 01 is CPU-intensive (parallel null model). With 1000 iterations on the full dataset (~8500 whistles), expect several minutes depending on core count.

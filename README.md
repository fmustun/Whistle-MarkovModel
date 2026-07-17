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

Script 05 (multi-loop spectrograms) is Python instead of R, and additionally
needs: `numpy`, `scipy`, `matplotlib`, `soundfile`.

## Data

See [data/README.md](data/README.md). Place `AllWhistlesSubClustering_final.csv` in `data/` (copy or symlink).

## Usage

From the repository root:

```bash
Rscript scripts/01_run_markov_model.R
Rscript scripts/02_node_network_measures.R
Rscript scripts/03_global_network_measures.R
Rscript scripts/04_multiloop_analysis.R
python3 scripts/05_plot_multiloop_spectrograms.py
Rscript scripts/06_plot_multiloop_network.R
Rscript scripts/07_plot_multiloop_type_network.R
Rscript scripts/08_plot_multiloop_transitions_network.R
Rscript scripts/09_plot_interevent_network.R
```

Or in R:

```r
setwd("/path/to/Whistle-MarkovModel")
source("scripts/01_run_markov_model.R")
source("scripts/02_node_network_measures.R")
source("scripts/03_global_network_measures.R")
source("scripts/04_multiloop_analysis.R")
source("scripts/06_plot_multiloop_network.R")
source("scripts/07_plot_multiloop_type_network.R")
source("scripts/08_plot_multiloop_transitions_network.R")
source("scripts/09_plot_interevent_network.R")
```

### Script overview

| Script | Purpose | Output |
|--------|---------|--------|
| `01_run_markov_model.R` | Build Markov model (null model + empirical p-values + BH correction), plot significant-transition network | `outputs/markov_model.rds`, `outputs/plots/markov_significant_network.pdf`, `outputs/bh_analysis/` |
| `02_node_network_measures.R` | Node metrics, turn-taking probabilities, bar/scatter plots | `outputs/node_measures.csv`, `outputs/plots/network_measures/*.pdf` |
| `03_global_network_measures.R` | Degree-preserving random networks, clustering p-value, small-world, communities | `outputs/global_null_metrics.rds`, `outputs/plots/*.pdf`, `outputs/plots/network_measures/*.pdf` |
| `04_multiloop_analysis.R` | Group temporally adjacent whistles into multi-loop chains (250ms IWI threshold), summarize the chain-length population and its overlap with the Markov transitions | `outputs/multiloop_analysis/`, `outputs/plots/multiloops/multiloop_combined_summary.pdf` |
| `05_plot_multiloop_spectrograms.py` | Plot grayscale spectrograms (2-20kHz) of example multi-loop chains from the source audio, one figure per chain with every whistle in the chain marked | `outputs/multiloop_analysis/spectrograms/*.png` |
| `06_plot_multiloop_network.R` | Plot a multi-loop-only network: script 01's node layout, but edges are multi-loop-internal transitions only, black border = same-type multi-loop repeat | `outputs/plots/multiloops/multiloop_network_overlay.pdf`, `outputs/multiloop_analysis/{multiloop_edge_counts,node_multiloop_participation}.csv` |
| `07_plot_multiloop_type_network.R` | Build a Markov-style transition network where each node is a distinct multi-loop composition (not an individual whistle); BH-corrected empirical shuffle-null significance test, only significant transitions plotted | `outputs/plots/multiloops/multiloop_type_network.pdf`, `outputs/multiloop_analysis/{multiloop_type_nodes,multiloop_type_edges}.csv` |
| `08_plot_multiloop_transitions_network.R` | Rebuild script 01's whistle-level network (same node universe, same BH-corrected significance test) with transition counts restricted to instances where both whistles ARE in the same multi-loop chain | `outputs/plots/multiloops/multiloop_only_significant_network.pdf`, `outputs/multiloop_analysis/{multiloop_only_network_nodes,multiloop_only_network_edges}.csv` |
| `09_plot_interevent_network.R` | The complement of script 08: same rebuild, but transition counts are restricted to instances where the two whistles are NOT in the same multi-loop chain | `outputs/plots/multiloops/interevent_significant_network.pdf`, `outputs/multiloop_analysis/{interevent_network_nodes,interevent_network_edges}.csv` |

Script 05 requires script 04's `outputs/multiloop_analysis/{multiloop_chains,multiloop_whistles}.csv` and read access to the raw `.wav` recordings (`--audio-dir`, default `/media/zfnews31/Dolphins/Sound`).

Plots are organized under `outputs/plots/`: multi-loop-related figures (scripts
04, 06-09) go in `outputs/plots/multiloops/`, per-node/global network-measure
figures (scripts 02-03) go in `outputs/plots/network_measures/`, and script
01's main network plot (plus script 03's community-detection plots) stay at
the top level. `plot_path()`/`ensure_plots_dir()`/`save_ggplot()`
(`R/plot_io.R`) all take an optional `subdir` argument for this.

Script 06 requires script 01's `outputs/markov_model.rds` (it computes multi-loop chains itself from the raw whistle CSV, like script 04).

Figures are written as PDFs under `outputs/plots/` (created automatically). Data files remain in `outputs/`.

Scripts 02, 03, and 06 require `outputs/markov_model.rds` from script 01. Scripts 04, 07, 08, and 09 are standalone (only need the raw whistle CSV), as is 06 aside from its script 01 dependency.

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
`multiloop_chains.csv`, plus a combined four-panel figure at
`outputs/plots/multiloops/multiloop_combined_summary.pdf`:

1. Chain-length distribution broken down by leading category (the category
   of the whistle that starts each chain).
2. Multi-loop chains split by whether all whistles in the chain share a
   single main category (e.g. all `SW_Neo`) or mix more than one — "main
   category" is the `LIST_NAMES`-level grouping, not the finer
   `whistle_type_chr` sub-category (e.g. `SW_Neo_Category_20`).
3. The Markov-transition/multi-loop overlap ratios (see below).
4. Sub-category pattern recurrence (see below).

A chain of length 1 has no neighbor within the IWI threshold, so it isn't a
multi-loop by definition. The CSVs retain these singleton rows (with
`chain_length = 1`) for completeness, and the summary reports their count/
percentage for context, but they are excluded from every multi-loop measure
(mean/median/max chain length, chain-length counts, category composition,
pattern recurrence) and from the figure.

### Do specific multi-loops recur?

Script 04 checks whether particular multi-loops show up more than once,
under two definitions of "the same multi-loop":

- **exact sequence**: the same sub-categories (`whistle_type_chr`) in the
  same order/position, e.g. two chains that are both
  `A -> A -> B`.
- **same composition, any order**: the same sub-categories regardless of
  position, e.g. `A -> A -> B` and `A -> B -> A` count as the same pattern
  here (but not under "exact sequence").

For this dataset, most patterns are unique, but there's a real recurring
one: `SW_Luna_Category_17 -> SW_Luna_Category_17` (a repeated 2-loop of the
same sub-category) appears in 73 different chains across 41 recordings —
by far the most common recurring multi-loop, under either definition. The
full per-pattern breakdown (pattern, occurrence count, and the specific
`chain_id`s/recordings it appears in) is written to
`outputs/multiloop_analysis/multiloop_pattern_occurrences.csv`; the top 5
recurring patterns under each definition are listed in
`multiloop_summary.txt`, and the fourth panel of the combined figure plots
the full distribution — number of chains sharing a pattern (x) vs. number
of distinct patterns with that many chains (y) — for both definitions
side by side.

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

### Example spectrograms

`scripts/05_plot_multiloop_spectrograms.py` plots a grayscale spectrogram
(2-20kHz, `scipy.signal.spectrogram` on a short audio segment read directly
from the recording's `.wav` file) for a small set of example multi-loop
chains, with every whistle in the chain drawn on the same figure (shaded
span + label) so the repeated/looped structure is visible in one image. By
default it picks up to 6 examples spread across the observed chain-length
range (one per length, using the longest-duration chain at that length) so
both common short loops and rare long ones are represented — pass
`--chain-ids <chain_id> ...` to plot specific chains instead, or
`--n-examples` to change the count. Output goes to
`outputs/multiloop_analysis/spectrograms/<chain_id>.png`. Audio files are
matched to the `recording` column by file stem (recursively under
`--audio-dir`); only the needed segment is read from disk.

### Multi-loop-only network

`scripts/06_plot_multiloop_network.R` draws a separate network that reuses
script 01's node set and layout (recomputed live from `gra1` via
`compute_markov_layout(gra1, PLOT_SEED)` — the same call
`markov_significant_network.pdf` uses internally — rather than a saved
coordinates file: `layout_with_fr` is force-directed, so its result depends
on the edge set it's run on, and this script's edges are about to be
replaced) but replaces its edges entirely:

- **Edges** are adjacent *different*-whistle-type pairs observed within
  multi-loop chains (IWI < `MULTILOOP_IWI_THRESHOLD`) —
  `compute_multiloop_edge_counts()` (`R/multiloop_utils.R`) — not the
  Markov model's wide-`TIME_WINDOW` transitions. Edge width is the count of
  that pair within multi-loop chains. A multi-loop pair whose whistle type
  falls outside script 01's significant-node set is dropped (reported via
  `message()`) rather than extending the network with new nodes.
- **Node border**: black if that whistle has a same-type multi-loop repeat
  (e.g. `A -> A` within a chain, `compute_self_repeat_counts()`), gray
  otherwise. This replaces the node's usual black/gray meaning from script
  01 (Markov self-transition significance) with a purely multi-loop-based
  one, since this plot no longer shows Markov edges at all.

Node fill color (category) is unchanged from script 01. Output goes to
`outputs/plots/multiloops/multiloop_network_overlay.pdf`, with the drawn edges and
per-node self-repeat flags written to
`outputs/multiloop_analysis/multiloop_edge_counts.csv` and
`node_multiloop_participation.csv`. `plot_markov_graph()`
(`R/graph_utils.R`) gained one optional `caption` parameter (default
`NULL`, unused by scripts 01/01b) to label the plot's now-different edge
and border meanings.

### Multi-loop *type* network

`scripts/07_plot_multiloop_type_network.R` builds a Markov-style network
the same way script 01 does (transition counts → row-normalized
probabilities via `ccber::CalcTransitionMatrix()` → directed graph, tested
against a null model, BH-corrected, only significant edges plotted), but at
a different unit: each **node is a distinct multi-loop composition** —
`compute_multiloop_patterns()`'s `unordered_pattern`, i.e. the same
sub-categories regardless of order/position (630 distinct compositions in
this dataset, from `SW_Luna_Category_17 + SW_Luna_Category_17` — the most
common, 73 chains — down to singletons). A candidate **edge** requires the
next multi-loop chain in the same recording to start within `TIME_WINDOW`
of the previous chain's end (same gap definition as `MarkovWhistle()`'s
`time_interval`) — 310 of the 967 consecutive same-recording chain pairs in
this dataset qualify. Multi-loop chains never overlap in time by
construction (a new chain requires a >= `MULTILOOP_IWI_THRESHOLD` gap), so
unlike `MarkovWhistle()` there's no need for its overlap/skip-to-next-candidate
handling — just a direct window check.

**Significance test.** `compute_multiloop_type_significance()`
(`R/multiloop_utils.R`) runs the same empirical-p / BH-correction test
script 01's `MarkovWhistle()` does, against the same `NULL_MODEL = "shuffle"`
null model (a global permutation of which composition label sits on which
multi-loop chain, keeping every chain's timing — and therefore which
chain-pairs fall within `TIME_WINDOW` — fixed), then
`adjust_transition_p_values()`/`select_transition_matrix()`
(`R/statistical_utils.R`, unchanged) apply the same `MULTIPLE_TESTING_CORRECTION`
and `GRAPH_ALPHA` as script 01. Only **significant** transitions
(`adjusted p < GRAPH_ALPHA`) are drawn as edges. The iteration count is
derived the same way script 01's `ITERATIONS` is documented to be derived —
`recommended_iterations(m, GRAPH_ALPHA, margin = 5)` where `m` is the
number of observed edges (~300 here, vs. ~1,200 at the whistle level) — so
it scales with this network's own hypothesis-family size rather than
reusing script 01's `ITERATIONS`. This network has far more nodes (630)
relative to observed edges (~300) than the whistle-level network does, so
building a dense 630×630 null-count matrix per iteration (as
`MarkovWhistle()` does at the whistle level) would be mostly zeros and,
at tens of thousands of iterations, too much memory/compute to hold;
`compute_multiloop_type_significance()` instead computes the empirical
p-value only at the observed-edge cells that can ever be significant
(`adjust_transition_p_values()`'s `family_mask` already restricts the
correction family to exactly those cells, so nothing downstream is
affected by skipping the rest).

Node color is the composition's single main category (same palette as the
whistle-type network) if every whistle in it shares one, or a fixed mid
gray if the composition spans more than one main category
(`n_categories > 1`, e.g. an `SW_Neo` whistle immediately followed by an
`SW_Luna` one within the same chain) — 361 single-category vs. 269
mixed-category compositions in this dataset (white was tried first but read
as too close to `NSW_9`'s very pale lavender, `#e6e6fa`). A **black node
border** follows script 01's convention exactly: it marks a composition
with a *significant* self-transition (a nonzero diagonal entry in the
selected/significant matrix) — drawn as a border rather than a self-loop
arc, which is then dropped from the edge set entirely (`delete_edges(...,
which_loop(...))`, same as `compute_graph_from_selection()` does for script
01's significant self-loops). The isolate filter (drop nodes with no
significant transition at all) runs *before* that self-loop removal, again
matching `compute_graph()`'s ordering — so a composition whose only
significant transition is a self-loop still appears in the plot, as an
isolated black-bordered dot with no edges, rather than being dropped (311
of 630 remain; the full node set is still in `multiloop_type_nodes.csv`).
This network gets its own independent layout (`compute_markov_layout()`,
not tied to script 01's node positions, since its nodes are entirely
different things). Output goes to `outputs/plots/multiloops/multiloop_type_network.pdf`,
with node/edge details in `outputs/multiloop_analysis/multiloop_type_nodes.csv`
and `multiloop_type_edges.csv` (the latter includes every observed edge,
not just the significant ones, with `empirical_p`/`bh_q`/`significant`
columns alongside `is_self_transition`).

### Multi-loop-only whistle-level network

`scripts/08_plot_multiloop_transitions_network.R` is a third network,
alongside script 01's full one and script 06's overlay: a **rebuild of
script 01's whole pipeline** (transition counts → row-normalized
probabilities → BH-corrected empirical shuffle-null test →
`select_transition_matrix()` → `compute_graph_from_selection()`) on the
**same node universe** — every whistle sub-category `MarkovWhistle()`
would consider (86 in this dataset), not just script 01's already-pruned
82-node subset — but with the transition counts feeding that test
restricted to Markov transition instances where **both whistles belong to
the same multi-loop chain**
(`compute_markov_transition_multiloop_overlap()`'s `within_multiloop_chain`
flag, IWI < `MULTILOOP_IWI_THRESHOLD`) — 970 of 4,671 instances (20.8%) in
this dataset qualify, the same figure script 04 reports. This differs from
script 06, which overlays raw (untested) multi-loop counts onto script 01's
*existing* significant-network node positions; script 08 instead asks
"what does the significant network look like if you only feed it
multi-loop-internal transitions?", with its own independent
`compute_markov_layout()` layout and its own significance test.

The significance test reuses the same permutation machinery as script 07
(`compute_permutation_transition_significance()`, `R/multiloop_utils.R`) —
refactored out of `compute_multiloop_type_significance()` once script 08
needed the identical logic at a different granularity (individual whistles
here, not multi-loop chains) — under `NULL_MODEL = "shuffle"`: a global
permutation of `whistle_type` across every whistle, keeping timing fixed.
Since multi-loop chain membership is purely a function of gaps between
start/end times, `within_multiloop_chain` (and therefore which transition
instances even count) is unaffected by that permutation, so one null
iteration only needs to permute the label vector and re-look-up the fixed
set of transition-instance row pairs
(`compute_markov_transition_multiloop_overlap()`'s new `from_row_id`/
`to_row_id` columns — each instance's row position in the input whistle
table, needed to look up a whistle's label under an arbitrary permutation).
The iteration count is derived the same way as script 07's, from this
network's own smaller hypothesis-family size (`m` ≈ 385 observed edges
here, vs. ~1,200 at the full whistle level) via
`recommended_iterations(m, GRAPH_ALPHA, margin = 5)`. `set.seed(SEED)` is
called immediately before the null model, matching script 01, so results
are reproducible run to run.

Node/edge construction, category colors, and the black-border
significant-self-transition convention are otherwise identical to script
01's (same `compute_graph()`/`compute_graph_from_selection()` calls) — 76
of 86 sub-categories retain at least one significant multi-loop-only
transition. Output goes to
`outputs/plots/multiloops/multiloop_only_significant_network.pdf`, with
node/edge details (including layout coordinates and every observed edge
with `empirical_p`/`bh_q`/`significant` columns, not just the significant
ones) in `outputs/multiloop_analysis/multiloop_only_network_nodes.csv` and
`multiloop_only_network_edges.csv`.

Filtering down to multi-loop-only edges frequently leaves the network not
fully connected — e.g. two sub-categories that only ever transition to
each other, nowhere near the main network in this dataset (76 nodes: one
74-node component plus an isolated pair). Forcing `compute_markov_layout()`'s
force-directed layout to place every component in one shared coordinate
space stretches the whole plot to make room for these outliers, at the
cost of legibility for the (much larger) main component. Instead,
`plot_network_by_component()` (`R/graph_utils.R`) splits the graph by
connected component (`igraph::components()`, weak/undirected sense) and
gives each its own independent layout and its own panel in the same
figure — the largest component as the main panel, every smaller component
stacked below it in its own (smaller) panel via `layout()` — so nothing is
dropped from the plot, but nothing distorts the main network either. It
returns the combined node table (with `component`/`component_size`
columns recording which panel each node ended up in) ready to write to
CSV, so a caller just needs a significant-transitions graph and a caption;
scripts 08 and 09 (below) both call it directly rather than duplicating
the layout/panel logic.

### Inter-event network

`scripts/09_plot_interevent_network.R` is script 08's exact complement:
same rebuild-of-script-01 pipeline, same node universe, same significance
test, but transition counts are restricted to instances where the two
whistles **do NOT** belong to the same multi-loop chain — "inter-event"
transitions, in the sense of crossing between separate multi-loop chains
(or between a chain and a singleton) rather than staying within one. Every
Markov transition instance is within a multi-loop chain (script 08) or is
not (script 09) — the two scripts' edge lists exactly partition
`compute_markov_transition_multiloop_overlap()`'s full output, and never
overlap. In this dataset 3,701 of 4,671 instances (79.2%) qualify — the
complement of script 08's 20.8% — giving a much larger hypothesis family
(`m` ≈ 1,065 vs. ~385) and a correspondingly larger derived iteration
count. 84 of 86 sub-categories retain a significant inter-event
transition, and (unlike script 08) the result is a single connected
component in this dataset, so `plot_network_by_component()` renders it as
one panel. Output goes to
`outputs/plots/multiloops/interevent_significant_network.pdf`, with
node/edge details in
`outputs/multiloop_analysis/interevent_network_nodes.csv` and
`interevent_network_edges.csv`.

### Betweenness vs. strength (multi-loop-only / inter-event)

Both scripts 08 and 09 also call `plot_betweenness_strength_scatter()`
(`R/graph_utils.R`) on their respective significant graph, producing
`outputs/plots/multiloops/multiloop_only_betweenness_vs_strength.pdf` and
`interevent_betweenness_vs_strength.pdf` — a scatter of each node's
betweenness against its total strength, each normalized to `[0, 1]` by
dividing by its own max within that network, with a marginal histogram
along each axis (plus a fitted exponential density curve, `rate = 1 /
mean` — the simplest single-parameter model for the right-skewed shape
these centrality measures typically have). Red dashed lines mark each
measure's empirical 90th percentile; a node is labeled with its id when it
exceeds *either* line (an outlier on betweenness, strength, or both) —
this is a different (and stricter) rule than script 02's existing
`betweenness_vs_strength_in` plot, which labels every node and uses
in-strength rather than total strength. This is a generic, reusable
function — not tied to multi-loop analysis — so it could equally be
pointed at `gra1` or any other significant-transitions graph.

## Lightweight validation

```bash
Rscript tests/test_statistical_utils.R
Rscript tests/test_bh_correction.R
Rscript tests/test_multiloop_utils.R
```

## Runtime

Script 01 is CPU-intensive (parallel null model). With 1000 iterations on the full dataset (~8500 whistles), expect several minutes depending on core count.

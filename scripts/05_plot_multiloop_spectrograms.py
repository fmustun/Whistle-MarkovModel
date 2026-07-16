#!/usr/bin/env python3
"""Plot spectrograms for a handful of example multi-loop whistle chains.

Reads outputs/multiloop_analysis/{multiloop_chains,multiloop_whistles}.csv
(written by scripts/04_multiloop_analysis.R) and, for a small set of example
multi-loop chains (chain_length >= 2), plots a single grayscale spectrogram
spanning all of that chain's whistles, with each whistle's time span marked.

Audio is read directly from the recording's .wav file under --audio-dir,
matched by the "recording" column (file name minus extension). Only the
needed segment is read (soundfile seeks rather than loading the full file).
"""
import argparse
import csv
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import soundfile as sf
from scipy.signal import spectrogram

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_AUDIO_DIR = Path("/media/zfnews31/Dolphins/Sound")
DEFAULT_CHAINS_CSV = (
    REPO_ROOT / "outputs" / "multiloop_analysis" / "multiloop_chains.csv"
)
DEFAULT_WHISTLES_CSV = (
    REPO_ROOT / "outputs" / "multiloop_analysis" / "multiloop_whistles.csv"
)
DEFAULT_OUTPUT_DIR = REPO_ROOT / "outputs" / "multiloop_analysis" / "spectrograms"

FREQ_MIN_HZ = 2000
FREQ_MAX_HZ = 20000
PADDING_SECONDS = 0.5
NPERSEG = 1024


def read_csv_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def build_audio_index(audio_dir):
    """Maps recording name (file stem) -> absolute .wav path, recursively."""
    index = {}
    for root, _dirs, files in os.walk(audio_dir):
        for fn in files:
            if fn.lower().endswith(".wav"):
                index.setdefault(os.path.splitext(fn)[0], os.path.join(root, fn))
    return index


def select_example_chains(chains, n):
    """Picks up to n example chains.

    First round: one chain per distinct chain length (the longest-duration
    chain at that length), spread evenly across the observed length range
    via linspace when n is smaller than the number of distinct lengths - so
    even a handful of examples covers both short and long multi-loops
    rather than clustering in the common (short) end.

    If n exceeds the number of distinct lengths, every length gets covered
    in the first round and additional rounds round-robin through all
    lengths again (ascending, skipping exhausted ones) for a 2nd/3rd/...
    example per length, up to the total number of multi-loop chains
    available."""
    multiloop = [c for c in chains if int(c["chain_length"]) >= 2]
    by_length = {}
    for c in multiloop:
        by_length.setdefault(int(c["chain_length"]), []).append(c)
    for candidates in by_length.values():
        candidates.sort(
            key=lambda c: float(c["end_time"]) - float(c["start_time"]), reverse=True
        )
    lengths_present = sorted(by_length)

    if n >= len(lengths_present):
        first_round_lengths = lengths_present
    else:
        idx = np.linspace(0, len(lengths_present) - 1, n)
        first_round_lengths = sorted({lengths_present[int(round(i))] for i in idx})

    selected = []
    cursor = {length: 0 for length in lengths_present}
    for length in first_round_lengths:
        selected.append(by_length[length][0])
        cursor[length] = 1

    while len(selected) < n:
        added = False
        for length in lengths_present:
            if len(selected) >= n:
                break
            i = cursor[length]
            candidates = by_length[length]
            if i < len(candidates):
                selected.append(candidates[i])
                cursor[length] = i + 1
                added = True
        if not added:
            break  # every multi-loop chain has been selected
    return selected


def plot_chain_spectrogram(chain, whistles_by_chain, audio_index, output_dir):
    chain_id = chain["chain_id"]
    recording = chain["recording"]
    whistles = sorted(
        whistles_by_chain[chain_id], key=lambda w: float(w["start_time"])
    )

    audio_path = audio_index.get(recording)
    if audio_path is None:
        raise FileNotFoundError(f"No .wav file found for recording '{recording}'")
    info = sf.info(audio_path)

    seg_start = max(0.0, float(whistles[0]["start_time"]) - PADDING_SECONDS)
    seg_end = min(info.duration, float(whistles[-1]["end_time"]) + PADDING_SECONDS)

    start_frame = int(seg_start * info.samplerate)
    n_frames = int((seg_end - seg_start) * info.samplerate)
    audio, sr = sf.read(
        audio_path, start=start_frame, frames=n_frames, dtype="float32",
        always_2d=False,
    )
    if audio.ndim > 1:
        audio = audio[:, 0]

    noverlap = int(NPERSEG * 0.75)
    freqs, times, sxx = spectrogram(
        audio, fs=sr, window="hann", nperseg=NPERSEG, noverlap=noverlap
    )
    freq_mask = (freqs >= FREQ_MIN_HZ) & (freqs <= FREQ_MAX_HZ)
    freqs = freqs[freq_mask]
    sxx = sxx[freq_mask, :]
    sxx_db = 10 * np.log10(sxx + 1e-12)
    vmin, vmax = np.percentile(sxx_db, [5, 99.5])

    fig_width = float(np.clip(1.2 * (seg_end - seg_start), 6, 20))
    fig, ax = plt.subplots(figsize=(fig_width, 4.5))
    mesh = ax.pcolormesh(
        times, freqs / 1000.0, sxx_db, cmap="Greys", vmin=vmin, vmax=vmax,
        shading="auto",
    )

    for whistle in whistles:
        w_start = float(whistle["start_time"]) - seg_start
        w_end = float(whistle["end_time"]) - seg_start
        ax.axvspan(w_start, w_end, color="red", alpha=0.12, lw=0)
        ax.axvline(w_start, color="red", alpha=0.4, lw=0.7, linestyle="--")
        ax.text(
            (w_start + w_end) / 2, FREQ_MAX_HZ / 1000.0 - 0.3,
            whistle["whistle_type_chr"], rotation=90, ha="center", va="top",
            fontsize=6, color="red",
        )

    ax.set_xlim(0, seg_end - seg_start)
    ax.set_ylim(FREQ_MIN_HZ / 1000.0, FREQ_MAX_HZ / 1000.0)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("frequency (kHz)")
    ax.set_title(
        f"{recording} | {chain_id} | {len(whistles)}-loop chain "
        f"({float(whistles[0]['start_time']):.2f}s - "
        f"{float(whistles[-1]['end_time']):.2f}s)"
    )
    fig.colorbar(mesh, ax=ax, label="power (dB)")
    fig.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{chain_id}.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--n-examples", type=int, default=6,
        help=(
            "Number of example multi-loop chains to plot (default: 6). "
            "Capped by the total number of multi-loop chains available."
        ),
    )
    parser.add_argument(
        "--chain-ids", nargs="*", default=None,
        help="Explicit chain_id(s) to plot instead of auto-selecting examples",
    )
    parser.add_argument("--audio-dir", default=str(DEFAULT_AUDIO_DIR))
    parser.add_argument("--chains-csv", default=str(DEFAULT_CHAINS_CSV))
    parser.add_argument("--whistles-csv", default=str(DEFAULT_WHISTLES_CSV))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    args = parser.parse_args()

    chains = read_csv_rows(args.chains_csv)
    whistles = read_csv_rows(args.whistles_csv)

    whistles_by_chain = {}
    for w in whistles:
        whistles_by_chain.setdefault(w["chain_id"], []).append(w)

    if args.chain_ids:
        chains_by_id = {c["chain_id"]: c for c in chains}
        selected = [chains_by_id[cid] for cid in args.chain_ids]
    else:
        selected = select_example_chains(chains, args.n_examples)

    print(f"Indexing .wav files under {args.audio_dir} ...")
    audio_index = build_audio_index(args.audio_dir)

    print(f"Plotting {len(selected)} example multi-loop chain(s)...")
    output_dir = Path(args.output_dir)
    for chain in selected:
        out_path = plot_chain_spectrogram(
            chain, whistles_by_chain, audio_index, output_dir
        )
        print(f"  {chain['chain_id']} ({chain['chain_length']}-loop): {out_path}")


if __name__ == "__main__":
    main()

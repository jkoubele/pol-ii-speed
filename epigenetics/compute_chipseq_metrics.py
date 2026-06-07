#!/usr/bin/env python3
"""Compute per-gene ChIP-seq metrics from bigWig and narrowPeak files.

Walks chipseq_data/{degron}/{mark}/{control,degron}/ and for each (degron, mark)
pair produces a TSV with signal and peak metrics over promoters and gene bodies.

Usage (from project root):
    python epigenetics/compute_chipseq_metrics.py
    python epigenetics/compute_chipseq_metrics.py --chipseq-dir /path/to/chipseq_data
"""

import argparse
import gzip
import logging
import numpy as np
import pandas as pd
import pyBigWig
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger(__name__)

# Added to both numerator and denominator before log2FC to avoid log(0).
# Fold-change-over-control bigWigs have values ~1-20 for enriched regions,
# so 1e-3 is well below the noise floor.
PSEUDOCOUNT = 1e-3


def load_bed(path):
    """Load BED file; add chr prefix to match bigWig/peak chromosome names."""
    df = pd.read_csv(
        path, sep="\t", header=None,
        names=["chrom", "start", "end", "gene_id", "score", "strand"],
        usecols=[0, 1, 2, 3, 4, 5],
        dtype={"chrom": str},
    )
    df["chrom"] = df["chrom"].apply(lambda c: c if c.startswith("chr") else "chr" + c)
    df = df.drop_duplicates(subset="gene_id", keep="first")
    return df.set_index("gene_id")


def load_peaks(path):
    """Load narrowPeak (optionally gzipped) into a DataFrame grouped by chrom."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        df = pd.read_csv(
            f, sep="\t", header=None,
            names=["chrom", "start", "end", "name", "score",
                   "strand", "signal_value", "p_value", "q_value", "summit"],
        )
    return {chrom: grp.reset_index(drop=True) for chrom, grp in df.groupby("chrom")}


def bw_stats(bw, chrom, start, end):
    """Return (mean, max) fold-change signal over region. Returns (0, 0) for missing/empty."""
    if chrom not in bw.chroms() or end <= start:
        return 0.0, 0.0
    try:
        mean = bw.stats(chrom, start, end, type="mean", nBins=1)[0]
        maxi = bw.stats(chrom, start, end, type="max", nBins=1)[0]
        return (mean or 0.0, maxi or 0.0)
    except RuntimeError:
        return 0.0, 0.0


def overlapping_peaks(peaks_by_chrom, chrom, start, end):
    """Return rows from peaks_by_chrom that overlap [start, end)."""
    chrom_peaks = peaks_by_chrom.get(chrom)
    if chrom_peaks is None or chrom_peaks.empty:
        return chrom_peaks if chrom_peaks is not None else pd.DataFrame()
    mask = (chrom_peaks["start"] < end) & (chrom_peaks["end"] > start)
    return chrom_peaks[mask]


def region_metrics(prefix, bw_ctrl, bw_deg, peaks_ctrl, peaks_deg, chrom, start, end):
    """Compute all metrics for one genomic region; returns a flat dict."""
    ctrl_mean, ctrl_max = bw_stats(bw_ctrl, chrom, start, end)
    deg_mean, deg_max = bw_stats(bw_deg, chrom, start, end)

    lfc_mean = np.log2((deg_mean + PSEUDOCOUNT) / (ctrl_mean + PSEUDOCOUNT))
    lfc_max = np.log2((deg_max + PSEUDOCOUNT) / (ctrl_max + PSEUDOCOUNT))

    ov_ctrl = overlapping_peaks(peaks_ctrl, chrom, start, end)
    ov_deg = overlapping_peaks(peaks_deg, chrom, start, end)

    def peak_stats(ov, tag):
        n = len(ov)
        sv = ov["signal_value"] if n > 0 else pd.Series(dtype=float)
        return {
            f"{prefix}_{tag}_n_peaks": n,
            f"{prefix}_{tag}_has_peak": int(n > 0),
            f"{prefix}_{tag}_max_signal": sv.max() if n > 0 else 0.0,
            f"{prefix}_{tag}_sum_signal": sv.sum() if n > 0 else 0.0,
            f"{prefix}_{tag}_mean_signal": sv.mean() if n > 0 else 0.0,
            f"{prefix}_{tag}_max_qvalue": ov["q_value"].max() if n > 0 else 0.0,
        }

    result = {
        f"{prefix}_ctrl_mean_fc": ctrl_mean,
        f"{prefix}_ctrl_max_fc": ctrl_max,
        f"{prefix}_deg_mean_fc": deg_mean,
        f"{prefix}_deg_max_fc": deg_max,
        f"{prefix}_lfc_mean_fc": lfc_mean,
        f"{prefix}_lfc_max_fc": lfc_max,
    }
    result.update(peak_stats(ov_ctrl, "ctrl"))
    result.update(peak_stats(ov_deg, "deg"))
    result[f"{prefix}_delta_n_peaks"] = len(ov_deg) - len(ov_ctrl)
    return result


def process_mark(degron, mark, mark_dir, promoters, gene_bodies, out_dir):
    ctrl_dir = mark_dir / "control"
    deg_dir = mark_dir / "degron"

    required = [
        ctrl_dir / "fold_change_over_control.bigWig",
        deg_dir / "fold_change_over_control.bigWig",
        ctrl_dir / "peaks.bed.gz",
        deg_dir / "peaks.bed.gz",
    ]
    missing = [p for p in required if not p.exists()]
    if missing:
        log.warning(f"Skipping {degron}/{mark}: missing {[p.name for p in missing]}")
        return

    log.info(f"Processing {degron}/{mark}")

    bw_ctrl = pyBigWig.open(str(ctrl_dir / "fold_change_over_control.bigWig"))
    bw_deg = pyBigWig.open(str(deg_dir / "fold_change_over_control.bigWig"))
    peaks_ctrl = load_peaks(ctrl_dir / "peaks.bed.gz")
    peaks_deg = load_peaks(deg_dir / "peaks.bed.gz")

    all_genes = promoters.index.union(gene_bodies.index)
    records = []

    for gene_id in all_genes:
        rec = {"gene_id": gene_id}

        if gene_id in promoters.index:
            row = promoters.loc[gene_id]
            rec.update(region_metrics(
                "promoter", bw_ctrl, bw_deg, peaks_ctrl, peaks_deg,
                row["chrom"], int(row["start"]), int(row["end"]),
            ))

        if gene_id in gene_bodies.index:
            row = gene_bodies.loc[gene_id]
            rec.update(region_metrics(
                "body", bw_ctrl, bw_deg, peaks_ctrl, peaks_deg,
                row["chrom"], int(row["start"]), int(row["end"]),
            ))

        records.append(rec)

    bw_ctrl.close()
    bw_deg.close()

    out_path = out_dir / degron / mark / "chipseq_metrics.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(out_path, sep="\t", index=False, float_format="%.6g")
    log.info(f"  -> {out_path} ({len(records)} genes)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chipseq-dir", default="epigenetics/chipseq_data",
                        help="Root directory containing {degron}/{mark}/{condition}/ subdirs")
    parser.add_argument("--intervals-dir", default="epigenetics/genomic_intervals",
                        help="Directory with .bed and gene_bodies.bed")
    parser.add_argument("--out-dir", default="epigenetics/epigenetic_marks",
                        help="Output root; TSVs written to {out-dir}/{degron}/{mark}/chipseq_metrics.tsv")
    args = parser.parse_args()

    promoters = load_bed(Path(args.intervals_dir) / "promoters.bed")
    gene_bodies = load_bed(Path(args.intervals_dir) / "gene_bodies.bed")
    log.info(f"Loaded {len(promoters)} promoters, {len(gene_bodies)} gene bodies")

    chipseq_dir = Path(args.chipseq_dir)
    out_dir = Path(args.out_dir)

    for degron_dir in sorted(chipseq_dir.iterdir()):
        if not degron_dir.is_dir():
            continue
        for mark_dir in sorted(degron_dir.iterdir()):
            if not mark_dir.is_dir():
                continue
            process_mark(degron_dir.name, mark_dir.name, mark_dir, promoters, gene_bodies, out_dir)


if __name__ == "__main__":
    main()
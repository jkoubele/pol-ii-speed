#!/usr/bin/env python3
"""Correlate per-gene ChIP-seq metrics with RNA Pol II pipeline LFCs.

For each (mark, epigenetic metric, LFC type) computes Pearson and Spearman
correlations, saves a summary TSV, and plots scatter plots for top hits.

Usage (from project root):
    python epigenetics/correlate_chipseq_vs_lfc.py
"""

import argparse
import logging
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger(__name__)

LFC_SOURCES = {
    "lfc_elongation_minus_degradation": {
        "file": "gene_specific_pol_2_model/test_results/test_results.tsv",
        "parameter_type": "lfc_elongation_minus_degradation",
    },
    "lfc_splicing_minus_degradation": {
        "file": "gene_specific_pol_2_model/test_results/test_results.tsv",
        "parameter_type": "lfc_splicing_minus_degradation",
    },
    "lfc_elongation_minus_splicing": {
        "file": "gene_specific_splicing_model/test_results/test_results.tsv",
        "parameter_type": "lfc_elongation_minus_splicing",
    },
}

# Metrics where higher = more active mark (for axis label orientation)
POSITIVE_ACTIVITY_METRICS = {
    "promoter_ctrl_mean_fc", "promoter_ctrl_max_fc",
    "promoter_ctrl_n_peaks", "promoter_ctrl_has_peak",
    "promoter_ctrl_max_signal", "promoter_ctrl_sum_signal",
    "body_ctrl_mean_fc", "body_ctrl_max_fc",
    "body_ctrl_n_peaks", "body_ctrl_has_peak",
    "body_ctrl_max_signal", "body_ctrl_sum_signal",
}


def load_lfcs(pipeline_dir):
    """Load all three LFC types; returns DataFrame with gene_name + lfc columns."""
    dfs = []
    for lfc_name, src in LFC_SOURCES.items():
        path = Path(pipeline_dir) / src["file"]
        df = pd.read_csv(path, sep="\t", usecols=["gene_name", "parameter_type", "l2fc_regularized"])
        df = df[df["parameter_type"] == src["parameter_type"]].copy()
        df = df.rename(columns={"gene_name": "gene_id", "l2fc_regularized": lfc_name})
        df = df[["gene_id", lfc_name]].drop_duplicates("gene_id")
        dfs.append(df)

    result = dfs[0]
    for df in dfs[1:]:
        result = result.merge(df, on="gene_id", how="outer")
    log.info(f"Loaded LFCs for {len(result)} genes")
    return result


def load_chipseq_metrics(epigenetics_dir):
    """Return dict {mark: DataFrame} from all chipseq_metrics.tsv files found."""
    metrics = {}
    epi_dir = Path(epigenetics_dir)
    for tsv in sorted(epi_dir.glob("*/chipseq_metrics.tsv")):
        mark = tsv.parent.name
        df = pd.read_csv(tsv, sep="\t")
        metrics[mark] = df
        log.info(f"Loaded {mark}: {len(df)} genes, {len(df.columns)-1} metrics")
    return metrics


def compute_correlations(epi_df, lfc_df, mark):
    """Compute Spearman and Pearson r for all (metric, lfc_type) pairs."""
    merged = epi_df.merge(lfc_df, on="gene_id", how="inner")
    log.info(f"  {mark}: {len(merged)} genes after merge")

    epi_cols = [c for c in epi_df.columns if c != "gene_id"]
    lfc_cols = list(LFC_SOURCES.keys())

    records = []
    for lfc_col in lfc_cols:
        lfc_vals = merged[lfc_col]
        for epi_col in epi_cols:
            epi_vals = merged[epi_col]
            valid = lfc_vals.notna() & epi_vals.notna() & np.isfinite(lfc_vals) & np.isfinite(epi_vals)
            n = valid.sum()
            if n < 20:
                continue
            x, y = epi_vals[valid].values, lfc_vals[valid].values
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sp_r, sp_p = stats.spearmanr(x, y)
                pe_r, pe_p = stats.pearsonr(x, y)
            # Fisher z CI for Spearman r (good approximation at n > ~20)
            z = np.arctanh(sp_r)
            se = 1.0 / np.sqrt(n - 3)
            ci_low = np.tanh(z - 1.96 * se)
            ci_high = np.tanh(z + 1.96 * se)
            records.append({
                "mark": mark,
                "metric": epi_col,
                "lfc_type": lfc_col,
                "spearman_r": sp_r,
                "spearman_ci_low": ci_low,
                "spearman_ci_high": ci_high,
                "spearman_p": sp_p,
                "pearson_r": pe_r,
                "pearson_p": pe_p,
                "n_genes": n,
            })

    return pd.DataFrame(records), merged


def plot_top_scatter(merged, mark, lfc_col, top_metrics, out_dir):
    """Scatter plot grid for the top correlated metrics vs one LFC type."""
    n = len(top_metrics)
    if n == 0:
        return
    ncols = min(n, 3)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = np.array(axes).flatten() if n > 1 else [axes]

    lfc_label = lfc_col.replace("_", " ")
    for ax, (metric, r_val) in zip(axes, top_metrics):
        valid = merged[lfc_col].notna() & merged[metric].notna() & \
                np.isfinite(merged[lfc_col]) & np.isfinite(merged[metric])
        x = merged.loc[valid, metric].values
        y = merged.loc[valid, lfc_col].values
        ax.scatter(x, y, s=4, alpha=0.3, rasterized=True, color="steelblue")
        ax.set_xlabel(metric, fontsize=8)
        ax.set_ylabel(lfc_label, fontsize=8)
        ax.set_title(f"Spearman r = {r_val:.3f}", fontsize=9)
        ax.tick_params(labelsize=7)

    for ax in axes[n:]:
        ax.set_visible(False)

    fig.suptitle(f"{mark}  ×  {lfc_label}", fontsize=11, y=1.01)
    fig.tight_layout()
    fname = out_dir / f"scatter_{mark}_{lfc_col}.pdf"
    fig.savefig(fname, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  -> {fname}")


def plot_correlation_heatmap(corr_df, out_dir, lfc_col, top_n=30):
    """Heatmap of Spearman r: top_n metrics (rows) × marks (columns)."""
    sub = corr_df[corr_df["lfc_type"] == lfc_col].copy()
    if sub.empty:
        return

    pivot = sub.pivot_table(index="metric", columns="mark", values="spearman_r")

    # Keep only top_n rows by max absolute r across marks
    top_idx = pivot.abs().max(axis=1).nlargest(top_n).index
    pivot = pivot.loc[top_idx]

    fig, ax = plt.subplots(figsize=(max(6, len(pivot.columns) * 1.2), max(8, len(pivot) * 0.35)))
    vmax = pivot.abs().max().max()
    im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax, interpolation="none")
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=30, ha="left", fontsize=9)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=8)
    plt.colorbar(im, ax=ax, label="Spearman r", shrink=0.6)
    ax.set_title(f"Spearman r  |  {lfc_col.replace('_', ' ')}  (top {top_n} metrics)", fontsize=10, pad=40)
    fig.tight_layout()
    fname = out_dir / f"heatmap_{lfc_col}.pdf"
    fig.savefig(fname, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  -> {fname}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pipeline-dir",
        default="/home/jakub/Desktop/pol-ii-speed/ENCODE_degron/pipeline_results/BRD4",
        help="Pipeline results dir containing gene_specific_pol_2_model/ and gene_specific_splicing_model/",
    )
    parser.add_argument(
        "--epigenetics-dir",
        default="/home/jakub/Desktop/pol-ii-speed/epigenetics/epigenetic_marks/BRD4",
        help="Directory containing {mark}/chipseq_metrics.tsv files",
    )
    parser.add_argument(
        "--out-dir",
        default="/home/jakub/Desktop/pol-ii-speed/epigenetics/results/BRD4/correlations",
        help="Output directory for correlation TSV and plots",
    )
    parser.add_argument(
        "--top-scatter", type=int, default=6,
        help="Number of top correlations to scatter-plot per LFC type per mark",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    lfc_df = load_lfcs(args.pipeline_dir)
    chipseq = load_chipseq_metrics(args.epigenetics_dir)

    all_corr = []
    merged_by_mark = {}

    for mark, epi_df in chipseq.items():
        log.info(f"Correlating {mark}")
        corr_df, merged = compute_correlations(epi_df, lfc_df, mark)
        all_corr.append(corr_df)
        merged_by_mark[mark] = merged

    corr_all = pd.concat(all_corr, ignore_index=True)
    corr_all["abs_spearman_r"] = corr_all["spearman_r"].abs()
    corr_all = corr_all.sort_values("abs_spearman_r", ascending=False)

    out_tsv = out_dir / "all_correlations.tsv"
    corr_all.drop(columns="abs_spearman_r").to_csv(out_tsv, sep="\t", index=False, float_format="%.6g")
    log.info(f"Saved {len(corr_all)} correlations -> {out_tsv}")

    # Scatter plots: top hits per (mark, lfc_type)
    for mark, merged in merged_by_mark.items():
        for lfc_col in LFC_SOURCES:
            sub = corr_all[(corr_all["mark"] == mark) & (corr_all["lfc_type"] == lfc_col)]
            top = sub.nlargest(args.top_scatter, "abs_spearman_r")[["metric", "spearman_r"]]
            plot_top_scatter(merged, mark, lfc_col, list(zip(top["metric"], top["spearman_r"])), out_dir)

    # Heatmaps: one per LFC type, marks as columns, top metrics as rows
    log.info("Plotting heatmaps")
    for lfc_col in LFC_SOURCES:
        plot_correlation_heatmap(corr_all, out_dir, lfc_col)


if __name__ == "__main__":
    main()

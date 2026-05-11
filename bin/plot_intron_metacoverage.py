#!/usr/bin/env python3

import argparse
import re
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from pathlib import Path
from tqdm import tqdm

NUM_BINS = 100
BIN_LOCATIONS = np.linspace(1 / (2 * NUM_BINS), 1 - 1 / (2 * NUM_BINS), NUM_BINS)


def fit_pi(coverage: np.ndarray) -> float:
    pi_values = np.linspace(0.0, 1.0, 101)[:, np.newaxis]

    trapezoid_density = 1 + pi_values - 2 * BIN_LOCATIONS[np.newaxis, :] * pi_values

    probability = trapezoid_density / trapezoid_density.sum(axis=1, keepdims=True)
    loss = np.dot(-np.log(probability), coverage)

    best_pi = pi_values[loss.argmin()][0]
    return best_pi


def plot_metacoverage(coverages_subset: dict[str, pd.DataFrame], introns_df: pd.DataFrame,
                      output_folder: Path) -> None:
    output_folder.mkdir(exist_ok=True, parents=True)

    introns_of_interest = set(introns_df['intron_name'])
    aggregated_coverage = sum(
        df[df.index.isin(introns_of_interest)].sum(axis=0).values
        for df in coverages_subset.values()
    )

    pi = fit_pi(aggregated_coverage)
    trapezoid_top = (1 + pi - 2 * BIN_LOCATIONS * pi) * aggregated_coverage.sum() / len(BIN_LOCATIONS)
    uniform_top = np.ones_like(BIN_LOCATIONS) * aggregated_coverage.sum() / len(BIN_LOCATIONS) * (1 - pi)

    fig, ax = plt.subplots()
    ax.plot(BIN_LOCATIONS, aggregated_coverage, color='#222222', linewidth=2.5, label='Coverage')
    ax.fill_between(BIN_LOCATIONS, trapezoid_top, y2=uniform_top, alpha=0.45, color='royalblue',
                    label='Transcribing introns')
    ax.fill_between(BIN_LOCATIONS, uniform_top, alpha=0.45, color='darkorange', label='Unspliced introns')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.set_xlabel("Position in intron (5' to 3')")
    ax.set_ylabel("Reads")
    ax.set_ylim(bottom=0)
    ax.set_title(f"Meta-coverage (all introns) \n est. fraction transcribing = {pi:.2f}")
    ax.legend()
    fig.savefig(output_folder / 'all_introns.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    special_ordinals = {1: 'st', 2: 'nd', 3: 'rd'}
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    legend_handles = None
    axes = axes.flatten()

    for ax, quantile in zip(axes, range(1, 5)):
        length_threshold_lower = np.quantile(introns_df['length'], (quantile - 1) / 4)
        length_threshold_upper = np.quantile(introns_df['length'], quantile / 4)

        introns_df_subset = introns_df[
            (length_threshold_lower <= introns_df['length']) & (introns_df['length'] <= length_threshold_upper)]
        quantile_introns = set(introns_df_subset['intron_name'])
        quantile_coverage = sum(
            df[df.index.isin(quantile_introns)].sum(axis=0).values
            for df in coverages_subset.values()
        )

        pi = fit_pi(quantile_coverage)
        trapezoid_top = (1 + pi - 2 * BIN_LOCATIONS * pi) * quantile_coverage.sum() / len(BIN_LOCATIONS)
        uniform_top = np.ones_like(BIN_LOCATIONS) * quantile_coverage.sum() / len(BIN_LOCATIONS) * (1 - pi)

        h1, = ax.plot(BIN_LOCATIONS, quantile_coverage, color='#222222', linewidth=2.5, label='Coverage')
        h2 = ax.fill_between(BIN_LOCATIONS, trapezoid_top, y2=uniform_top, alpha=0.45, color='royalblue',
                              label='Transcribing introns')
        h3 = ax.fill_between(BIN_LOCATIONS, uniform_top, alpha=0.45, color='darkorange', label='Unspliced introns')

        if legend_handles is None:
            legend_handles = [h1, h2, h3]

        ax.set_xlabel("Position in intron (5' to 3')", fontsize=12)
        ax.set_ylabel("Reads", fontsize=12)
        ax.set_ylim(bottom=0)
        ax.set_xlim(left=0, right=1.0)
        ordinal = special_ordinals.get(quantile, 'th')
        ax.set_title(
            f"{quantile}{ordinal} length quantile ({int(length_threshold_lower)} - {int(length_threshold_upper)} bp) \n est. fraction transcribing = {pi:.2f}",
            fontsize=14)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

    fig.suptitle('Meta-coverages (stratified by intron length)', fontsize=18)
    fig.legend(handles=legend_handles, loc='center left', bbox_to_anchor=(0.84, 0.5), frameon=False, fontsize=14)
    fig.subplots_adjust(left=0.08, right=0.82, top=0.87, bottom=0.08, hspace=0.35, wspace=0.25)
    fig.savefig(output_folder / 'introns_stratified_by_length.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def extract_formula_columns(formula: str) -> list[str]:
    formula = formula.lstrip('~').strip()
    tokens = re.split(r'[+:*]', formula)
    seen = set()
    result = []
    for t in tokens:
        t = t.strip()
        if t and t not in seen:
            seen.add(t)
            result.append(t)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--coverage_data_folder',
                        type=Path,
                        required=True,
                        help='Path to the folder with .parquet files with intron coverage.')
    parser.add_argument('--introns_bed_file',
                        type=Path,
                        required=True,
                        help='Path to the .bed files with introns.')
    parser.add_argument('--gene_names',
                        type=Path,
                        required=True,
                        help='Path to the input CSV file with gene names.')
    parser.add_argument('--output_folder',
                        type=Path,
                        default=Path('.'),
                        help='Output folder.')
    parser.add_argument('--samplesheet',
                        type=Path,
                        default=None,
                        help='Path to samplesheet CSV. When provided together with --design_formula, '
                             'produces plots stratified by design columns.')
    parser.add_argument('--design_formula',
                        type=str,
                        default=None,
                        help='R-style design formula (e.g. ~genotype+treatment). Columns present in '
                             'both the formula and the samplesheet are used for stratification.')

    args = parser.parse_args()
    output_folder = args.output_folder
    output_folder.mkdir(exist_ok=True, parents=True)

    genes_df = pd.read_csv(args.gene_names)
    genes_of_interest = set(genes_df['gene_id'])

    introns_df = pd.read_csv(args.introns_bed_file,
                             sep='\t',
                             names=['chromosome', 'start', 'end', 'intron_name', 'score', 'strand'],
                             dtype={"chromosome": "string", "strand": "string"})
    introns_df['gene'] = introns_df['intron_name'].apply(lambda x: x.split('_')[0])
    introns_df['length'] = introns_df['end'] - introns_df['start']
    introns_df = introns_df[introns_df['gene'].isin(genes_of_interest)]

    coverages_by_sample: dict[str, pd.DataFrame] = {}
    for file in tqdm(args.coverage_data_folder.iterdir(), desc='Loading coverage data'):
        if file.suffix != '.parquet':
            continue
        coverage_df = pd.read_parquet(file).set_index('intron_name')
        coverages_by_sample[file.stem] = coverage_df

    if args.samplesheet is not None and args.design_formula is not None:
        samplesheet_df = pd.read_csv(args.samplesheet)
        required_cols = {'sample', 'fastq_1', 'fastq_2', 'strandedness'}

        for col in extract_formula_columns(args.design_formula):
            if col not in samplesheet_df.columns:
                print(f"Warning: formula term '{col}' not found in samplesheet columns, skipping.")
                continue
            if col in required_cols:
                continue
            for level in samplesheet_df[col].dropna().unique():
                samples_in_level = set(samplesheet_df.loc[samplesheet_df[col] == level, 'sample'].astype(str))
                subset = {k: v for k, v in coverages_by_sample.items() if k in samples_in_level}
                if not subset:
                    continue
                level_str = re.sub(r'[^\w\-]', '_', str(level))
                plot_metacoverage(subset, introns_df, output_folder / col / level_str)
    else:
        plot_metacoverage(coverages_by_sample, introns_df, output_folder)
import pandas as pd
import argparse
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
from tqdm import tqdm
from typing import NamedTuple, Optional, Iterator


class IntronSetSpecification(NamedTuple):
    title: str
    file_basename: str
    length_quantile_lower: float
    length_quantile_upper: float


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

    # For development
    import sys

    results_folder = Path('/cellfile/projects/pol_ii_speed/jkoubele/analysis/senescent_cells/results/')
    # model_run_name = 'model_run_07_Apr_2026_15_33_32'
    sys.argv = [
        'script_name.py',
        '--coverage_data_folder',
        str(results_folder / 'preprocessing' / 'rescaled_coverage'),
        '--introns_bed_file',
        str(results_folder / 'preprocessing' / 'genomic_features' / 'introns.bed'),
        '--gene_names',
        str(results_folder / 'preprocessing' / 'genomic_features' / 'protein_coding_genes.csv')
    ]

    args = parser.parse_args()

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
        coverage_df = pd.read_parquet(file).set_index('intron_name')
        coverages_by_sample[file.stem] = coverage_df

    intron_specifications = [
        IntronSetSpecification(title='All introns',
                               file_basename='all_introns',
                               length_quantile_lower=0.0,
                               length_quantile_upper=1.0),
        IntronSetSpecification(title='Introns 1st length quartile',
                               file_basename='introns_quartile_1',
                               length_quantile_lower=0.0,
                               length_quantile_upper=0.25),
        IntronSetSpecification(title='Introns 2nd length quartile',
                               file_basename='introns_quartile_2',
                               length_quantile_lower=0.25,
                               length_quantile_upper=0.5),
        IntronSetSpecification(title='Introns 3rd length quartile',
                               file_basename='introns_quartile_3',
                               length_quantile_lower=0.5,
                               length_quantile_upper=0.75),
        IntronSetSpecification(title='Introns 4th length quartile',
                               file_basename='introns_quartile_4',
                               length_quantile_lower=0.75,
                               length_quantile_upper=1.0),
    ]
    for intron_specification in intron_specifications:
        length_threshold_lower = np.quantile(introns_df['length'], intron_specification.length_quantile_lower)
        length_threshold_upper = np.quantile(introns_df['length'], intron_specification.length_quantile_upper)

        introns_df_subset = introns_df[
            (length_threshold_lower <= introns_df['length']) & (introns_df['length'] <= length_threshold_upper)]
        introns_of_interest = set(introns_df_subset['intron_name'])

        aggregated_coverage = sum(
            [df[df.index.isin(introns_of_interest)].sum(axis=0).values for df in coverages_by_sample.values()])

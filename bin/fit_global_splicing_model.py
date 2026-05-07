#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from tqdm import tqdm



from pol_ii_speed_modeling.load_dataset import load_dataset_metadata
from pol_ii_speed_modeling.train import get_splicing_model_results

# --- Dev block: comment out when running via Nextflow ---
results_folder = Path('/home/jakub/Desktop/c_elegans_test/results/')
model_run_timestamp = 'model_run_26_Apr_2026_11_47_29'
# sys.argv = [
#     'fit_global_splicing_model.py',
#     '--modeled_introns', str(results_folder / f'modeling/{model_run_timestamp}/modeled_genes/modeled_introns.tsv'),
#     '--design_matrix', str(results_folder / f'modeling/{model_run_timestamp}/design_matrices/design_matrix.tsv'),
#     '--library_size_factors', str(results_folder / 'preprocessing/aggregated_counts/library_size_factors.tsv'),
#     '--lrt_metadata', str(results_folder / f'modeling/{model_run_timestamp}/design_matrices/lrt_metadata.tsv'),
#     '--reduced_matrices_folder', str(results_folder / f'modeling/{model_run_timestamp}/design_matrices/reduced_design_matrices'),
#     '--coverage_data_folder', str(results_folder / 'preprocessing/rescaled_coverage/'),
#     '--output_folder', '/home/jakub/Desktop/c_elegans_test/results/test_splicing_model/',
# ]
# --- End dev block ---

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--modeled_introns',
                        type=Path,
                        required=True,
                        help='Path to the input TSV file with modeled introns.')
    parser.add_argument('--design_matrix',
                        type=Path,
                        required=True,
                        help='Path to the input TSV file with design matrix.')
    parser.add_argument('--library_size_factors',
                        type=Path,
                        required=True,
                        help='Path to the input TSV file with library size factors.')
    parser.add_argument('--lrt_metadata',
                        type=Path,
                        required=True,
                        help='Path to the TSV file with LRT metadata.')
    parser.add_argument('--reduced_matrices_folder',
                        type=Path,
                        required=True,
                        help='Folder with TSV files with reduced matrices for LRT.')
    parser.add_argument('--coverage_data_folder',
                        type=Path,
                        required=True,
                        help='Path to the folder with .parquet files with intron coverage.')
    parser.add_argument('--output_folder',
                        type=Path,
                        default=Path('.'),
                        help='Output folder.')
    args = parser.parse_args()

    args.output_folder.mkdir(exist_ok=True, parents=True)

    dataset_metadata = load_dataset_metadata(
        design_matrix_file=args.design_matrix,
        library_size_factors_file=args.library_size_factors,
        lrt_metadata_file=args.lrt_metadata,
        reduced_matrices_folder=args.reduced_matrices_folder,
    )

    modeled_introns_df = pd.read_csv(args.modeled_introns, sep='\t')
    intron_names = modeled_introns_df['intron_id'].tolist()

    coverage_df_by_sample: dict[str, pd.DataFrame] = {}
    for sample in tqdm(dataset_metadata.sample_names, desc='Loading coverage data'):
        sample_coverage_df = pd.read_parquet(
            args.coverage_data_folder / f'{sample}.parquet'
        ).set_index('intron_name')
        coverage_df_by_sample[sample] = sample_coverage_df.loc[intron_names]

    coverage = torch.tensor(
        np.stack([coverage_df_by_sample[sample].values for sample in dataset_metadata.sample_names]),
        dtype=torch.float32,
    )

    model_param_df, test_results_df, _ = get_splicing_model_results(
        coverage=coverage,
        dataset_metadata=dataset_metadata,
        intron_names=intron_names,
        verbose=True,
        compute_wald_test=False,
    )

    model_param_df.reset_index(drop=True).to_csv(
        args.output_folder / 'model_parameters.tsv', sep='\t', index=False
    )
    test_results_df.to_csv(args.output_folder / 'test_results.tsv', sep='\t', index=False)

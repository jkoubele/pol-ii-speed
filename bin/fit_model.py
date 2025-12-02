#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from pol_ii_speed_modeling.load_dataset import load_dataset_metadata, load_gene_data_list
from pol_ii_speed_modeling.train import get_results_for_gene


def parse_bool_in_argparse(argument: str) -> bool:
    if argument.lower() == 'true':
        return True
    elif argument.lower() == 'false':
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected. Use either 'true' or 'false'.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--design_matrix',
                        type=Path,
                        required=True,
                        help='Path to the input CSV file with design matrix.')
    parser.add_argument('--gene_names',
                        type=Path,
                        required=True,
                        help='Path to the input CSV file with gene names.')
    parser.add_argument('--exon_counts',
                        type=Path,
                        required=True,
                        help='Path to the input .tsv file with exon counts.')
    parser.add_argument('--intron_counts',
                        type=Path,
                        required=True,
                        help='Path to the input .tsv file with intron counts.')
    parser.add_argument('--library_size_factors',
                        type=Path,
                        required=True,
                        help='Path to the input .tsv file with library size factors.')
    parser.add_argument('--isoform_length_factors',
                        type=Path,
                        required=True,
                        help='Path to the input .tsv file with isoform length factors.')
    parser.add_argument('--coverage_data_folder',
                        type=Path,
                        required=True,
                        help='Path to the folder with .parquet files with intron coverage.')
    parser.add_argument('--lrt_metadata',
                        type=Path,
                        required=True,
                        help='Path to the .csv file with LRT metadata.')
    parser.add_argument('--reduced_matrices_folder',
                        type=Path,
                        required=True,
                        help='Folder with CSV files with reduced matrices for LRT.')
    parser.add_argument('--intron_specific_lfc',
                        type=parse_bool_in_argparse,
                        required=True,
                        help='Whether to use intron-specific LFC (true or false).')
    parser.add_argument('--output_folder',
                        type=Path,
                        default=Path('.'),
                        help='Output folder.')
    parser.add_argument('--output_name_suffix',
                        type=str,
                        default='',
                        help='Suffix for the output CSV files.')
    # For development
    import sys
    #
    sys.argv = [
        'script_name.py',
        '--design_matrix',
        '/cellfile/projects/pol_ii_speed/jkoubele/analysis/EU_seq_Joris/results/modeling/model_run_01_Dec_2025_01_01_40/design_matrices/design_matrix.csv',
        '--gene_names',
        '/cellfile/projects/pol_ii_speed/jkoubele/analysis/EU_seq_Joris/results/preprocessing/gene_names/test_10.csv',
        '--exon_counts',
        '/cellfile/projects/pol_ii_speed/jkoubele/analysis/EU_seq_Joris/results/preprocessing/aggregated_counts/exon_counts.tsv',
        '--intron_counts',
        '/cellfile/projects/pol_ii_speed/jkoubele/analysis/EU_seq_Joris/results/preprocessing/aggregated_counts/intron_counts.tsv',
        '--library_size_factors',
        '/cellfile/projects/pol_ii_speed/jkoubele/analysis/EU_seq_Joris/results/preprocessing/aggregated_counts/library_size_factors.tsv',
        '--isoform_length_factors',
        '/cellfile/projects/pol_ii_speed/jkoubele/analysis/EU_seq_Joris/results/preprocessing/aggregated_counts/isoform_length_factors.tsv',
        '--coverage_data_folder',
        '/cellfile/projects/pol_ii_speed/jkoubele/analysis/EU_seq_Joris/results/preprocessing/rescaled_coverage',
        '--reduced_matrices_folder',
        '/cellfile/projects/pol_ii_speed/jkoubele/analysis/EU_seq_Joris/results/modeling/model_run_01_Dec_2025_01_01_40/design_matrices/reduced_design_matrices',
        '--lrt_metadata',
        '/cellfile/projects/pol_ii_speed/jkoubele/analysis/EU_seq_Joris/results/modeling/model_run_01_Dec_2025_01_01_40/design_matrices/lrt_metadata.csv',
        '--intron_specific_lfc', 'false',
        '--output_folder', '/cellfile/projects/pol_ii_speed/jkoubele/pol-ii-speed/test_regularization',
        '--output_name_suffix', '_dev'
    ]

    args = parser.parse_args()

    output_folder = args.output_folder
    output_folder.mkdir(exist_ok=True, parents=True)
    dataset_metadata = load_dataset_metadata(design_matrix_file=args.design_matrix,
                                             library_size_factors_file=args.library_size_factors,
                                             lrt_metadata_file=args.lrt_metadata,
                                             reduced_matrices_folder=args.reduced_matrices_folder)

    gene_data_list = load_gene_data_list(gene_names_file=args.gene_names,
                                         exon_counts_file=args.exon_counts,
                                         intron_counts_file=args.intron_counts,
                                         isoform_length_factors_file=args.isoform_length_factors,
                                         coverage_folder=args.coverage_data_folder,
                                         sample_names=dataset_metadata.sample_names,
                                         log_output_folder=Path("./logs"),
                                         log_output_name_suffix=args.output_name_suffix)

    if gene_data_list:  # gene_data_list may be empty if all genes are filtered out in the load_gene_data_list()
        result_list = [get_results_for_gene(gene_data=gene_data,
                                            dataset_metadata=dataset_metadata,
                                            intron_specific_lfc=args.intron_specific_lfc)
                       for gene_data in tqdm(gene_data_list)]
        all_model_param_df = pd.concat([result[0] for result in result_list]).reset_index(drop=True)
        all_test_results_df = pd.concat([result[1] for result in result_list]).reset_index(drop=True)

        all_model_param_df.to_csv(output_folder / f"model_parameters{args.output_name_suffix}.csv", index=False)
        all_test_results_df.to_csv(output_folder / f"test_results{args.output_name_suffix}.csv", index=False)

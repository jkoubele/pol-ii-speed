#!/usr/bin/env python3

import argparse
from pathlib import Path
from modeling.load_dataset import load_dataset_matedata

import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--design_matrix',
                        type=Path,
                        required=True,
                        help='Path to the input .csv file with design matrix.')
    parser.add_argument('--gene_names',
                        type=Path,
                        required=True,
                        help='Path to the input .csv file with gene names.')
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
    parser.add_argument('--coverage_data_folder',
                        type=Path,
                        required=True,
                        help='Path to the folder with .parquet files with intron coverage.')
    parser.add_argument('--output_folder',
                        type=Path,
                        default=Path('.'),
                        help='Output folder.')
    parser.add_argument('--output_basename',
                        type=str,
                        default='model_results',
                        help='Basename of the output csv file.')
    # For development
    # import sys
    # sys.argv = [
    #         'script_name.py',
    #         '--design_matrix', '/cellfile/datapublic/jkoubele/drosophila_mutants/results/design_matrix/design_matrix.csv',
    #         '--gene_names', '/cellfile/datapublic/jkoubele/drosophila_mutants/results/gene_names_chunks/gene_names_chunk_01.csv',
    #         '--exon_counts', '/cellfile/datapublic/jkoubele/drosophila_mutants/results/aggregated_counts/exon_counts.tsv',
    #         '--intron_counts', '/cellfile/datapublic/jkoubele/drosophila_mutants/results/aggregated_counts/intron_counts.tsv',
    #         '--library_size_factors', '/cellfile/datapublic/jkoubele/drosophila_mutants/results/aggregated_counts/library_size_factors.tsv',
    #         '--coverage_data_folder', '/cellfile/datapublic/jkoubele/drosophila_mutants/results/rescaled_coverage',
    #         '--output_folder', '/cellfile/datapublic/jkoubele/drosophila_mutants/results/model_results',
    #         '--output_basename', 'chunk_001'
    #     ]

    args = parser.parse_args()    
    
    output_folder = args.output_folder
    output_folder.mkdir(exist_ok=True, parents=True)    
    dataset_matedata = load_dataset_matedata(design_matrix_file=args.design_matrix,
                                             library_size_factors_file=args.library_size_factors)

    design_matrix_df = pd.read_csv(args.design_matrix)
    gene_names_df = pd.read_csv(args.gene_names)
    exon_counts_df = pd.read_csv(args.exon_counts, sep='\t')
    intron_counts_df = pd.read_csv(args.intron_counts, sep='\t')
    library_size_factors_df = pd.read_csv(args.library_size_factors, sep='\t')
    
    df_out = gene_names_df.copy()
    df_out['dummy_results'] = df_out['gene_id'].apply(lambda x: f"Result for gene {x}")
    
    df_out.to_csv(output_folder / f"{args.output_basename}.csv", index=False)


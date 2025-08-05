#!/usr/bin/env python3

import argparse
from pathlib import Path

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

    args = parser.parse_args()

    design_matrix_df = pd.read_csv(args.design_matrix)

    design_matrix_df = pd.read_csv(args.design_matrix)
    gene_names_df = pd.read_csv(args.gene_names)
    exon_counts_df = pd.read_csv(args.exon_counts, sep='\t')
    intron_counts_df = pd.read_csv(args.intron_counts, sep='\t')
    library_size_factors_df = pd.read_csv(args.library_size_factors, sep='\t')

    with open("success.txt", "w") as f:
        print("All data loaded!", file=f)

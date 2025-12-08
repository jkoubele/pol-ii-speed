#!/usr/bin/env python3

import argparse
import pandas as pd
import torch
from pathlib import Path
from tqdm import tqdm

from pol_ii_speed_modeling.train import get_regularized_model_results, CacheForRegularization, StateDict
from pol_ii_speed_modeling.pol_ii_model import GeneData, DatasetMetadata

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--cache_for_regularization',
                        type=Path,
                        required=True,
                        help='Path to the pickle file with CacheForRegularization object.')
    parser.add_argument('--regularization_coefficients',
                        type=Path,
                        required=True,
                        help='Path to the CSV file with regularization coefficients.')
    parser.add_argument('--output_folder',
                        type=Path,
                        default=Path('.'),
                        help='Output folder.')
    parser.add_argument('--output_name_suffix',
                        type=str,
                        default='',
                        help='Suffix for the output CSV files.')

    args = parser.parse_args()
    output_folder = args.output_folder
    output_folder.mkdir(exist_ok=True, parents=True)

    regularization_coefficients_df = pd.read_csv(args.regularization_coefficients)

    device = 'cpu'
    cache_for_regularization = torch.load(cache_for_regularization_file,
                                          weights_only=False,
                                          map_location=device)
    dataset_metadata = cache_for_regularization.dataset_metadata
    intron_specific_lfc = cache_for_regularization.intron_specific_lfc
    regularized_model_params_list: list[pd.DataFrame] = [get_regularized_model_results(gene_data=gene_data,
                                                                                       dataset_metadata=dataset_metadata,
                                                                                       intron_specific_lfc=intron_specific_lfc,
                                                                                       hot_start_state_dict=model_state_dict,
                                                                                       regularization_coefficients_df=None,
                                                                                       device=device) for
                                                         gene_data, model_state_dict in
                                                         tqdm(cache_for_regularization.training_input_per_gene)]
    all_regularized_model_param_df = pd.concat(regularized_model_params_list).reset_index(drop=True)

    all_regularized_model_param_df.to_csv(output_folder / f"regularized_model_parameters{args.output_name_suffix}.csv",
                                          index=False)

import argparse
from pathlib import Path
import torch
import numpy as np

import pandas as pd
from tqdm import tqdm

from pol_ii_speed_modeling.load_dataset import load_dataset_metadata
from pol_ii_speed_modeling.pol_ii_model import SplicingModel, DatasetMetadata


results_folder = Path('/cellfile/projects/pol_ii_speed/jkoubele/analysis/c_elegans_test/results')
model_run_timestamp = 'model_run_20_Apr_2026_19_05_29'

dataset_metadata = load_dataset_metadata(
    design_matrix_file= results_folder / f"modeling/{model_run_timestamp}/design_matrices/design_matrix.tsv",
    library_size_factors_file=results_folder / "preprocessing/aggregated_counts/library_size_factors.tsv",
    lrt_metadata_file= results_folder / f'modeling/{model_run_timestamp}/design_matrices/lrt_metadata.tsv',
    reduced_matrices_folder=results_folder / f'modeling/{model_run_timestamp}/design_matrices/reduced_design_matrices')

modeled_introns_df = pd.read_csv(results_folder / f'modeling/{model_run_timestamp}/modeled_genes/modeled_introns.tsv', sep='\t')
modeled_introns_df = modeled_introns_df.head(100)

rescaled_coverage_folder = results_folder / 'preprocessing/rescaled_coverage/'

coverage_df_by_sample: dict[str, pd.DataFrame] = {}
for sample in tqdm(dataset_metadata.sample_names, desc = 'Loading coverage data'):
    sample_coverage_df = pd.read_parquet(rescaled_coverage_folder / f'{sample}.parquet').set_index('intron_name')
    coverage_df_by_sample[sample] = sample_coverage_df.loc[modeled_introns_df['intron_id']]
    
coverage = torch.tensor(np.stack([coverage_df_by_sample[sample].values for sample in dataset_metadata.sample_names]), dtype=torch.float32)



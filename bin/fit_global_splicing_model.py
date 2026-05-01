from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy import stats
from tqdm import tqdm

from pol_ii_speed_modeling.load_dataset import load_dataset_metadata
from pol_ii_speed_modeling.pol_ii_model import DatasetMetadata, SplicingModel
from pol_ii_speed_modeling.train import make_lbfgs_optimizer, train_model, make_splicing_closures

# results_folder = Path('/cellfile/projects/pol_ii_speed/jkoubele/analysis/c_elegans_test/results')
results_folder = Path('/home/jakub/Desktop/c_elegans_test/results/')
model_run_timestamp = 'model_run_26_Apr_2026_11_47_29'
# model_run_timestamp = 'model_run_20_Apr_2026_19_05_29'

dataset_metadata = load_dataset_metadata(
    design_matrix_file=results_folder / f"modeling/{model_run_timestamp}/design_matrices/design_matrix.tsv",
    library_size_factors_file=results_folder / "preprocessing/aggregated_counts/library_size_factors.tsv",
    lrt_metadata_file=results_folder / f'modeling/{model_run_timestamp}/design_matrices/lrt_metadata.tsv',
    reduced_matrices_folder=results_folder / f'modeling/{model_run_timestamp}/design_matrices/reduced_design_matrices')

modeled_introns_df = pd.read_csv(results_folder / f'modeling/{model_run_timestamp}/modeled_genes/modeled_introns.tsv',
                                 sep='\t')
# modeled_introns_df = modeled_introns_df.head(1000)

rescaled_coverage_folder = results_folder / 'preprocessing/rescaled_coverage/'

coverage_df_by_sample: dict[str, pd.DataFrame] = {}
for sample in tqdm(dataset_metadata.sample_names, desc='Loading coverage data'):
    sample_coverage_df = pd.read_parquet(rescaled_coverage_folder / f'{sample}.parquet').set_index('intron_name')
    coverage_df_by_sample[sample] = sample_coverage_df.loc[modeled_introns_df['intron_id']]

coverage = torch.tensor(
    np.stack([coverage_df_by_sample[sample].values for sample in dataset_metadata.sample_names]),
    dtype=torch.float32
)
intron_names = modeled_introns_df['intron_id'].tolist()


def get_splicing_model_results(
        coverage: torch.Tensor,
        dataset_metadata: DatasetMetadata,
        intron_names: list[str],
        device: str = 'cpu',
) -> tuple[pd.DataFrame, pd.DataFrame]:
    coverage = coverage.to(device)
    design_matrix = dataset_metadata.design_matrix.to(device)
    model_full = SplicingModel(
        feature_names=dataset_metadata.feature_names,
        intron_names=intron_names,
    ).to(device)
    model_full.initialize_theta(coverage)

    optimizer_full = make_lbfgs_optimizer(model_full)
    closure_full, evaluate_loss_full = make_splicing_closures(
        model_full, optimizer_full, coverage, design_matrix,
    )
    model_full, results_full = train_model(model_full, optimizer_full, closure_full, evaluate_loss_full)

    model_param_df = model_full.get_param_df()
    model_param_df['loss_full_model'] = results_full.final_loss
    model_param_df['training_diverged_full_model'] = results_full.training_diverged
    model_param_df['training_converged_within_max_epochs_full_model'] = results_full.converged_within_max_epochs
    model_param_df = model_param_df.set_index(['parameter_type', 'feature_name', 'intron_name'], drop=False)

    test_results_list: list[dict] = []
    for _, lrt_row in dataset_metadata.lrt_metadata.iterrows():
        reduced_matrix = dataset_metadata.reduced_matrices[lrt_row['test_id']].to(device)
        num_reduced_features = reduced_matrix.shape[1]
        placeholder_names = [f'reduced_feature_{i}' for i in range(num_reduced_features)]

        model_reduced = SplicingModel(
            feature_names=placeholder_names,
            intron_names=intron_names,
        ).to(device)

        with torch.no_grad():
            model_reduced.theta.data.copy_(model_full.theta.data)
            if num_reduced_features > 0:
                full_lfc_contribution = design_matrix @ model_full.lfc
                model_reduced.lfc.data.copy_(
                    torch.linalg.lstsq(reduced_matrix, full_lfc_contribution).solution
                )

        optimizer_reduced = make_lbfgs_optimizer(model_reduced)
        closure_reduced, evaluate_loss_reduced = make_splicing_closures(
            model_reduced, optimizer_reduced, coverage, reduced_matrix,
        )
        model_reduced, results_reduced = train_model(
            model_reduced, optimizer_reduced, closure_reduced, evaluate_loss_reduced,
        )

        lfc_positive = 0.0 if pd.isna(lrt_row['lfc_column_positive']) else \
            model_param_df.loc[('lfc', lrt_row['lfc_column_positive'], None)]['value']
        lfc_negative = 0.0 if pd.isna(lrt_row['lfc_column_negative']) else \
            model_param_df.loc[('lfc', lrt_row['lfc_column_negative'], None)]['value']

        chi2_stat = 2 * (results_reduced.final_loss - results_full.final_loss)
        test_result = lrt_row.to_dict()
        test_result['lfc'] = lfc_positive - lfc_negative
        test_result['loss_full_model'] = results_full.final_loss
        test_result['loss_reduced_model'] = results_reduced.final_loss
        test_result['chi2_test_statistics'] = chi2_stat
        test_result['p_value'] = 1 - stats.chi2.cdf(chi2_stat, df=lrt_row['lrt_df'])
        test_result['training_diverged_reduced_model'] = results_reduced.training_diverged
        test_result['training_converged_within_max_epochs_reduced_model'] = results_reduced.converged_within_max_epochs
        test_results_list.append(test_result)

    return model_param_df, pd.DataFrame(test_results_list)


model_param_df, test_results_df = get_splicing_model_results(
    coverage=coverage,
    dataset_metadata=dataset_metadata,
    intron_names=intron_names,
)

print(model_param_df)
print(test_results_df)

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
import scipy
import torch
from scipy import stats
from torch import optim
from torch.func import functional_call, hessian

from pol_ii_speed_modeling.pol_ii_model import GeneData, DatasetMetadata, Pol2TotalLoss, Pol2Model


@dataclass
class TrainingResults:
    final_loss: float
    converged_within_max_epochs: bool
    num_epochs: int
    losses: list[float]
    training_diverged: bool


def train_model(gene_data: GeneData,
                dataset_metadata: DatasetMetadata,
                pol_2_total_loss: Pol2TotalLoss,
                device='cpu',
                max_epochs=200,
                max_patience=5,
                loss_change_tolerance=1e-6,
                model: Optional[Pol2Model] = None,
                intron_specific_lfc=True
                ) -> tuple[Pol2Model, TrainingResults]:
    if model is None:
        model = Pol2Model(feature_names=dataset_metadata.feature_names,
                          intron_names=gene_data.intron_names,
                          intron_specific_lfc=intron_specific_lfc).to(device)
        model.initialize_intercepts(gene_data, dataset_metadata.library_sizes)

    optimizer = optim.LBFGS(model.parameters(),
                            lr=0.05,
                            max_iter=20,
                            tolerance_change=1e-09,
                            tolerance_grad=1e-07,
                            history_size=100,
                            line_search_fn='strong_wolfe')

    def closure():
        optimizer.zero_grad()
        predicted_reads_exon, predicted_reads_intron, phi = model(dataset_metadata.design_matrix,
                                                                  dataset_metadata.log_library_sizes)
        loss = pol_2_total_loss(reads_exon=gene_data.exon_reads,
                                reads_introns=gene_data.intron_reads,
                                coverage=gene_data.coverage,
                                predicted_reads_exon=predicted_reads_exon,
                                predicted_reads_intron=predicted_reads_intron,
                                phi=phi)
        loss.backward()
        return loss

    previous_loss = None
    patience_counter = 0

    best_loss = np.inf
    best_state_dict = {k: v.detach().clone() for k, v in model.state_dict().items()}

    converged_within_max_epochs = False
    training_diverged = False

    losses: list[float] = []
    for epoch in range(max_epochs):
        loss = optimizer.step(closure).item()
        losses.append(loss)
        if np.isnan(loss):
            training_diverged = True
            break

        if loss < best_loss:
            best_loss = loss
            best_state_dict = {k: v.detach().clone() for k, v in model.state_dict().items()}

        if previous_loss is not None:
            relative_loss_change = abs(loss - previous_loss) / (abs(previous_loss) + 1e-10)
            if relative_loss_change < loss_change_tolerance:
                patience_counter += 1
                if patience_counter >= max_patience:
                    converged_within_max_epochs = True
                    break
            else:
                patience_counter = 0

        previous_loss = loss

    if best_loss < loss or np.isnan(loss):
        model.load_state_dict(best_state_dict)

    training_results = TrainingResults(final_loss=loss,
                                       converged_within_max_epochs=converged_within_max_epochs,
                                       num_epochs=epoch + 1,
                                       losses=losses,
                                       training_diverged=training_diverged
                                       )

    return model, training_results


def flatten_hessian_dict(hessian_dict, model_parameters):
    param_names = list(model_parameters.keys())
    param_sizes = {name: model_parameters[name].numel() for name in param_names}
    param_offsets = {name: sum(param_sizes[n] for n in param_names[:i])
                     for i, name in enumerate(param_names)}
    total_param_size = sum(param_sizes.values())
    hessian_matrix = torch.zeros((total_param_size, total_param_size))

    for param_1 in param_names:
        for param_2 in param_names:
            block = hessian_dict[param_1][param_2].detach()
            block = block.reshape(param_sizes[param_1], param_sizes[param_2])
            hessian_matrix[
            param_offsets[param_1]:param_offsets[param_1] + param_sizes[param_1],
            param_offsets[param_2]:param_offsets[param_2] + param_sizes[param_2]
            ] = block

    return hessian_matrix


def get_fisher_information_matrix(model: Pol2Model,
                                  pol_2_total_loss: Pol2TotalLoss,
                                  gene_data: GeneData,
                                  dataset_metadata: DatasetMetadata) -> torch.Tensor:
    model_parameters = dict(model.named_parameters())

    def loss_by_model_parameters(model_parameters):
        outputs = functional_call(model, model_parameters, (dataset_metadata.design_matrix,
                                                            dataset_metadata.log_library_sizes))
        predicted_reads_exon, predicted_reads_intron, phi = outputs
        return pol_2_total_loss(reads_exon=gene_data.exon_reads,
                                reads_introns=gene_data.intron_reads,
                                coverage=gene_data.coverage,
                                predicted_reads_exon=predicted_reads_exon,
                                predicted_reads_intron=predicted_reads_intron,
                                phi=phi)

    hessian_dict = hessian(loss_by_model_parameters)(model_parameters)
    hessian_matrix = flatten_hessian_dict(hessian_dict, model_parameters)
    return hessian_matrix


def add_wald_test_results(df_param: pd.DataFrame, hessian_matrix: torch.Tensor) -> pd.DataFrame:
    """
    Adds results of the Wald test to the dataframe with model parameters.
    """
    df_param = df_param.copy()
    rank = torch.linalg.matrix_rank(hessian_matrix)

    # Rank-revealing QR is currently not available in Pytorch (see https://github.com/pytorch/pytorch/issues/10454),
    # so we are going to use SciPy implementation.
    _, _, qr_pivots = scipy.linalg.qr(hessian_matrix.cpu().numpy(), pivoting=True)
    identifiable_indices = np.sort(qr_pivots[:rank])
    hessian_subset = hessian_matrix[identifiable_indices][:, identifiable_indices]

    assert set(qr_pivots) == set(df_param.index)

    # hessian_subset should have a full rank by the way it is constructed:
    # assert hessian_subset.shape[0] == torch.linalg.matrix_rank(hessian_subset)    

    cov_matrix = torch.linalg.pinv(hessian_subset).cpu()
    standard_errors = torch.sqrt(torch.diag(cov_matrix)).numpy()

    df_param['SE'] = np.nan
    df_param['z_score'] = np.nan
    df_param['p_value_wald'] = np.nan
    df_param['identifiable'] = df_param.index.isin(identifiable_indices)

    df_param.loc[identifiable_indices, 'SE'] = standard_errors
    df_param.loc[identifiable_indices, 'z_score'] = df_param.loc[identifiable_indices, 'value'] / standard_errors
    df_param.loc[identifiable_indices, 'p_value_wald'] = 2 * (
            1 - stats.norm.cdf(np.abs(df_param.loc[identifiable_indices, 'z_score'])))

    return df_param


def get_results_for_gene(gene_data: GeneData,
                         dataset_metadata: DatasetMetadata,
                         device='cpu',
                         perform_wald_test=True,
                         perform_lrt=False,
                         intron_specific_lfc=True) -> pd.DataFrame:
    gene_data = gene_data.to(device)
    dataset_metadata = dataset_metadata.to(device)
    pol_2_total_loss = Pol2TotalLoss().to(device)

    model, training_results = train_model(gene_data=gene_data,
                                          dataset_metadata=dataset_metadata,
                                          pol_2_total_loss=pol_2_total_loss,
                                          device=device,
                                          intron_specific_lfc=intron_specific_lfc)

    param_df = model.get_param_df()
    param_df['gene_name'] = gene_data.gene_name
    param_df['loss_unrestricted'] = training_results.final_loss

    if perform_wald_test:
        fisher_information_matrix = get_fisher_information_matrix(model=model,
                                                                  pol_2_total_loss=pol_2_total_loss,
                                                                  gene_data=gene_data,
                                                                  dataset_metadata=dataset_metadata)

        param_df = add_wald_test_results(param_df, fisher_information_matrix)

    if perform_lrt:
        loss_unrestricted = training_results.final_loss

        loss_differences: list[Optional[float]] = []

        for _, row in param_df.iterrows():
            if row['parameter_type'] not in ('alpha', 'beta', 'gamma'):
                loss_differences.append(None)
            else:
                model_restricted = Pol2Model(feature_names=dataset_metadata.feature_names,
                                             intron_names=gene_data.intron_names,
                                             intron_specific_lfc=intron_specific_lfc).to(device)
                model_restricted.load_state_dict(model.state_dict())

                model_restricted.set_parameter_mask(
                    param_name=row['parameter_type'],
                    feature_name=row['feature_name'],
                    intron_name=None if (not intron_specific_lfc or row['parameter_type'] == 'alpha') else row[
                        'intron_name'],
                    value=0.0)
                model_restricted, training_results = train_model(gene_data=gene_data,
                                                                 dataset_metadata=dataset_metadata,
                                                                 pol_2_total_loss=pol_2_total_loss,
                                                                 device=device,
                                                                 model=model_restricted,
                                                                 intron_specific_lfc=intron_specific_lfc)
                if training_results.converged_within_max_epochs:
                    loss_differences.append(training_results.final_loss - loss_unrestricted)
                else:
                    loss_differences.append(None)

        param_df['loss_differences'] = loss_differences
        param_df['loss_restricted'] = param_df['loss_unrestricted'] + param_df['loss_differences']
        param_df['p_value_lrt'] = 1 - stats.chi2.cdf(2 * param_df['loss_differences'], df=1)
    return param_df

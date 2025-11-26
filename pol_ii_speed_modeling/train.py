import numpy as np
import pandas as pd
import scipy
import torch
from dataclasses import dataclass
from scipy import stats
from torch import optim
from torch.func import functional_call, hessian
from typing import Optional

from pol_ii_speed_modeling.pol_ii_model import GeneData, DatasetMetadata, Pol2TotalLoss, Pol2Model, TestableParameters, \
    LRTSpecification

LOSS_CLAMP_VALUE = 1e30


@dataclass
class TrainingResults:
    final_loss: float
    converged_within_max_epochs: bool
    num_epochs: int
    losses: list[float]
    training_diverged: bool


def train_model(model: Pol2Model,
                gene_data: GeneData,
                dataset_metadata: DatasetMetadata,
                pol_2_total_loss: Pol2TotalLoss,
                reduced_matrix_name: Optional[str] = None,
                device='cpu',
                max_epochs=200,
                max_patience=5,
                loss_change_tolerance=1e-6) -> tuple[Pol2Model, TrainingResults]:
    reduced_matrix = None if reduced_matrix_name is None else dataset_metadata.reduced_matrices[reduced_matrix_name]

    optimizer = optim.LBFGS(model.parameters(),
                            lr=0.05,
                            max_iter=20,
                            tolerance_change=1e-09,
                            tolerance_grad=1e-07,
                            history_size=100,
                            line_search_fn='strong_wolfe')

    def closure():
        optimizer.zero_grad()
        predicted_reads_exon, predicted_reads_intron, pi = model(design_matrix=dataset_metadata.design_matrix,
                                                                 log_library_sizes=dataset_metadata.log_library_sizes,
                                                                 isoform_length_offset=gene_data.isoform_length_offset,
                                                                 reduced_design_matrix=reduced_matrix)
        loss = pol_2_total_loss(reads_exon=gene_data.exon_reads,
                                reads_introns=gene_data.intron_reads,
                                coverage=gene_data.coverage,
                                predicted_reads_exon=predicted_reads_exon,
                                predicted_reads_intron=predicted_reads_intron,
                                pi=pi)
        if not torch.isfinite(loss):
            return torch.as_tensor(LOSS_CLAMP_VALUE, dtype=loss.dtype, device=loss.device)
        loss.backward()
        return loss

    def evaluate_loss():
        with torch.no_grad():
            predicted_reads_exon, predicted_reads_intron, pi = model(design_matrix=dataset_metadata.design_matrix,
                                                                     log_library_sizes=dataset_metadata.log_library_sizes,
                                                                     isoform_length_offset=gene_data.isoform_length_offset,
                                                                     reduced_design_matrix=reduced_matrix)
            return pol_2_total_loss(reads_exon=gene_data.exon_reads,
                                    reads_introns=gene_data.intron_reads,
                                    coverage=gene_data.coverage,
                                    predicted_reads_exon=predicted_reads_exon,
                                    predicted_reads_intron=predicted_reads_intron,
                                    pi=pi)

    previous_loss = None
    patience_counter = 0

    best_loss = np.inf
    best_state_dict = {k: v.detach().clone() for k, v in model.state_dict().items()}

    converged_within_max_epochs = False
    training_diverged = False

    losses: list[float] = []
    for epoch in range(max_epochs):
        optimizer.step(closure)
        loss_tensor = evaluate_loss()
        loss = loss_tensor.item()
        losses.append(loss)
        if not torch.isfinite(loss_tensor):
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
                                                            dataset_metadata.log_library_sizes,
                                                            gene_data.isoform_length_offset))
        predicted_reads_exon, predicted_reads_intron, pi = outputs
        return pol_2_total_loss(reads_exon=gene_data.exon_reads,
                                reads_introns=gene_data.intron_reads,
                                coverage=gene_data.coverage,
                                predicted_reads_exon=predicted_reads_exon,
                                predicted_reads_intron=predicted_reads_intron,
                                pi=pi)

    hessian_dict = hessian(loss_by_model_parameters)(model_parameters)
    hessian_matrix = flatten_hessian_dict(hessian_dict, model_parameters)
    return hessian_matrix


def add_wald_test_results(df_param: pd.DataFrame, hessian_matrix: torch.Tensor) -> pd.DataFrame:
    """
    Adds results of the Wald test to the dataframe with model parameters.
    """
    df_param = df_param.copy()
    if torch.isfinite(hessian_matrix).all().item():
        df_param['loss_hessian_is_finite'] = True
    else:
        df_param['SE'] = np.nan
        df_param['z_score'] = np.nan
        df_param['p_value_wald'] = np.nan
        df_param['identifiable'] = False
        df_param['loss_hessian_is_finite'] = False
        return df_param

    rank = torch.linalg.matrix_rank(hessian_matrix).item()

    # Rank-revealing QR is currently not available in Pytorch (see https://github.com/pytorch/pytorch/issues/10454),
    # so we are going to use SciPy implementation.
    _, _, qr_pivots = scipy.linalg.qr(hessian_matrix.cpu().numpy(), pivoting=True)
    identifiable_indices = np.sort(qr_pivots[:rank])
    hessian_subset = hessian_matrix[identifiable_indices][:, identifiable_indices]

    if not set(qr_pivots) == set(df_param.index):
        raise ValueError(f'Unexpected value of df_param.index: {df_param=}')

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
                         intron_specific_lfc: bool,
                         device='cpu'
                         ) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_data = gene_data.to(device)
    dataset_metadata = dataset_metadata.to(device)
    pol_2_total_loss = Pol2TotalLoss().to(device)

    model_full = Pol2Model(feature_names=dataset_metadata.feature_names,
                           intron_names=gene_data.intron_names,
                           intron_specific_lfc=intron_specific_lfc).to(device)
    model_full.initialize_intercepts(gene_data, dataset_metadata.library_sizes)

    model_full, training_results_full = train_model(model=model_full,
                                                    gene_data=gene_data,
                                                    dataset_metadata=dataset_metadata,
                                                    pol_2_total_loss=pol_2_total_loss,
                                                    device=device)

    model_param_df = model_full.get_param_df()
    model_param_df['gene_name'] = gene_data.gene_name
    model_param_df['loss_full_model'] = training_results_full.final_loss
    model_param_df['training_diverged_full_model'] = training_results_full.training_diverged
    model_param_df[
        'training_converged_within_max_epochs_full_model'] = training_results_full.converged_within_max_epochs

    # Wald test
    fisher_information_matrix = get_fisher_information_matrix(model=model_full,
                                                              pol_2_total_loss=pol_2_total_loss,
                                                              gene_data=gene_data,
                                                              dataset_metadata=dataset_metadata)
    model_param_df = add_wald_test_results(model_param_df, fisher_information_matrix)
    model_param_df = model_param_df.set_index(["parameter_type", "feature_name", "intron_name"], drop=False)

    test_results_list: list[dict] = []
    for _, lrt_metadata_row in dataset_metadata.lrt_metadata.iterrows():
        reduced_design_matrix = dataset_metadata.reduced_matrices[lrt_metadata_row['test_name']]
        for tested_parameter in TestableParameters:
            intron_names = gene_data.intron_names if intron_specific_lfc else [None]
            for intron_name in intron_names:
                lrt_specification = LRTSpecification(num_features_reduced_matrix=reduced_design_matrix.shape[1],
                                                     tested_parameter=tested_parameter,
                                                     tested_intron=intron_name)

                model_reduced = Pol2Model(feature_names=dataset_metadata.feature_names,
                                          intron_names=gene_data.intron_names,
                                          intron_specific_lfc=intron_specific_lfc,
                                          lrt_specification=lrt_specification).to(device)

                # TODO: hot-start the restricted model from the unrestricted model
                state_dict_model_full = model_full.state_dict()
                hot_start_state_dict = model_reduced.state_dict()
                for key, value in state_dict_model_full.items():
                    hot_start_state_dict[key] = value
                    
                model_reduced.load_state_dict(hot_start_state_dict)
                model_reduced, training_results_reduced = train_model(model=model_reduced,
                                                                      gene_data=gene_data,
                                                                      dataset_metadata=dataset_metadata,
                                                                      pol_2_total_loss=pol_2_total_loss,
                                                                      device=device,
                                                                      reduced_matrix_name=lrt_metadata_row[
                                                                          'test_name'])

                test_result = lrt_metadata_row.to_dict()
                test_result['tested_parameter'] = lrt_specification.tested_parameter
                test_result['tested_intron'] = lrt_specification.tested_intron

                test_result['training_diverged_reduced_model'] = training_results_reduced.training_diverged
                test_result[
                    'training_converged_within_max_epochs_reduced_model'] = training_results_reduced.converged_within_max_epochs

                lfc_value_positive = 0 if pd.isna(test_result['lfc_column_positive']) else model_param_df.loc[
                    (lrt_specification.tested_parameter, test_result['lfc_column_positive'],
                     lrt_specification.tested_intron)]['value']

                lfc_value_negative = 0 if pd.isna(test_result['lfc_column_negative']) else model_param_df.loc[
                    (lrt_specification.tested_parameter, test_result['lfc_column_negative'],
                     lrt_specification.tested_intron)]['value']

                test_result['lfc'] = lfc_value_positive - lfc_value_negative
                test_result['loss_full_model'] = training_results_full.final_loss
                test_result['loss_reduced_model'] = training_results_reduced.final_loss
                test_result['chi_2_test_statistics'] = 2 * (
                        training_results_reduced.final_loss - training_results_full.final_loss)
                test_result['p_value_lrt'] = 1 - stats.chi2.cdf(test_result['chi_2_test_statistics'],
                                                            df=test_result['lrt_df'])

                test_results_list.append(test_result)

    test_results_df = pd.DataFrame(test_results_list)
    return model_param_df, test_results_df

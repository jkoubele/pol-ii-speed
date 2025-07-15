from pathlib import Path
import pickle
import torch
from torch import nn
from torch import optim
import pandas as pd
import numpy as np

from tqdm import trange, tqdm
from dataclasses import dataclass
from pol_ii_model_refactor import GeneData, DatasetMetadata, Pol2TotalLoss, Pol2Model
from torch.func import functional_call, hessian
from time import time
from scipy import stats

import scipy

from concurrent.futures import ProcessPoolExecutor, as_completed


def train_model(gene_data: GeneData,
                dataset_metadata: DatasetMetadata,
                pol_2_total_loss: Pol2TotalLoss,
                device: str,
                max_epochs=100,
                max_patience=5,
                loss_change_tolerance=1e-6
                ) -> tuple[Pol2Model, float]:
    model = Pol2Model(feature_names=feature_names, intron_names=gene_data.intron_names).to(device)
    model.initialize_intercepts(gene_data, library_sizes)

    optimizer = optim.LBFGS(model.parameters(),
                            lr=1.0,
                            max_iter=100,
                            tolerance_change=1e-09,
                            tolerance_grad=1e-07,
                            history_size=100,
                            line_search_fn=None)

    def closure():
        optimizer.zero_grad()
        predicted_log_reads_exon, predicted_reads_intron, phi = model(dataset_metadata.design_matrix,
                                                                      dataset_metadata.log_library_sizes)
        loss = pol_2_total_loss(reads_exon=gene_data.exon_reads,
                                reads_introns=gene_data.intron_reads,
                                coverage=gene_data.coverage,
                                predicted_log_reads_exon=predicted_log_reads_exon,
                                predicted_reads_intron=predicted_reads_intron,
                                phi=phi)
        loss.backward()
        return loss

    previous_loss = None
    patience_counter = 0

    for epoch in range(max_epochs):
        loss = optimizer.step(closure).item()

        if previous_loss is not None:
            relative_loss_change = abs(loss - previous_loss) / (abs(previous_loss) + 1e-10)
            if relative_loss_change < loss_change_tolerance:
                patience_counter += 1
                if patience_counter >= max_patience:
                    break
            else:
                patience_counter = 0

        previous_loss = loss

    return model, loss


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
        predicted_log_reads_exon, predicted_reads_intron, phi = outputs
        return pol_2_total_loss(reads_exon=gene_data.exon_reads,
                                reads_introns=gene_data.intron_reads,
                                coverage=gene_data.coverage,
                                predicted_log_reads_exon=predicted_log_reads_exon,
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
                         perform_wald_test=True) -> pd.DataFrame:
    gene_data = gene_data.to(device)
    dataset_metadata = dataset_metadata.to(device)
    pol_2_total_loss = Pol2TotalLoss().to(device)

    model, loss = train_model(gene_data=gene_data,
                              dataset_metadata=dataset_metadata,
                              pol_2_total_loss=pol_2_total_loss,
                              device=device)

    param_df = model.get_param_df()
    param_df['gene_name'] = gene_data.gene_name
    if perform_wald_test:
        fisher_information_matrix = get_fisher_information_matrix(model=model,
                                                                  pol_2_total_loss=pol_2_total_loss,
                                                                  gene_data=gene_data,
                                                                  dataset_metadata=dataset_metadata)

        param_df = add_wald_test_results(param_df, fisher_information_matrix)
    return param_df


if __name__ == "__main__":
    # Load training data
    input_folder = Path('/home/jakub/Desktop/dev-pol-ii-analysis/data/drosophila_mutants/refactored_train_data/')
    with open(input_folder / 'train_data.pkl', 'rb') as file:
        train_data = pickle.load(file)

    gene_data_list = [GeneData(**gene_data_dict) for gene_data_dict in train_data]

    # Load or create dataset metadata    
    design_matrix_df = pd.DataFrame(data={'genotype_rp2': 6 * [0] + 6 * [1],
                                          'dummy': 3 * [1] + 6 * [0] + 3 * [1]})
    # design_matrix_df = pd.DataFrame(data={'genotype_rp2': 6 * [0] + 6 * [1]})    

    feature_names = design_matrix_df.columns.to_list()
    design_matrix = torch.tensor(design_matrix_df.to_numpy()).float()

    library_sizes = torch.ones(len(design_matrix))

    dataset_metadata = DatasetMetadata(design_matrix=design_matrix,
                                       library_sizes=library_sizes,
                                       log_library_sizes=torch.log(library_sizes),
                                       feature_names=feature_names)

    device = 'cpu'
    # %%
    # all_results = []
    # for gene_data in tqdm(gene_data_list):
    #     all_results.append(get_results_for_gene(gene_data=gene_data, dataset_metadata=dataset_metadata))

    # %%

    gene_data = gene_data_list[0]
    gene_data = gene_data.to(device)
    dataset_metadata = dataset_metadata.to(device)

    pol_2_total_loss = Pol2TotalLoss().to(device)

    model, loss = train_model(gene_data=gene_data,
                              dataset_metadata=dataset_metadata,
                              pol_2_total_loss=pol_2_total_loss,
                              device=device)

    param_df = model.get_param_df()
    fisher_information_matrix = get_fisher_information_matrix(model=model,
                                                              pol_2_total_loss=pol_2_total_loss,
                                                              gene_data=gene_data,
                                                              dataset_metadata=dataset_metadata)

    param_df = add_wald_test_results(param_df, fisher_information_matrix)

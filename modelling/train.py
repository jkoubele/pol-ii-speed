import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scipy
import torch
from scipy import stats
from torch import optim
from torch.func import functional_call, hessian
from tqdm import tqdm, trange

from pol_ii_model import GeneData, DatasetMetadata, Pol2TotalLoss, Pol2Model


@dataclass
class TrainingResults:
    final_loss: float
    converged_within_max_epochs: bool
    num_epochs: int
    losses: list[float]
    gradient_norms: list[float]
    post_adam_loss: Optional[float] = None


def train_model(gene_data: GeneData,
                dataset_metadata: DatasetMetadata,
                pol_2_total_loss: Pol2TotalLoss,
                device: str,
                max_epochs=100,
                max_patience=5,
                loss_change_tolerance=1e-6,
                model: Optional[Pol2Model] = None,
                adam_steps=1
                ) -> tuple[Pol2Model, TrainingResults]:
    if model is None:
        model = Pol2Model(feature_names=feature_names, intron_names=gene_data.intron_names).to(device)
        model.initialize_intercepts(gene_data, library_sizes)

    optimizer = optim.LBFGS(model.parameters(),
                            lr=1.0,
                            max_iter=50,
                            tolerance_change=1e-09,
                            tolerance_grad=1e-07,
                            history_size=100,
                            line_search_fn='strong_wolfe') # 'strong_wolfe'
    
    # optimizer = optim.Adam(model.parameters(),lr=1e-2)

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
        
        grad_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1e20)
        closure.grad_norm = grad_norm  # save to function attribute
    
        return loss
    
    adam_optimizer = optim.Adam(model.parameters(),lr=1e-2)
    for _ in trange(adam_steps):  # You can adjust 200 to 100–500
        adam_optimizer.zero_grad()
        predicted_reads_exon, predicted_reads_intron, phi = model(dataset_metadata.design_matrix,
                                                                   dataset_metadata.log_library_sizes)
        loss = pol_2_total_loss(reads_exon=gene_data.exon_reads,
                                reads_introns=gene_data.intron_reads,
                                coverage=gene_data.coverage,
                                predicted_reads_exon=predicted_reads_exon,
                                predicted_reads_intron=predicted_reads_intron,
                                phi=phi)
        loss.backward()
        adam_optimizer.step()
    
    post_adam_loss = loss.item()

    previous_loss = None
    patience_counter = 0

    converged = False

    losses: list[float] = []
    gradient_norms: list[float] = []
    for epoch in range(max_epochs):
        loss = optimizer.step(closure).item()
        grad_norm = closure.grad_norm  # retrieve from closure

        losses.append(loss)
        gradient_norms.append(grad_norm)
        if np.isnan(loss):
            break

        if previous_loss is not None:
            relative_loss_change = abs(loss - previous_loss) / (abs(previous_loss) + 1e-10)
            if relative_loss_change < loss_change_tolerance:
                patience_counter += 1
                if patience_counter >= max_patience:
                    converged = True                    
                    break
            else:
                patience_counter = 0

        previous_loss = loss
        
    
    
    
    training_results = TrainingResults(final_loss=loss,
                                       converged_within_max_epochs=converged,
                                       num_epochs=epoch + 1,
                                       losses=losses,
                                       gradient_norms=gradient_norms,
                                       post_adam_loss=post_adam_loss)

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
                         perform_lrt=False) -> pd.DataFrame:
    gene_data = gene_data.to(device)
    dataset_metadata = dataset_metadata.to(device)
    pol_2_total_loss = Pol2TotalLoss().to(device)

    model, training_results = train_model(gene_data=gene_data,
                                          dataset_metadata=dataset_metadata,
                                          pol_2_total_loss=pol_2_total_loss,
                                          device=device)

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

        for _, row in tqdm(param_df.iterrows()):
            if row['parameter_type'] not in ('alpha', 'beta', 'gamma'):
                loss_differences.append(None)
            else:
                model_restricted = Pol2Model(feature_names=feature_names, intron_names=gene_data.intron_names).to(
                    device)
                model_restricted.load_state_dict(model.state_dict())

                model_restricted.set_parameter_mask(param_name=row['parameter_type'],
                                                    feature_name=row['feature_name'],
                                                    intron_name=None if row['parameter_type'] == 'alpha' else row[
                                                        'intron_name'],
                                                    value=0.0)
                model_restricted, training_results = train_model(gene_data=gene_data,
                                                                 dataset_metadata=dataset_metadata,
                                                                 pol_2_total_loss=pol_2_total_loss,
                                                                 device=device,
                                                                 model=model_restricted)
                if training_results.converged_within_max_epochs:
                    loss_differences.append(training_results.final_loss - loss_unrestricted)
                else:
                    loss_differences.append(None)

                    # if training_results.final_loss > 1e6:
                #     assert False, 'Model diverged'
                #     break

        
        param_df['loss_differences'] = loss_differences
        param_df['loss_restricted'] = param_df['loss_unrestricted'] + param_df['loss_differences']
        param_df['p_value_lrt'] = 1 - stats.chi2.cdf(2 * param_df['loss_differences'], df=1)
    return param_df


if __name__ == "__main__":
    # Load training data
    input_folder = Path('/cellfile/datapublic/jkoubele/data_pol_ii/drosophila_mutants/train_data')
    with open(input_folder / 'gene_data.pkl', 'rb') as file:
        gene_data_list = pickle.load(file)

    # gene_data_list = [GeneData(**gene_data_dict) for gene_data_dict in train_data]

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
    # for gene_data in tqdm(gene_data_list[:5]):
    #     all_results.append(get_results_for_gene(gene_data=gene_data, 
    #                                             dataset_metadata=dataset_metadata,
    #                                             perform_lrt=True))
    # df_all = pd.concat(all_results)

    # %%
    
    # gene_name	FBgn0000120 - large loss but may be ok (32k loss) - is ok as the read counts is quite high
	# gene_name FBgn0000250 - negative loss (thats ok since continuous NLL can be both positive and negative)
    # gene_name FBgn0000036 diverging training


    gene_data = [x for x in gene_data_list if x.gene_name == 'FBgn0000036'][0]
    gene_data = gene_data.to(device)
    dataset_metadata = dataset_metadata.to(device)

    pol_2_total_loss = Pol2TotalLoss().to(device)

    model, training_results = train_model(gene_data=gene_data,
                                          dataset_metadata=dataset_metadata,
                                          pol_2_total_loss=pol_2_total_loss,
                                          device=device,
                                          adam_steps=1)

    print(f"{training_results=}")
    

    param_df = model.get_param_df()
    
 
    
    fisher_information_matrix = get_fisher_information_matrix(model=model,
                                                              pol_2_total_loss=pol_2_total_loss,
                                                              gene_data=gene_data,
                                                              dataset_metadata=dataset_metadata)

    param_df = add_wald_test_results(param_df, fisher_information_matrix)

    # Logic for LRT

    loss_unrestricted = training_results.final_loss

    loss_differences: list[Optional[float]] = []
    for _, row in tqdm(param_df.iterrows()):
        if row['parameter_type'] not in ('alpha', 'beta', 'gamma'):
            loss_differences.append(None)
        else:
            model_restricted = Pol2Model(feature_names=feature_names, intron_names=gene_data.intron_names).to(device)
            model_restricted.load_state_dict(model.state_dict())

            model_restricted.set_parameter_mask(param_name=row['parameter_type'],
                                                feature_name=row['feature_name'],
                                                intron_name=None if row['parameter_type'] == 'alpha' else row[
                                                    'intron_name'],
                                                value=0.0)
            model_restricted, training_results = train_model(gene_data=gene_data,
                                                             dataset_metadata=dataset_metadata,
                                                             pol_2_total_loss=pol_2_total_loss,
                                                             device=device,
                                                             model=model_restricted)
            if training_results.converged_within_max_epochs:
                loss_differences.append(training_results.final_loss - loss_unrestricted)
            else:
                loss_differences.append(None)

            if training_results.final_loss > 1e6:
                assert False, 'Model diverged'
                break

    param_df['loss_unrestricted'] = loss_unrestricted
    param_df['loss_differences'] = loss_differences
    param_df['loss_restricted'] = param_df['loss_unrestricted'] + param_df['loss_differences']
    param_df['p_value_lrt'] = 1 - stats.chi2.cdf(2 * param_df['loss_differences'], df=1)

import pickle
from pathlib import Path

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import scipy
import torch
import torch.nn as nn
import torch.optim as optim
from scipy import stats
from torch.func import functional_call, hessian
from tqdm import trange
from scipy.stats import chi2

from typing import Optional
from pol_ii_model_old import Pol2Model, GeneData, Pol2TotalLoss, ParameterMask
from read_location_model import estimate_phi

HessianDict = dict[str, dict[str, torch.Tensor]]  # structure outputed by Pytorch hessian() function


def train_model(gene_data: GeneData,
                X: torch.Tensor,
                library_sizes: torch.Tensor,
                pol_2_total_loss: Pol2TotalLoss,
                device: str,
                max_epochs=100,
                mask_alpha: Optional[ParameterMask] = None,
                mask_beta: Optional[ParameterMask] = None,
                mask_gamma: Optional[ParameterMask] = None
                ) -> Pol2Model:
    log_library_size = torch.log(library_sizes).to(device)
    model = Pol2Model(num_features=X.shape[1],
                      intron_names=gene_data.intron_names,
                      exon_intercept_init=float(torch.log((gene_data.exon_reads / library_sizes).mean())),
                      intron_intercepts_init={intron_name: float(torch.log((intron_reads / library_sizes).mean()))
                                              for intron_name, intron_reads in gene_data.intron_reads.items()},
                      mask_alpha=mask_alpha,
                      mask_beta=mask_beta,
                      mask_gamma=mask_gamma)
    model = model.to(device)
    pol_2_total_loss = pol_2_total_loss.to(device)
    gene_data = gene_data.to(device)

    optimizer = optim.LBFGS(model.parameters(),
                            lr=1.0,
                            max_iter=100,
                            tolerance_change=1e-09,
                            tolerance_grad=1e-07,
                            history_size=100,
                            line_search_fn=None)

    # optimizer = optim.Adam(model.parameters(),lr=1e-1)
    # When using Adam, set max_epochs to about 1-10k

    def closure():
        optimizer.zero_grad()
        predicted_log_reads_exon, predicted_reads_intron, phi = model(X, log_library_size)
        loss = pol_2_total_loss(gene_data, predicted_log_reads_exon, predicted_reads_intron, phi)
        loss.backward()
        return loss

    for epoch in trange(max_epochs):
        loss = optimizer.step(closure)
    return model


def fit_analytical_solution(gene_data: GeneData,
                            X: torch.Tensor,
                            library_sizes: torch.Tensor,
                            device: str
                            ) -> Pol2Model:
    assert X.shape[1] == 1
    assert set(X.flatten().tolist()) == {0, 1}
    threshold_phi = 0.001  # to assert that phi lies in (threshold_phi, 1-threshold_phi)

    model_analytical_fit = Pol2Model(num_features=1,
                                     intron_names=gene_data.intron_names)

    mask_group_1 = X.flatten() == 0
    mask_group_2 = X.flatten() == 1

    mu_1 = gene_data.exon_reads[mask_group_1].sum() / library_sizes[mask_group_1].sum()
    mu_2 = gene_data.exon_reads[mask_group_2].sum() / library_sizes[mask_group_2].sum()

    assert mu_1 > 0 and mu_2 > 0

    model_analytical_fit.intercept_exon = nn.Parameter(torch.log(mu_1))
    model_analytical_fit.alpha = nn.Parameter(torch.log(mu_2 / mu_1).unsqueeze(0))

    for intron_name in gene_data.intron_names:
        intron_reads = gene_data.intron_reads[intron_name]
        coverage_density = gene_data.coverage_density[intron_name]

        nu_1 = intron_reads[mask_group_1].sum() / library_sizes[mask_group_1].sum()
        nu_2 = intron_reads[mask_group_2].sum() / library_sizes[mask_group_2].sum()

        assert nu_1 > 0 and nu_2 > 0

        coverage_group_1 = coverage_density[mask_group_1]
        aggregate_density_1 = (coverage_group_1 * intron_reads[mask_group_1].unsqueeze(1)).sum(dim=0)
        aggregate_density_1 /= aggregate_density_1.sum()

        coverage_group_2 = coverage_density[mask_group_2]
        aggregate_density_2 = (coverage_group_2 * intron_reads[mask_group_2].unsqueeze(1)).sum(dim=0)
        aggregate_density_2 /= aggregate_density_2.sum()

        phi_1 = torch.tensor(estimate_phi(aggregate_density_1.numpy()), dtype=torch.float32)
        phi_2 = torch.tensor(estimate_phi(aggregate_density_2.numpy()), dtype=torch.float32)

        phi_1 = torch.clip(phi_1, min=threshold_phi, max=1 - threshold_phi)
        phi_2 = torch.clip(phi_2, min=threshold_phi, max=1 - threshold_phi)

        model_analytical_fit.intercept_intron[intron_name] = nn.Parameter(torch.log((1 - phi_1) * nu_1))
        model_analytical_fit.log_phi_zero[intron_name] = nn.Parameter(torch.log(phi_1 / (1 - phi_1)))

        model_analytical_fit.beta[intron_name] = nn.Parameter(
            (torch.log(mu_2 / mu_1) - torch.log(nu_2 / nu_1) - torch.log(phi_2 / phi_1)).unsqueeze(0))

        model_analytical_fit.gamma[intron_name] = nn.Parameter(
            (-torch.log(mu_2 / mu_1) + torch.log(nu_2 / nu_1) + torch.log((1 - phi_2) / (1 - phi_1))).unsqueeze(0))

    model_analytical_fit = model_analytical_fit.to(device)
    return model_analytical_fit


def get_param_names_and_sizes(model_parameters: dict[str, nn.Parameter]
                              ) -> tuple[list[str], list[int]]:
    parameter_names: list[str] = []
    parameter_sizes: list[int] = []
    for name, parameter in model_parameters.items():
        parameter_names.append(name)
        parameter_sizes.append(parameter.numel())
    return parameter_names, parameter_sizes


def hessian_matrix_from_hessian_dict(hessian_dict: HessianDict,
                                     parameter_names: list[str],
                                     parameter_sizes: list[int],
                                     device: str) -> torch.Tensor:
    # Flatten the Hessian into a 2D matrix
    index_map = {}
    start = 0
    for name, size in zip(parameter_names, parameter_sizes):
        index_map[name] = (start, start + size)
        start += size
    total_params = start

    hessian_matrix = torch.zeros((total_params, total_params), device=device)

    for name1, block_row in hessian_dict.items():
        idx1_start, idx1_end = index_map[name1]
        for name2, block in block_row.items():
            idx2_start, idx2_end = index_map[name2]
            hessian_matrix[idx1_start:idx1_end, idx2_start:idx2_end] = block.detach().reshape(
                idx1_end - idx1_start, idx2_end - idx2_start
            )

    hessian_matrix = hessian_matrix.detach()
    return hessian_matrix


def compute_hessian_matrix(model: Pol2Model,
                           pol_2_total_loss: Pol2TotalLoss,
                           X: torch.Tensor,
                           log_library_size: torch.Tensor,
                           gene_data: GeneData,
                           device: str) -> torch.Tensor:
    model_parameters = dict(model.named_parameters())
    parameter_names, parameter_sizes = get_param_names_and_sizes(model_parameters)

    def loss_by_model_parameters(model_parameters):
        outputs = functional_call(model, model_parameters, (X, log_library_size))
        predicted_log_reads_exon, predicted_reads_intron, phi = outputs
        return pol_2_total_loss(gene_data, predicted_log_reads_exon, predicted_reads_intron, phi)

    hessian_dict = hessian(loss_by_model_parameters)(model_parameters)
    hessian_matrix = hessian_matrix_from_hessian_dict(hessian_dict, parameter_names, parameter_sizes, device)
    return hessian_matrix


def get_param_df(model: Pol2Model,
                 feature_names: list[str]) -> pd.DataFrame:
    model_parameters = dict(model.named_parameters())
    parameter_names, parameter_sizes = get_param_names_and_sizes(model_parameters)
    parameter_data: list[dict] = []
    for name, size in zip(parameter_names, parameter_sizes):
        if '.' in name:
            parameter_type = name.split('.', 1)[0]
            intron_name = name.split('.', 1)[1]
        else:
            parameter_type = name.split('.', 1)[0]
            intron_name = None

        parameter_value_list = [model_parameters[name].item()] if model_parameters[name].numel() == 1 \
            else model_parameters[name].tolist()
        for feature_index in range(size):
            feature_name = None
            if parameter_type in ('alpha', 'beta', 'gamma'):
                feature_name = feature_names[feature_index]
            parameter_data.append({'parameter_type': parameter_type,
                                   'intron_name': intron_name,
                                   'feature_name': feature_name,
                                   'value': parameter_value_list[feature_index]})

    df_param = pd.DataFrame(data=parameter_data)
    return df_param


def compare_model_parameters(model_analytical: Pol2Model, model_numerical: Pol2Model, rtol=1e-2) -> None:
    analytical_parameters = dict(model_analytical.named_parameters())
    numerical_parameters = dict(model_numerical.named_parameters())
    assert analytical_parameters.keys() == numerical_parameters.keys()

    for name, param_analytical in analytical_parameters.items():
        param_numerical = numerical_parameters[name]
        if not torch.allclose(param_numerical, param_analytical, rtol=rtol):
            print("------OUT OF TOLERANCE--------")
            print(f"{name=}")
            print(f"{param_numerical=}")
            print(f"{param_analytical=}")
            abs_diff = torch.abs(param_numerical - param_analytical)
            rel_diff = abs_diff / param_analytical
            print(f"{abs_diff=}")
            print(f"{rel_diff=}")


def compare_model_predictions(model_analytical: Pol2Model, 
                              model_numerical: Pol2Model,
                              X:torch.Tensor, 
                              log_library_size: torch.Tensor,
                              rtol=1e-3) -> None:
    predicted_log_reads_exon_numerical, predicted_reads_intron_numerical, phi_numerical = model_numerical(X,
                                                                                                          log_library_size)
    predicted_log_reads_exon_analytical, predicted_reads_intron_analytical, phi_analytical = model_analytical(X,
                                                                                                              log_library_size)

    if not torch.allclose(predicted_log_reads_exon_numerical, predicted_log_reads_exon_analytical, rtol=rtol):
        print(50 * "-")
        print(f"{predicted_log_reads_exon_numerical=}")
        print(f"{predicted_log_reads_exon_analytical=}")

    assert predicted_reads_intron_numerical.keys() == predicted_reads_intron_analytical.keys()
    for intron_name in predicted_reads_intron_numerical.keys():
        reads_numerical = predicted_reads_intron_numerical[intron_name]
        reads_analytical = predicted_reads_intron_analytical[intron_name]

        if not torch.allclose(reads_numerical, reads_analytical, rtol=rtol):
            print(50 * "-")
            print(f"{reads_numerical=}")
            print(f"{reads_analytical=}")

        phi_intron_numerical = phi_numerical[intron_name]
        phi_intron_analytical = phi_analytical[intron_name]
        if not torch.allclose(phi_intron_numerical, phi_intron_analytical, rtol=rtol):
            print(50 * "-")
            print(f"{phi_intron_numerical=}")
            print(f"{phi_intron_analytical=}")


def add_wald_test_results(df_param: pd.DataFrame, hessian_matrix: torch.Tensor) -> pd.DataFrame:
    """
    Adds results of the Wald test to the dataframe with model parameters.
    """
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


# %%
if __name__ == "__main__":
    project_path = Path('/home/jakub/Desktop/dev-pol-ii-analysis/data/drosophila_mutants')
    # project_path = Path('/cellfile/datapublic/jkoubele/data_pol_ii/drosophila_mutants/')
    with open(project_path / 'train_data' / 'data_train_with_solutions.pkl', 'rb') as file:
        data_train_with_solutions = pickle.load(file)

    gene_data = data_train_with_solutions['gene_data_with_solutions'][0].gene_data

    design_matrix = pd.DataFrame(data={'genotype_rp2': 6 * [0] + 6 * [1],
                                       'dummy': 3 * [1] + 6 * [0] + 3 * [1]})

    design_matrix = pd.DataFrame(data={'genotype_rp2': 6 * [0] + 6 * [1]})
    
    use_cuda = True
    # Above is the input to the 'main' per-gene function

    feature_names = design_matrix.columns.to_list()
    device = torch.device("cuda" if torch.cuda.is_available() and use_cuda else "cpu")
    device = 'cpu'

    X = torch.tensor(design_matrix.to_numpy()).float().to(device)
    library_sizes = torch.ones(len(design_matrix)).float().to(device)
    library_sizes /= library_sizes.max()  # This would be good to move to data preparation maybe
    log_library_size = torch.log(library_sizes).to(device)
    pol_2_total_loss = Pol2TotalLoss().to(device)

    model_numerical = train_model(gene_data, X, library_sizes, pol_2_total_loss, device)
    df_param_numerical = get_param_df(model_numerical, feature_names)

    hessian_matrix_numerical = compute_hessian_matrix(model=model_numerical,
                                                      pol_2_total_loss=pol_2_total_loss,
                                                      X=X,
                                                      log_library_size=log_library_size,
                                                      gene_data=gene_data,
                                                      device=device)

    df_param_numerical = add_wald_test_results(df_param_numerical, hessian_matrix_numerical)
    
    # LRT
    predicted_log_reads_exon, predicted_reads_intron, phi = model_numerical(X, log_library_size)
    loss_unrestricted = pol_2_total_loss(gene_data, predicted_log_reads_exon, predicted_reads_intron, phi)    
    df_param_numerical['LRT_statistics'] = None
    df_param_numerical['p_value_LRT'] = None
    
    #TODO: Separate LRT to a function
    for feature_index, feature_name in enumerate(feature_names):        
        # Test for alpha
        model_restricted = train_model(gene_data, X, library_sizes, pol_2_total_loss, device,
                                       mask_alpha=ParameterMask(feature_index=feature_index, value=0.0))    
        predicted_log_reads_exon_restricted, predicted_reads_intron_restricted, phi_restricted = model_restricted(X, log_library_size)
        loss_restricted = pol_2_total_loss(gene_data,
                                           predicted_log_reads_exon_restricted, 
                                           predicted_reads_intron_restricted,
                                           phi_restricted)
        lrt_statistics = 2 * (loss_restricted-loss_unrestricted).item()
        p_value_lrt = 1 - chi2.cdf(lrt_statistics, df=1)
        
        dataframe_row_mask = (df_param_numerical['parameter_type'] == 'alpha') & (df_param_numerical['feature_name'] == feature_name)
        df_param_numerical.loc[dataframe_row_mask, 'LRT_statistics'] = lrt_statistics
        df_param_numerical.loc[dataframe_row_mask, 'p_value_LRT'] = p_value_lrt
        
        for intron_name in gene_data.intron_names:
            # Test for beta
            model_restricted = train_model(gene_data, X, library_sizes, pol_2_total_loss, device,
                                           mask_beta=ParameterMask(feature_index=feature_index, value=0.0, intron_name=intron_name))    
            predicted_log_reads_exon_restricted, predicted_reads_intron_restricted, phi_restricted = model_restricted(X, log_library_size)
            loss_restricted = pol_2_total_loss(gene_data,
                                               predicted_log_reads_exon_restricted, 
                                               predicted_reads_intron_restricted,
                                               phi_restricted)
            lrt_statistics = 2 * (loss_restricted-loss_unrestricted).item()
            p_value_lrt = 1 - chi2.cdf(lrt_statistics, df=1)
            
            dataframe_row_mask = (df_param_numerical['parameter_type'] == 'beta') & (df_param_numerical['feature_name'] == feature_name) & (df_param_numerical['intron_name'] == intron_name)
            df_param_numerical.loc[dataframe_row_mask, 'LRT_statistics'] = lrt_statistics
            df_param_numerical.loc[dataframe_row_mask, 'p_value_LRT'] = p_value_lrt
            
            # Test for gamma
            model_restricted = train_model(gene_data, X, library_sizes, pol_2_total_loss, device,
                                           mask_gamma=ParameterMask(feature_index=feature_index, value=0.0, intron_name=intron_name))    
            predicted_log_reads_exon_restricted, predicted_reads_intron_restricted, phi_restricted = model_restricted(X, log_library_size)
            loss_restricted = pol_2_total_loss(gene_data,
                                               predicted_log_reads_exon_restricted, 
                                               predicted_reads_intron_restricted,
                                               phi_restricted)
            lrt_statistics = 2 * (loss_restricted-loss_unrestricted).item()
            p_value_lrt = 1 - chi2.cdf(lrt_statistics, df=1)
            
            dataframe_row_mask = (df_param_numerical['parameter_type'] == 'gamma') & (df_param_numerical['feature_name'] == feature_name) & (df_param_numerical['intron_name'] == intron_name)
            df_param_numerical.loc[dataframe_row_mask, 'LRT_statistics'] = lrt_statistics
            df_param_numerical.loc[dataframe_row_mask, 'p_value_LRT'] = p_value_lrt
    df_param_numerical['Wald_minus_LRT'] = df_param_numerical['p_value_wald'] - df_param_numerical['p_value_LRT']
    print(f"{p_value_lrt}")
    # %%
    model_analytical = fit_analytical_solution(gene_data, X, library_sizes, device)
    compare_model_parameters(model_analytical=model_analytical, model_numerical=model_numerical)
    compare_model_predictions(model_analytical=model_analytical,
                              model_numerical=model_numerical, X=X, 
                              log_library_size=log_library_size)

    df_param_analytical = get_param_df(model_analytical, feature_names)

    hessian_matrix_analytical = compute_hessian_matrix(model=model_analytical,
                                                       pol_2_total_loss=pol_2_total_loss,
                                                       X=X,
                                                       log_library_size=log_library_size,
                                                       gene_data=gene_data,
                                                       device=device)

    df_param_analytical = add_wald_test_results(df_param_analytical, hessian_matrix_analytical)

import numpy as np
import pandas as pd
import scipy
import torch
from dataclasses import dataclass
from scipy import stats
from torch import nn, optim
from torch.func import functional_call, hessian
from typing import Callable, Optional, OrderedDict

from pol_ii_speed_modeling.pol_ii_model import GeneData, IntronData, DatasetMetadata, Pol2TotalLoss, Pol2Model, \
    TestableParameters, LRTSpecification, SplicingModel, CoverageLoss, \
    GlobalGeneData, GlobalPol2Model

LOSS_CLAMP_VALUE = 1e30

StateDict = OrderedDict[str, torch.Tensor]


@dataclass
class TrainingResults:
    final_loss: float
    converged_within_max_epochs: bool
    num_epochs: int
    losses: list[float]
    training_diverged: bool


@dataclass
class CacheForRegularization:
    training_input_per_gene: list[tuple[GeneData, StateDict]]
    dataset_metadata: DatasetMetadata
    intron_specific_lfc: Optional[bool] = None  # Pol2 model only
    intron_specific_splicing: bool = False


def make_lbfgs_optimizer(model: nn.Module) -> optim.LBFGS:
    return optim.LBFGS(
        model.parameters(),
        lr=1.0,
        max_iter=20,
        tolerance_change=1e-9,
        tolerance_grad=1e-7,
        history_size=100,
        line_search_fn='strong_wolfe',
    )


def train_model(
        model: nn.Module,
        optimizer: optim.Optimizer,
        closure: Callable[[], torch.Tensor],
        evaluate_loss: Callable[[], torch.Tensor],
        max_epochs=200,
        max_patience=3,
        loss_change_tolerance=1e-6,
        verbose: bool = False,
) -> tuple[nn.Module, TrainingResults]:
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
        if verbose:
            print(f'Epoch {epoch + 1}/{max_epochs}: loss={loss:.6f}', flush=True)
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
    training_results = TrainingResults(
        final_loss=loss,
        converged_within_max_epochs=converged_within_max_epochs,
        num_epochs=epoch + 1,
        losses=losses,
        training_diverged=training_diverged,
    )

    return model, training_results


def _make_pol2_closures(
        model: Pol2Model,
        optimizer: optim.Optimizer,
        gene_data: GeneData,
        dataset_metadata: DatasetMetadata,
        reduced_matrix: Optional[torch.Tensor] = None,
        regularization_coefficients_beta: Optional[torch.Tensor] = None,
        regularization_coefficients_gamma: Optional[torch.Tensor] = None,
) -> tuple[Callable[[], torch.Tensor], Callable[[], torch.Tensor]]:
    pol_2_total_loss = Pol2TotalLoss().to(next(model.parameters()).device)

    def _compute_loss() -> torch.Tensor:
        predicted_reads_exon, predicted_reads_intron, pi = model(
            design_matrix=dataset_metadata.design_matrix,
            log_library_sizes=dataset_metadata.log_library_sizes,
            isoform_length_offset=gene_data.isoform_length_offset,
            reduced_design_matrix=reduced_matrix,
        )
        loss = pol_2_total_loss(
            reads_exon=gene_data.exon_reads,
            reads_introns=gene_data.intron_reads,
            coverage=gene_data.coverage,
            predicted_reads_exon=predicted_reads_exon,
            predicted_reads_intron=predicted_reads_intron,
            pi=pi,
        )
        if regularization_coefficients_beta is not None:
            loss += torch.sum(regularization_coefficients_beta * torch.square(model.beta))
        if regularization_coefficients_gamma is not None:
            loss += torch.sum(regularization_coefficients_gamma * torch.square(model.gamma))
        return loss

    def closure() -> torch.Tensor:
        optimizer.zero_grad()
        loss = _compute_loss()
        if not torch.isfinite(loss):
            return torch.as_tensor(LOSS_CLAMP_VALUE, dtype=loss.dtype, device=loss.device)
        loss.backward()
        return loss

    def evaluate_loss() -> torch.Tensor:
        with torch.no_grad():
            return _compute_loss()

    return closure, evaluate_loss


def make_splicing_closures(
        model: SplicingModel,
        optimizer: optim.Optimizer,
        coverage: torch.Tensor,
        design_matrix: torch.Tensor,
        regularization_coefficients_lfc: Optional[torch.Tensor] = None,
) -> tuple[Callable[[], torch.Tensor], Callable[[], torch.Tensor]]:
    coverage_loss_fn = CoverageLoss(num_position_coverage=coverage.shape[2]).to(coverage.device)

    def _compute_loss() -> torch.Tensor:
        loss = coverage_loss_fn(model(design_matrix), coverage)
        if regularization_coefficients_lfc is not None:
            loss += torch.sum(regularization_coefficients_lfc * torch.square(model.lfc))
        return loss

    def closure() -> torch.Tensor:
        optimizer.zero_grad()
        loss = _compute_loss()
        if not torch.isfinite(loss):
            return torch.as_tensor(LOSS_CLAMP_VALUE, dtype=loss.dtype, device=loss.device)
        loss.backward()
        return loss

    def evaluate_loss() -> torch.Tensor:
        with torch.no_grad():
            return _compute_loss()

    return closure, evaluate_loss


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


def get_fisher_information_matrix(model: nn.Module,
                                  loss_fn: Callable[[StateDict], torch.Tensor]) -> torch.Tensor:
    model_parameters = dict(model.named_parameters())
    hessian_dict = hessian(loss_fn)(model_parameters)
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
    standard_errors = torch.sqrt(torch.clamp(torch.diag(cov_matrix), min=0.0)).numpy()

    df_param['SE'] = np.nan
    df_param['z_score'] = np.nan
    df_param['p_value_wald'] = np.nan
    df_param['identifiable'] = df_param.index.isin(identifiable_indices)

    df_param.loc[identifiable_indices, 'SE'] = standard_errors
    df_param.loc[identifiable_indices, 'z_score'] = df_param.loc[identifiable_indices, 'value'] / standard_errors
    df_param.loc[identifiable_indices, 'p_value_wald'] = 2 * (
            1 - stats.norm.cdf(np.abs(df_param.loc[identifiable_indices, 'z_score'])))

    return df_param


def get_model_results(gene_data: GeneData,
                      dataset_metadata: DatasetMetadata,
                      intron_specific_lfc: bool,
                      device='cpu'
                      ) -> tuple[pd.DataFrame, pd.DataFrame, StateDict]:
    gene_data = gene_data.to(device)
    dataset_metadata = dataset_metadata.to(device)

    model_full = Pol2Model(feature_names=dataset_metadata.feature_names,
                           intron_names=gene_data.intron_names,
                           intron_specific_lfc=intron_specific_lfc).to(device)
    model_full.initialize_parameters(gene_data, dataset_metadata.library_sizes, dataset_metadata.design_matrix)

    optimizer_full = make_lbfgs_optimizer(model_full)
    closure_full, evaluate_loss_full = _make_pol2_closures(
        model_full, optimizer_full, gene_data, dataset_metadata,
    )
    model_full, training_results_full = train_model(model_full, optimizer_full, closure_full, evaluate_loss_full)

    model_param_df = model_full.get_param_df()
    model_param_df['gene_name'] = gene_data.gene_name
    model_param_df['loss_full_model'] = training_results_full.final_loss
    model_param_df['training_diverged_full_model'] = training_results_full.training_diverged
    model_param_df[
        'training_converged_within_max_epochs_full_model'] = training_results_full.converged_within_max_epochs

    # Wald test
    pol_2_total_loss = Pol2TotalLoss().to(device)

    def pol2_loss_by_params(params):
        predicted_reads_exon, predicted_reads_intron, pi = functional_call(
            model_full, params,
            (dataset_metadata.design_matrix, dataset_metadata.log_library_sizes, gene_data.isoform_length_offset),
        )
        return pol_2_total_loss(reads_exon=gene_data.exon_reads,
                                reads_introns=gene_data.intron_reads,
                                coverage=gene_data.coverage,
                                predicted_reads_exon=predicted_reads_exon,
                                predicted_reads_intron=predicted_reads_intron,
                                pi=pi)

    fisher_information_matrix = get_fisher_information_matrix(model_full, pol2_loss_by_params)
    model_param_df = add_wald_test_results(model_param_df, fisher_information_matrix)
    model_param_df = model_param_df.set_index(["parameter_type", "feature_name", "intron_name"], drop=False)

    test_results_list: list[dict] = []
    for _, lrt_metadata_row in dataset_metadata.lrt_metadata.iterrows():
        reduced_design_matrix = dataset_metadata.reduced_matrices[lrt_metadata_row['test_id']]
        for tested_parameter in TestableParameters:
            intron_names = gene_data.intron_names if intron_specific_lfc and tested_parameter in (
                TestableParameters.BETA, TestableParameters.GAMMA) else [None]
            for intron_name in intron_names:
                lrt_specification = LRTSpecification(num_features_reduced_matrix=reduced_design_matrix.shape[1],
                                                     tested_parameter=tested_parameter,
                                                     tested_intron=intron_name)

                model_reduced = Pol2Model(feature_names=dataset_metadata.feature_names,
                                          intron_names=gene_data.intron_names,
                                          intron_specific_lfc=intron_specific_lfc,
                                          lrt_specification=lrt_specification).to(device)

                state_dict_model_full = model_full.state_dict()
                hot_start_state_dict = model_reduced.state_dict()
                for key, value in state_dict_model_full.items():
                    hot_start_state_dict[key] = value

                if reduced_design_matrix.shape[1] > 0:
                    if lrt_specification.tested_parameter == TestableParameters.ALPHA:
                        lfc_full_model = state_dict_model_full[lrt_specification.tested_parameter]
                    else:
                        intron_index = 0 if not model_full.intron_specific_lfc else model_full.intron_names.index(
                            lrt_specification.tested_intron)
                        lfc_full_model = state_dict_model_full[lrt_specification.tested_parameter][:, intron_index]
                    # Initialize LFC in reduced model via least squares
                    hot_start_state_dict['reduced_lfc'] = torch.linalg.lstsq(reduced_design_matrix,
                                                                             dataset_metadata.design_matrix @ lfc_full_model).solution

                model_reduced.load_state_dict(hot_start_state_dict)

                optimizer_reduced = make_lbfgs_optimizer(model_reduced)
                closure_reduced, evaluate_loss_reduced = _make_pol2_closures(
                    model_reduced, optimizer_reduced, gene_data, dataset_metadata,
                    reduced_matrix=reduced_design_matrix,
                )
                model_reduced, training_results_reduced = train_model(
                    model_reduced, optimizer_reduced, closure_reduced, evaluate_loss_reduced,
                )

                test_result = lrt_metadata_row.to_dict()
                test_result['tested_parameter'] = lrt_specification.tested_parameter
                test_result['gene_name'] = gene_data.gene_name
                test_result['intron_name'] = lrt_specification.tested_intron

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
                test_result['chi2_test_statistics'] = 2 * (
                        training_results_reduced.final_loss - training_results_full.final_loss)
                test_result['p_value'] = 1 - stats.chi2.cdf(test_result['chi2_test_statistics'],
                                                            df=test_result['lrt_df'])

                test_results_list.append(test_result)

    test_results_df = pd.DataFrame(test_results_list)
    return model_param_df, test_results_df, model_full.state_dict()


def get_regularized_model_results(gene_data: GeneData,
                                  dataset_metadata: DatasetMetadata,
                                  intron_specific_lfc: bool,
                                  hot_start_state_dict: StateDict,
                                  regularization_coefficients_df: pd.DataFrame,
                                  device='cpu'
                                  ) -> pd.DataFrame:
    gene_data = gene_data.to(device)
    dataset_metadata = dataset_metadata.to(device)
    regularization_coefficients_df = regularization_coefficients_df.set_index(['parameter_type', 'feature_name'])

    model_regularized = Pol2Model(feature_names=dataset_metadata.feature_names,
                                  intron_names=gene_data.intron_names,
                                  intron_specific_lfc=intron_specific_lfc).to(device)
    model_regularized.load_state_dict(hot_start_state_dict)

    regularization_coefficients_beta = torch.zeros((len(model_regularized.feature_names), 1)).to(device)
    regularization_coefficients_gamma = torch.zeros((len(model_regularized.feature_names), 1)).to(device)
    for feature_index, feature_name in enumerate(model_regularized.feature_names):
        regularization_coefficients_beta[feature_index] = regularization_coefficients_df.loc[('beta', feature_name)][
            'lambda']
        regularization_coefficients_gamma[feature_index] = regularization_coefficients_df.loc[('gamma', feature_name)][
            'lambda']

    optimizer_regularized = make_lbfgs_optimizer(model_regularized)
    closure_regularized, evaluate_loss_regularized = _make_pol2_closures(
        model_regularized, optimizer_regularized, gene_data, dataset_metadata,
        regularization_coefficients_beta=regularization_coefficients_beta,
        regularization_coefficients_gamma=regularization_coefficients_gamma,
    )
    model_regularized, training_results_regularized = train_model(
        model_regularized, optimizer_regularized, closure_regularized, evaluate_loss_regularized,
    )

    model_param_df = model_regularized.get_param_df()
    model_param_df['gene_name'] = gene_data.gene_name
    model_param_df['loss_regularized_model'] = training_results_regularized.final_loss
    model_param_df['training_diverged_regularized_model'] = training_results_regularized.training_diverged
    model_param_df[
        'training_converged_within_max_epochs_regularized_model'] = training_results_regularized.converged_within_max_epochs
    return model_param_df


def get_splicing_model_results(
        coverage: torch.Tensor,
        dataset_metadata: DatasetMetadata,
        intron_names: list[str],
        device: str = 'cpu',
        verbose: bool = False,
        compute_wald_test: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    coverage = coverage.to(device)
    dataset_metadata = dataset_metadata.to(device)

    model_full = SplicingModel(
        feature_names=dataset_metadata.feature_names,
        intron_names=intron_names,
    ).to(device)
    model_full.initialize_theta(coverage)

    optimizer_full = make_lbfgs_optimizer(model_full)
    closure_full, evaluate_loss_full = make_splicing_closures(
        model_full, optimizer_full, coverage, dataset_metadata.design_matrix,
    )
    model_full, results_full = train_model(model_full, optimizer_full, closure_full, evaluate_loss_full,
                                           verbose=verbose)

    model_param_df = model_full.get_param_df()
    model_param_df['loss_full_model'] = results_full.final_loss
    model_param_df['training_diverged_full_model'] = results_full.training_diverged
    model_param_df['training_converged_within_max_epochs_full_model'] = results_full.converged_within_max_epochs

    coverage_loss_fn = CoverageLoss(num_position_coverage=coverage.shape[2]).to(device)

    def splicing_loss_by_params(params):
        pi = functional_call(model_full, params, (dataset_metadata.design_matrix,))
        return coverage_loss_fn(pi, coverage)

    if compute_wald_test:
        fisher_information_matrix = get_fisher_information_matrix(model_full, splicing_loss_by_params)
        model_param_df = add_wald_test_results(model_param_df, fisher_information_matrix)
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
                full_lfc_contribution = dataset_metadata.design_matrix @ model_full.lfc
                model_reduced.lfc.data.copy_(
                    torch.linalg.lstsq(reduced_matrix, full_lfc_contribution).solution
                )

        optimizer_reduced = make_lbfgs_optimizer(model_reduced)
        closure_reduced, evaluate_loss_reduced = make_splicing_closures(
            model_reduced, optimizer_reduced, coverage, reduced_matrix,
        )
        model_reduced, results_reduced = train_model(
            model_reduced, optimizer_reduced, closure_reduced, evaluate_loss_reduced,
            verbose=verbose,
        )

        lfc_positive = 0.0 if pd.isna(lrt_row['lfc_column_positive']) else \
            model_param_df.loc[('lfc', lrt_row['lfc_column_positive'], None)]['value']
        lfc_negative = 0.0 if pd.isna(lrt_row['lfc_column_negative']) else \
            model_param_df.loc[('lfc', lrt_row['lfc_column_negative'], None)]['value']

        chi2_stat = 2 * (results_reduced.final_loss - results_full.final_loss)
        test_result = lrt_row.to_dict()
        test_result['tested_parameter'] = 'lfc'
        test_result['intron_name'] = None
        test_result['lfc'] = lfc_positive - lfc_negative
        test_result['loss_full_model'] = results_full.final_loss
        test_result['loss_reduced_model'] = results_reduced.final_loss
        test_result['chi2_test_statistics'] = chi2_stat
        test_result['p_value'] = 1 - stats.chi2.cdf(chi2_stat, df=lrt_row['lrt_df'])
        test_result['training_diverged_reduced_model'] = results_reduced.training_diverged
        test_result['training_converged_within_max_epochs_reduced_model'] = results_reduced.converged_within_max_epochs
        test_result['num_epochs_full_model'] = results_full.num_epochs
        test_result['num_epochs_reduced_model'] = results_reduced.num_epochs
        test_results_list.append(test_result)

    return model_param_df, pd.DataFrame(test_results_list), model_full.state_dict()


def make_global_pol2_closures(
        model: GlobalPol2Model,
        optimizer: optim.Optimizer,
        global_gene_data: GlobalGeneData,
        dataset_metadata: DatasetMetadata,
        reduced_matrix: Optional[torch.Tensor] = None,
) -> tuple[Callable[[], torch.Tensor], Callable[[], torch.Tensor]]:
    pol_2_total_loss = Pol2TotalLoss().to(next(model.parameters()).device)

    def _compute_loss() -> torch.Tensor:
        predicted_reads_exon, predicted_reads_intron, pi = model(
            design_matrix=dataset_metadata.design_matrix,
            log_library_sizes=dataset_metadata.log_library_sizes,
            isoform_length_offset=global_gene_data.isoform_length_offset,
            reduced_design_matrix=reduced_matrix,
        )
        return pol_2_total_loss(
            reads_exon=global_gene_data.exon_reads,
            reads_introns=global_gene_data.intron_reads,
            coverage=global_gene_data.coverage,
            predicted_reads_exon=predicted_reads_exon,
            predicted_reads_intron=predicted_reads_intron,
            pi=pi,
        )

    def closure() -> torch.Tensor:
        optimizer.zero_grad()
        loss = _compute_loss()
        if not torch.isfinite(loss):
            return torch.as_tensor(LOSS_CLAMP_VALUE, dtype=loss.dtype, device=loss.device)
        loss.backward()
        return loss

    def evaluate_loss() -> torch.Tensor:
        with torch.no_grad():
            return _compute_loss()

    return closure, evaluate_loss


def _train_two_stage(
        model: nn.Module,
        global_param_names: set[str],
        make_closures: Callable[[optim.Optimizer], tuple[Callable, Callable]],
        label: str = '',
        verbose: bool = False,
) -> tuple[nn.Module, TrainingResults]:
    prefix = f'[{label}] ' if label else ''

    # Stage 1: freeze per-gene/per-intron params; optimize only shared params
    print(f'{prefix}Stage 1: optimizing shared params ({", ".join(sorted(global_param_names))})', flush=True)
    for name, param in model.named_parameters():
        param.requires_grad_(name in global_param_names)
    optimizer_stage_1 = optim.LBFGS(
        [p for p in model.parameters() if p.requires_grad],
        lr=1.0, max_iter=20, tolerance_change=1e-9, tolerance_grad=1e-7,
        history_size=100, line_search_fn='strong_wolfe',
    )
    closure_stage_1, eval_stage_1 = make_closures(optimizer_stage_1)
    model, _ = train_model(model, optimizer_stage_1, closure_stage_1, eval_stage_1, max_epochs=500, verbose=verbose)

    # Stage 2: unfreeze all, refine jointly from warm start
    print(f'{prefix}Stage 2: joint optimization (all params)', flush=True)
    for param in model.parameters():
        param.requires_grad_(True)
    optimizer_stage_2 = make_lbfgs_optimizer(model)
    closure_stage_2, eval_stage_2 = make_closures(optimizer_stage_2)
    model, training_results = train_model(model, optimizer_stage_2, closure_stage_2, eval_stage_2, max_epochs=500, verbose=verbose)
    return model, training_results


def get_global_pol2_model_results(
        global_gene_data: GlobalGeneData,
        dataset_metadata: DatasetMetadata,
        device: str = 'cpu',
        verbose: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    global_gene_data = global_gene_data.to(device)
    dataset_metadata = dataset_metadata.to(device)

    model = GlobalPol2Model(
        feature_names=dataset_metadata.feature_names,
        gene_names=global_gene_data.gene_names,
        intron_names=global_gene_data.intron_names,
        gene_idx=global_gene_data.gene_idx,
    ).to(device)
    model.initialize_parameters(global_gene_data, dataset_metadata.library_sizes, dataset_metadata.design_matrix)

    model, training_results = _train_two_stage(
        model,
        global_param_names={'beta', 'gamma'},
        make_closures=lambda opt: make_global_pol2_closures(model, opt, global_gene_data, dataset_metadata),
        label='Global Pol2 | full model',
        verbose=verbose,
    )

    feature_names = dataset_metadata.feature_names
    param_rows: list[dict] = []
    for param_name in ('beta', 'gamma'):
        param_val = getattr(model, param_name).detach()
        for feature_index, feature_name in enumerate(feature_names):
            param_rows.append({
                'parameter_type': param_name,
                'feature_name': feature_name,
                'intron_name': None,
                'value': param_val[feature_index].item(),
                'loss_full_model': training_results.final_loss,
                'training_diverged_full_model': training_results.training_diverged,
                'training_converged_within_max_epochs_full_model': training_results.converged_within_max_epochs,
            })
    model_param_df = pd.DataFrame(param_rows)

    test_results_list: list[dict] = []
    for _, lrt_row in dataset_metadata.lrt_metadata.iterrows():
        reduced_matrix = dataset_metadata.reduced_matrices[lrt_row['test_id']].to(device)

        for tested_parameter in (TestableParameters.BETA, TestableParameters.GAMMA):
            lrt_spec = LRTSpecification(
                num_features_reduced_matrix=reduced_matrix.shape[1],
                tested_parameter=tested_parameter,
            )
            model_reduced = GlobalPol2Model(
                feature_names=feature_names,
                gene_names=global_gene_data.gene_names,
                intron_names=global_gene_data.intron_names,
                gene_idx=global_gene_data.gene_idx,
                lrt_specification=lrt_spec,
            ).to(device)

            full_state = model.state_dict()
            reduced_state = model_reduced.state_dict()
            for key in full_state:
                reduced_state[key] = full_state[key]

            if reduced_matrix.shape[1] > 0:
                full_lfc = model.beta if tested_parameter == TestableParameters.BETA else model.gamma
                full_contribution = dataset_metadata.design_matrix @ full_lfc.detach()
                reduced_state['reduced_lfc'] = torch.linalg.lstsq(reduced_matrix, full_contribution).solution

            model_reduced.load_state_dict(reduced_state)

            model_reduced, results_reduced = _train_two_stage(
                model_reduced,
                global_param_names={'beta', 'gamma', 'reduced_lfc'},
                make_closures=lambda opt: make_global_pol2_closures(
                    model_reduced, opt, global_gene_data, dataset_metadata, reduced_matrix,
                ),
                label=f'Global Pol2 | LRT {lrt_row["test_id"]} | {tested_parameter}',
                verbose=verbose,
            )

            full_lfc_vals = model.beta.detach() if tested_parameter == TestableParameters.BETA else model.gamma.detach()
            lfc_positive = (0.0 if pd.isna(lrt_row['lfc_column_positive'])
                            else full_lfc_vals[feature_names.index(lrt_row['lfc_column_positive'])].item())
            lfc_negative = (0.0 if pd.isna(lrt_row['lfc_column_negative'])
                            else full_lfc_vals[feature_names.index(lrt_row['lfc_column_negative'])].item())

            chi2_stat = 2 * (results_reduced.final_loss - training_results.final_loss)
            test_result = lrt_row.to_dict()
            test_result['tested_parameter'] = tested_parameter
            test_result['lfc'] = lfc_positive - lfc_negative
            test_result['loss_full_model'] = training_results.final_loss
            test_result['loss_reduced_model'] = results_reduced.final_loss
            test_result['chi2_test_statistics'] = chi2_stat
            test_result['p_value'] = 1 - stats.chi2.cdf(chi2_stat, df=lrt_row['lrt_df'])
            test_result['training_diverged_reduced_model'] = results_reduced.training_diverged
            test_result[
                'training_converged_within_max_epochs_reduced_model'] = results_reduced.converged_within_max_epochs
            test_result['num_epochs_full_model'] = training_results.num_epochs
            test_result['num_epochs_reduced_model'] = results_reduced.num_epochs
            test_results_list.append(test_result)

    return model_param_df, pd.DataFrame(test_results_list)


def get_regularized_splicing_model_results(
        input_data: GeneData | IntronData,
        dataset_metadata: DatasetMetadata,
        hot_start_state_dict: StateDict,
        regularization_coefficients_df: pd.DataFrame,
        device: str = 'cpu',
) -> pd.DataFrame:
    coverage = input_data.coverage.to(device)
    dataset_metadata = dataset_metadata.to(device)
    regularization_coefficients_df = regularization_coefficients_df.set_index(['parameter_type', 'feature_name'])

    model_regularized = SplicingModel(
        feature_names=dataset_metadata.feature_names,
        intron_names=input_data.intron_names,
    ).to(device)
    model_regularized.load_state_dict(hot_start_state_dict)

    regularization_coefficients_lfc = torch.zeros((len(dataset_metadata.feature_names), 1)).to(device)
    for feature_index, feature_name in enumerate(dataset_metadata.feature_names):
        regularization_coefficients_lfc[feature_index] = regularization_coefficients_df.loc[
            ('lfc', feature_name)]['lambda']

    optimizer_regularized = make_lbfgs_optimizer(model_regularized)
    closure_regularized, evaluate_loss_regularized = make_splicing_closures(
        model_regularized, optimizer_regularized, coverage, dataset_metadata.design_matrix,
        regularization_coefficients_lfc=regularization_coefficients_lfc,
    )
    model_regularized, training_results_regularized = train_model(
        model_regularized, optimizer_regularized, closure_regularized, evaluate_loss_regularized,
    )

    model_param_df = model_regularized.get_param_df()
    model_param_df['gene_name'] = input_data.gene_name
    model_param_df['loss_regularized_model'] = training_results_regularized.final_loss
    model_param_df['training_diverged_regularized_model'] = training_results_regularized.training_diverged
    model_param_df[
        'training_converged_within_max_epochs_regularized_model'] = training_results_regularized.converged_within_max_epochs
    return model_param_df

import numpy as np
import numpy as np
import pandas as pd
import statsmodels.api as sm
import torch
import torch.nn as nn
from dataclasses import dataclass
from enum import StrEnum
from typing import Optional, NamedTuple


def _init_alpha_poisson(X: np.ndarray, y: np.ndarray, offset: np.ndarray) -> np.ndarray:
    """
    Initialise alpha via Poisson GLM (log link, statsmodels IRLS).
    X: (S, F), y: (S,) or (S, G), offset: same shape as y.
    Returns alpha of shape (F,) or (G, F).
    Falls back to zeros for any gene where the GLM fails (e.g. all-zero counts).
    """
    if y.ndim == 1:
        try:
            return sm.GLM(y, X, family=sm.families.Poisson(), offset=offset).fit(disp=False).params
        except Exception:
            return np.zeros(X.shape[1])
    G = y.shape[1]
    alpha = np.zeros((G, X.shape[1]))
    for g in range(G):
        try:
            alpha[g] = sm.GLM(y[:, g], X, family=sm.families.Poisson(),
                               offset=offset[:, g]).fit(disp=False).params
        except Exception:
            pass
    return alpha


def safe_exp(x: torch.Tensor, output_threshold: float = 1e20) -> torch.Tensor:
    output_threshold_tensor = torch.tensor(output_threshold, dtype=x.dtype, device=x.device)
    input_threshold = torch.log(output_threshold_tensor)
    return torch.where(
        x <= input_threshold,
        torch.exp(x),
        output_threshold_tensor + output_threshold_tensor * (x - input_threshold)
    )


@dataclass
class IntronData:
    intron_name: str
    gene_name: str
    coverage: torch.Tensor  # shape (num_samples, 1, num_coverage_bins)

    @property
    def intron_names(self) -> list[str]:
        return [self.intron_name]

    def to(self, device):
        self.coverage = self.coverage.to(device)
        return self


@dataclass
class GeneData:
    gene_name: str
    intron_names: list[str]
    exon_reads: torch.Tensor  # shape (num_samples,)
    intron_reads: torch.Tensor  # shape (num_samples, num_introns)
    coverage: torch.Tensor  # shape (num_samples, num_introns, num_coverage_bins)
    isoform_length_offset: torch.Tensor  # shape (num_samples,)

    def to(self, device):
        self.exon_reads = self.exon_reads.to(device)
        self.intron_reads = self.intron_reads.to(device)
        self.coverage = self.coverage.to(device)
        self.isoform_length_offset = self.isoform_length_offset.to(device)
        return self


@dataclass
class DatasetMetadata:
    design_matrix: torch.Tensor
    library_sizes: torch.Tensor
    log_library_sizes: torch.Tensor
    feature_names: list[str]
    sample_names: list[str]
    lrt_metadata: pd.DataFrame
    reduced_matrices: dict[str, torch.Tensor]
    num_coverage_bins: int = 100

    def to(self, device):
        self.design_matrix = self.design_matrix.to(device)
        self.library_sizes = self.library_sizes.to(device)
        self.log_library_sizes = self.log_library_sizes.to(device)
        for name, matrix in self.reduced_matrices.items():
            self.reduced_matrices[name] = matrix.to(device)
        return self


class TestableParameters(StrEnum):
    ALPHA = 'alpha'
    BETA = 'beta'
    GAMMA = 'gamma'


class LRTSpecification(NamedTuple):
    num_features_reduced_matrix: int
    tested_parameter: TestableParameters
    tested_intron: Optional[str] = None


class Pol2Model(nn.Module):

    def __init__(self,
                 feature_names: list[str],
                 intron_names: list[str],
                 intron_specific_lfc: bool,
                 lrt_specification: Optional[LRTSpecification] = None
                 ):
        super().__init__()
        self.feature_names = feature_names
        self.intron_names = intron_names

        num_features = len(feature_names)
        num_introns = len(intron_names)

        self.alpha = nn.Parameter(torch.zeros(num_features))
        self.intercept_exon = nn.Parameter(torch.zeros(1))

        if intron_specific_lfc:
            self.beta = nn.Parameter(torch.zeros(num_features, num_introns))
            self.gamma = nn.Parameter(torch.zeros(num_features, num_introns))
        else:
            self.beta = nn.Parameter(torch.zeros(num_features, 1))
            self.gamma = nn.Parameter(torch.zeros(num_features, 1))
        self.intron_specific_lfc = intron_specific_lfc

        self.intercept_intron = nn.Parameter(torch.zeros(num_introns))
        self.theta = nn.Parameter(torch.zeros(num_introns))

        self.lrt_specification = lrt_specification
        self.tested_intron_index: Optional[int] = None
        if self.lrt_specification is not None:
            self.reduced_lfc = nn.Parameter(torch.zeros(self.lrt_specification.num_features_reduced_matrix))
            self.tested_intron_index: Optional[
                int] = None if not intron_specific_lfc or self.lrt_specification.tested_parameter == TestableParameters.ALPHA else self.intron_names.index(
                self.lrt_specification.tested_intron)

    def initialize_parameters(self,
                              gene_data: GeneData,
                              library_sizes: torch.Tensor,
                              design_matrix: torch.Tensor,
                              pi_eps: float = 0.05,
                              num_pi_grid_points: int = 20) -> None:
        with torch.no_grad():
            intercept_exon_scalar = torch.log(gene_data.exon_reads.mean() / library_sizes.mean())
            self.intercept_exon.data[:] = intercept_exon_scalar  # preserve shape [1]

            offset = (intercept_exon_scalar + torch.log(library_sizes) + gene_data.isoform_length_offset).cpu().numpy()
            alpha_init = _init_alpha_poisson(design_matrix.cpu().numpy(), gene_data.exon_reads.cpu().numpy(), offset)
            self.alpha.data.copy_(torch.from_numpy(alpha_init).to(self.alpha.dtype))

            intercept_intron_vector = torch.log(gene_data.intron_reads.mean(dim=0) / library_sizes.mean() / 2)
            self.intercept_intron.data.copy_(intercept_intron_vector)

            aggregated_coverage = gene_data.coverage.sum(dim=0)
            pi_grid = torch.linspace(
                pi_eps,
                1 - pi_eps,
                num_pi_grid_points,
                device=aggregated_coverage.device,
                dtype=aggregated_coverage.dtype,
            )
            coverage_loss = CoverageLoss(num_position_coverage=aggregated_coverage.shape[1])
            coverage_loss_grid = coverage_loss.loss_for_pi_grid(pi_grid, aggregated_coverage)
            best_pi = pi_grid[coverage_loss_grid.argmin(dim=0)]
            self.theta.data.copy_(torch.logit(best_pi, eps=pi_eps))

    def forward(self,
                design_matrix: torch.Tensor,
                log_library_sizes: torch.Tensor,
                isoform_length_offset: torch.Tensor,
                reduced_design_matrix: Optional[torch.Tensor] = None):

        if self.lrt_specification is None:
            gene_expression_term = design_matrix @ self.alpha
            speed_term = design_matrix @ self.beta
            splicing_term = design_matrix @ self.gamma
        else:
            if reduced_design_matrix is None:
                raise ValueError(
                    "reduced_design_matrix must be provided in the LRT mode (i.e. when lrt_specification is not None.")
            gene_expression_term = reduced_design_matrix @ self.reduced_lfc if self.lrt_specification.tested_parameter == TestableParameters.ALPHA else design_matrix @ self.alpha

            if self.intron_specific_lfc:
                speed_term = design_matrix @ self.beta
                splicing_term = design_matrix @ self.gamma
                if self.lrt_specification.tested_parameter == TestableParameters.BETA:
                    speed_term[:, self.tested_intron_index] = reduced_design_matrix @ self.reduced_lfc
                elif self.lrt_specification.tested_parameter == TestableParameters.GAMMA:
                    splicing_term[:, self.tested_intron_index] = reduced_design_matrix @ self.reduced_lfc

            else:
                speed_term = (reduced_design_matrix @ self.reduced_lfc).unsqueeze(
                    1) if self.lrt_specification.tested_parameter == TestableParameters.BETA else design_matrix @ self.beta
                splicing_term = (reduced_design_matrix @ self.reduced_lfc).unsqueeze(
                    1) if self.lrt_specification.tested_parameter == TestableParameters.GAMMA else design_matrix @ self.gamma

        predicted_log_reads_exon = self.intercept_exon + log_library_sizes + isoform_length_offset + gene_expression_term

        pi = torch.sigmoid(self.theta - speed_term - splicing_term)

        intron_gene_expression_term = self.intercept_intron + log_library_sizes.unsqueeze(
            1) + gene_expression_term.unsqueeze(1)
        reads_intronic_polymerases = safe_exp(intron_gene_expression_term + self.theta - speed_term)
        reads_unspliced_transcripts = safe_exp(intron_gene_expression_term + splicing_term)
        predicted_reads_intron = reads_intronic_polymerases + reads_unspliced_transcripts

        return safe_exp(predicted_log_reads_exon), predicted_reads_intron, pi

    def get_param_df(self) -> pd.DataFrame:
        model_parameters = dict(self.named_parameters())
        parameter_data: list[dict] = []
        for param_name, param_value in model_parameters.items():
            if param_name == 'intercept_exon':
                parameter_data.append({'parameter_type': param_name,
                                       'intron_name': None,
                                       'feature_name': None,
                                       'value': param_value.item()})
            elif param_name == 'alpha':
                for feature_index, feature_name in enumerate(self.feature_names):
                    parameter_data.append({'parameter_type': param_name,
                                           'intron_name': None,
                                           'feature_name': feature_name,
                                           'value': param_value[feature_index].item()})
            elif param_name in ('beta', 'gamma'):
                for feature_index, feature_name in enumerate(self.feature_names):
                    if self.intron_specific_lfc:
                        for intron_index, intron_name in enumerate(self.intron_names):
                            parameter_data.append({'parameter_type': param_name,
                                                   'intron_name': intron_name,
                                                   'feature_name': feature_name,
                                                   'value': param_value[feature_index, intron_index].item()})
                    else:
                        parameter_data.append({'parameter_type': param_name,
                                               'intron_name': None,
                                               'feature_name': feature_name,
                                               'value': param_value[feature_index, 0].item()})

            elif param_name in ('intercept_intron', 'theta'):
                for intron_index, intron_name in enumerate(self.intron_names):
                    parameter_data.append({'parameter_type': param_name,
                                           'intron_name': intron_name,
                                           'feature_name': None,
                                           'value': param_value[intron_index].item()})
            else:
                raise RuntimeError(f"Unexpected parameter name: {param_name}")
        df_param = pd.DataFrame(data=parameter_data)
        return df_param


class SplicingModel(nn.Module):

    def __init__(self,
                 feature_names: list[str],
                 intron_names: list[str]):
        super().__init__()
        self.feature_names = feature_names
        self.intron_names = intron_names

        num_features = len(feature_names)
        num_introns = len(intron_names)

        self.lfc = nn.Parameter(torch.zeros(num_features, 1))
        self.theta = nn.Parameter(torch.zeros(num_introns))

    def initialize_theta(self,
                         coverage: torch.Tensor,
                         pi_eps: float = 0.05,
                         num_pi_grid_points: int = 20) -> None:
        aggregated_coverage = coverage.sum(dim=0)
        pi_grid = torch.linspace(
            pi_eps,
            1 - pi_eps,
            num_pi_grid_points,
            device=aggregated_coverage.device,
            dtype=aggregated_coverage.dtype,
        )
        coverage_loss = CoverageLoss(num_position_coverage=aggregated_coverage.shape[1])
        coverage_loss_grid = coverage_loss.loss_for_pi_grid(pi_grid, aggregated_coverage)
        best_pi = pi_grid[coverage_loss_grid.argmin(dim=0)]
        self.theta.data.copy_(torch.logit(best_pi, eps=pi_eps))

    def forward(self, design_matrix: torch.Tensor):
        logit_term = design_matrix @ self.lfc
        pi = torch.sigmoid(self.theta + logit_term)
        return pi

    def get_param_df(self) -> pd.DataFrame:
        model_parameters = dict(self.named_parameters())
        parameter_data: list[dict] = []
        for param_name, param_value in model_parameters.items():
            if param_name == 'lfc':
                for feature_index, feature_name in enumerate(self.feature_names):
                    parameter_data.append({'parameter_type': param_name,
                                           'intron_name': None,
                                           'feature_name': feature_name,
                                           'value': param_value[feature_index, 0].item()})
            elif param_name == 'theta':
                for intron_index, intron_name in enumerate(self.intron_names):
                    parameter_data.append({'parameter_type': param_name,
                                           'intron_name': intron_name,
                                           'feature_name': None,
                                           'value': param_value[intron_index].item()})
            else:
                raise RuntimeError(f"Unexpected parameter name: {param_name}")
        df_param = pd.DataFrame(data=parameter_data)
        return df_param


@dataclass
class GlobalGeneData:
    gene_names: list[str]
    intron_names: list[str]
    exon_reads: torch.Tensor          # (S, G)
    intron_reads: torch.Tensor        # (S, N)
    coverage: torch.Tensor            # (S, N, B)
    isoform_length_offset: torch.Tensor  # (S, G)
    gene_idx: torch.Tensor            # (N,) int64 — maps intron index → gene index

    def to(self, device):
        self.exon_reads = self.exon_reads.to(device)
        self.intron_reads = self.intron_reads.to(device)
        self.coverage = self.coverage.to(device)
        self.isoform_length_offset = self.isoform_length_offset.to(device)
        self.gene_idx = self.gene_idx.to(device)
        return self


def concat_gene_data_list(gene_data_list: list[GeneData]) -> GlobalGeneData:
    gene_idx = torch.cat([
        torch.full((len(g.intron_names),), i, dtype=torch.long)
        for i, g in enumerate(gene_data_list)
    ])
    return GlobalGeneData(
        gene_names=[g.gene_name for g in gene_data_list],
        intron_names=[name for g in gene_data_list for name in g.intron_names],
        exon_reads=torch.stack([g.exon_reads for g in gene_data_list], dim=1),
        intron_reads=torch.cat([g.intron_reads for g in gene_data_list], dim=1),
        coverage=torch.cat([g.coverage for g in gene_data_list], dim=1),
        isoform_length_offset=torch.stack([g.isoform_length_offset for g in gene_data_list], dim=1),
        gene_idx=gene_idx,
    )


class GlobalPol2Model(nn.Module):

    def __init__(self,
                 feature_names: list[str],
                 gene_names: list[str],
                 intron_names: list[str],
                 gene_idx: torch.Tensor,
                 lrt_specification: Optional[LRTSpecification] = None,
                 ):
        super().__init__()
        self.feature_names = feature_names
        self.gene_names = gene_names
        self.intron_names = intron_names

        num_features = len(feature_names)
        num_genes = len(gene_names)
        num_introns = len(intron_names)

        self.alpha = nn.Parameter(torch.zeros(num_genes, num_features))
        self.beta = nn.Parameter(torch.zeros(num_features))
        self.gamma = nn.Parameter(torch.zeros(num_features))
        self.intercept_exon = nn.Parameter(torch.zeros(num_genes))
        self.intercept_intron = nn.Parameter(torch.zeros(num_introns))
        self.theta = nn.Parameter(torch.zeros(num_introns))

        self.register_buffer('gene_idx', gene_idx)

        self.lrt_specification = lrt_specification
        if lrt_specification is not None:
            if lrt_specification.tested_parameter not in (TestableParameters.BETA, TestableParameters.GAMMA):
                raise ValueError(
                    f"GlobalPol2Model LRT only supports BETA and GAMMA, got {lrt_specification.tested_parameter!r}.")
            self.reduced_lfc = nn.Parameter(torch.zeros(lrt_specification.num_features_reduced_matrix))

    def initialize_parameters(self,
                              global_gene_data: GlobalGeneData,
                              library_sizes: torch.Tensor,
                              design_matrix: torch.Tensor,
                              pi_eps: float = 0.05,
                              num_pi_grid_points: int = 20) -> None:
        with torch.no_grad():
            self.intercept_exon.data.copy_(
                torch.log(global_gene_data.exon_reads.mean(dim=0) / library_sizes.mean())
            )

            log_lib = torch.log(library_sizes)
            offset = (self.intercept_exon.unsqueeze(0) + log_lib.unsqueeze(1)
                      + global_gene_data.isoform_length_offset).cpu().numpy()  # (S, G)
            alpha_init = _init_alpha_poisson(
                design_matrix.cpu().numpy(), global_gene_data.exon_reads.cpu().numpy(), offset
            )  # (G, F)
            self.alpha.data.copy_(torch.from_numpy(alpha_init).to(self.alpha.dtype))

            self.intercept_intron.data.copy_(
                torch.log(global_gene_data.intron_reads.mean(dim=0) / library_sizes.mean() / 2)
            )
            aggregated_coverage = global_gene_data.coverage.sum(dim=0)  # (N, B)
            pi_grid = torch.linspace(pi_eps, 1 - pi_eps, num_pi_grid_points,
                                     device=aggregated_coverage.device,
                                     dtype=aggregated_coverage.dtype)
            coverage_loss = CoverageLoss(num_position_coverage=aggregated_coverage.shape[1])
            coverage_loss_grid = coverage_loss.loss_for_pi_grid(pi_grid, aggregated_coverage)
            best_pi = pi_grid[coverage_loss_grid.argmin(dim=0)]
            self.theta.data.copy_(torch.logit(best_pi, eps=pi_eps))

    def forward(self,
                design_matrix: torch.Tensor,
                log_library_sizes: torch.Tensor,
                isoform_length_offset: torch.Tensor,
                reduced_design_matrix: Optional[torch.Tensor] = None):

        gene_expression_term = design_matrix @ self.alpha.T   # (S, G)

        predicted_log_reads_exon = (
            self.intercept_exon
            + log_library_sizes.unsqueeze(1)
            + isoform_length_offset
            + gene_expression_term
        )

        if self.lrt_specification is None:
            speed_term = (design_matrix @ self.beta).unsqueeze(1)
            splicing_term = (design_matrix @ self.gamma).unsqueeze(1)
        else:
            if reduced_design_matrix is None:
                raise ValueError("reduced_design_matrix must be provided in LRT mode.")
            reduced_term = (reduced_design_matrix @ self.reduced_lfc).unsqueeze(1)
            speed_term = (
                reduced_term if self.lrt_specification.tested_parameter == TestableParameters.BETA
                else (design_matrix @ self.beta).unsqueeze(1)
            )
            splicing_term = (
                reduced_term if self.lrt_specification.tested_parameter == TestableParameters.GAMMA
                else (design_matrix @ self.gamma).unsqueeze(1)
            )

        gene_expr_per_intron = gene_expression_term[:, self.gene_idx]  # (S, N)

        intron_gene_expression_term = (
            self.intercept_intron
            + log_library_sizes.unsqueeze(1)
            + gene_expr_per_intron
        )

        pi = torch.sigmoid(self.theta - speed_term - splicing_term)
        reads_intronic_polymerases = safe_exp(intron_gene_expression_term + self.theta - speed_term)
        reads_unspliced_transcripts = safe_exp(intron_gene_expression_term + splicing_term)
        predicted_reads_intron = reads_intronic_polymerases + reads_unspliced_transcripts

        return safe_exp(predicted_log_reads_exon), predicted_reads_intron, pi


class CoverageLoss(nn.Module):

    def __init__(self, num_position_coverage: int = 100):
        super().__init__()
        locations = torch.linspace(start=1 / (2 * num_position_coverage),
                                   end=1 - 1 / (2 * num_position_coverage),
                                   steps=num_position_coverage)
        location_term = 1 - 2 * locations
        self.register_buffer("location_term", location_term)

    def forward(self, pi, coverage):
        loss_per_location = -torch.log(1 + pi.unsqueeze(2) * self.location_term)
        return torch.sum(loss_per_location * coverage)

    def loss_for_pi_grid(self, candidate_pi: torch.Tensor, aggregated_coverage: torch.Tensor) -> torch.Tensor:
        """
        Compute coverage loss for a grid of candidate pi values.
        Args:
            candidate_pi: Tensor of shape (k,), candidate values of pi.
            aggregated_coverage: Tensor of shape (i, b), coverage summed over samples.
        Returns:
            Tensor of shape (k, i) containing the total loss for each candidate pi
            and each intron.
        """
        loss_terms = -torch.log(
            1 + candidate_pi[:, None] * self.location_term[None, :]
        )  # shape (k, b)
        return torch.einsum("kb,ib->ki", loss_terms, aggregated_coverage)


class Pol2TotalLoss(nn.Module):
    def __init__(self, num_position_coverage: int = 100):
        super().__init__()
        self.loss_function_exon = nn.PoissonNLLLoss(log_input=False, full=True, reduction='sum')
        self.loss_function_intron = nn.PoissonNLLLoss(log_input=False, full=True, reduction='sum')
        self.loss_function_coverage = CoverageLoss(num_position_coverage=num_position_coverage)

    def forward(self,
                reads_exon: torch.Tensor,
                reads_introns: torch.Tensor,
                coverage: torch.Tensor,
                predicted_reads_exon: torch.Tensor,
                predicted_reads_intron: torch.Tensor,
                pi: torch.Tensor):
        loss_exon = self.loss_function_exon(predicted_reads_exon,
                                            reads_exon)
        loss_intron = self.loss_function_intron(predicted_reads_intron, reads_introns)
        loss_coverage = self.loss_function_coverage(pi, coverage)

        total_loss = loss_exon + loss_intron + loss_coverage
        return total_loss

import pandas as pd
import torch
import torch.nn as nn
from dataclasses import dataclass
from enum import StrEnum
from typing import Optional, NamedTuple


def safe_exp(x: torch.Tensor, output_threshold: float = 1e20) -> torch.Tensor:
    output_threshold_tensor = torch.tensor(output_threshold, dtype=x.dtype, device=x.device)
    input_threshold = torch.log(output_threshold_tensor)
    return torch.where(
        x <= input_threshold,
        torch.exp(x),
        output_threshold_tensor + output_threshold_tensor * (x - input_threshold)
    )


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

    def initialize_intercepts(self, gene_data: GeneData, library_sizes: torch.Tensor) -> None:
        with torch.no_grad():
            intercept_exon_scalar = torch.log(gene_data.exon_reads.mean() / library_sizes.mean())
            self.intercept_exon.data[:] = intercept_exon_scalar  # preserve shape [1]

            intercept_intron_vector = torch.log(gene_data.intron_reads.mean(dim=0) / library_sizes.mean() / 2)
            self.intercept_intron.data.copy_(intercept_intron_vector)

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

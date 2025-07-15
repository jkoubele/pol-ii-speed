from dataclasses import dataclass
from typing import Optional, NamedTuple

import torch
import torch.nn as nn


@dataclass
class GeneData:
    exon_reads: torch.Tensor
    intron_names: list[str]
    intron_reads: dict[str, torch.Tensor]
    coverage_density: dict[str, torch.Tensor]

    def to(self, device):
        self.exon_reads = self.exon_reads.to(device)
        self.intron_reads = {key: value.to(device) for key, value in self.intron_reads.items()}
        self.coverage_density = {key: value.to(device) for key, value in self.coverage_density.items()}
        return self


@dataclass
class GeneDataWithSolution:
    gene_data: GeneData
    mu_1: float
    mu_2: float
    nu_1: dict[str, float]
    nu_2: dict[str, float]
    phi_1: dict[str, float]
    phi_2: dict[str, float]


class ParameterMask(NamedTuple):
    feature_index: int
    value: float
    intron_name: Optional[str] = None


class Pol2ModelOptimized(nn.Module):

    def __init__(self,
                 feature_names: list[str],
                 intron_names: list[str],
                 ):
        super().__init__()
        self.feature_names = feature_names
        self.intron_names = intron_names

        num_features = len(feature_names)
        num_introns = len(intron_names)

        self.alpha = nn.Parameter(torch.zeros(num_features))
        self.intercept_exon = nn.Parameter(torch.zeros(1))

        self.beta = nn.Parameter(torch.zeros(num_features, num_introns))
        self.gamma = nn.Parameter(torch.zeros(num_features, num_introns))

        self.intercept_intron = nn.Parameter(torch.zeros(num_introns))
        self.log_phi_zero = nn.Parameter(torch.zeros(num_introns))

        self.mask_alpha = None
        self.mask_beta = None
        self.mask_gamma = None

    def initialize_intercepts(self, gene_data: GeneData, library_sizes: torch.Tensor) -> None:
        with torch.no_grad():
            device = self.intercept_exon.device  # use model’s parameter device
            exon_reads = gene_data.exon_reads.to(device)
            library_sizes = library_sizes.to(device)
            self.intercept_exon.data = torch.log(exon_reads.mean() / library_sizes.mean())
            for i, intron_name in enumerate(self.intron_names):
                intron_reads = gene_data.intron_reads[intron_name].to(device)
                mean_intron_value = intron_reads.mean() / library_sizes.mean()
                self.intercept_intron.data[i] = torch.log(mean_intron_value)

    def forward(self, X: torch.Tensor, log_library_sizes: torch.Tensor):
        alpha = self.alpha
        if self.mask_alpha is not None:
            pass
            # alpha = self.alpha.clone()
            # alpha[self.mask_alpha.feature_index] = self.mask_alpha.value

        gene_expression_term = X @ alpha
        predicted_log_reads_exon = self.intercept_exon + log_library_sizes + gene_expression_term

        beta = self.beta
        gamma = self.gamma
        # TODO: handle masking for LRT

        speed_term = X @ beta
        splicing_term = X @ gamma

        phi = torch.sigmoid(self.log_phi_zero - speed_term - splicing_term)  # (num_samples, num_introns)

        intron_gene_expression_term = self.intercept_intron + log_library_sizes.unsqueeze(
            1) + gene_expression_term.unsqueeze(1)
        reads_intronic_polymerases = torch.exp(intron_gene_expression_term + self.log_phi_zero - speed_term)
        reads_unspliced_transcripts = torch.exp(intron_gene_expression_term + splicing_term)
        predicted_reads_intron = reads_intronic_polymerases + reads_unspliced_transcripts

        return predicted_log_reads_exon, predicted_reads_intron, phi


class Pol2Model(nn.Module):

    def __init__(self,
                 num_features: int,
                 intron_names: list[str],
                 exon_intercept_init: Optional[float] = None,
                 intron_intercepts_init: Optional[dict[str, float]] = None,
                 mask_alpha: Optional[ParameterMask] = None,
                 mask_beta: Optional[ParameterMask] = None,
                 mask_gamma: Optional[ParameterMask] = None):
        super().__init__()
        self.intron_names = intron_names

        self.alpha = nn.Parameter(torch.zeros(num_features))
        self.intercept_exon = nn.Parameter(
            torch.tensor(exon_intercept_init if exon_intercept_init is not None else 0.0))

        self.beta = nn.ParameterDict({
            intron: nn.Parameter(torch.zeros(num_features))
            for intron in intron_names})
        self.gamma = nn.ParameterDict({
            intron: nn.Parameter(torch.zeros(num_features))
            for intron in intron_names})

        self.intercept_intron = nn.ParameterDict({
            intron_name: nn.Parameter(
                torch.tensor(intron_intercepts_init[intron_name] if intron_intercepts_init is not None else 0.0))
            for intron_name in intron_names})
        self.log_phi_zero = nn.ParameterDict({
            intron: nn.Parameter(torch.tensor(0.0))
            for intron in intron_names})

        self.mask_alpha = mask_alpha
        self.mask_beta = mask_beta
        self.mask_gamma = mask_gamma

    def forward(self, X: torch.Tensor, log_library_sizes: torch.Tensor):
        alpha = self.alpha
        if self.mask_alpha is not None:
            alpha = self.alpha.clone()
            alpha[self.mask_alpha.feature_index] = self.mask_alpha.value
        gene_expression_term = X @ alpha
        predicted_log_reads_exon = self.intercept_exon + log_library_sizes + gene_expression_term

        predicted_reads_intron = {}
        phi = {}
        for intron_name in self.intron_names:
            beta = self.beta[intron_name]
            if self.mask_beta is not None and self.mask_beta.intron_name == intron_name:
                beta = self.beta[intron_name].clone()
                beta[self.mask_beta.feature_index] = self.mask_beta.value

            gamma = self.gamma[intron_name]
            if self.mask_gamma is not None and self.mask_gamma.intron_name == intron_name:
                gamma = self.gamma[intron_name].clone()
                gamma[self.mask_gamma.feature_index] = self.mask_gamma.value

            speed_term = X @ beta
            splicing_term = X @ gamma

            phi[intron_name] = torch.sigmoid(self.log_phi_zero[intron_name] - speed_term - splicing_term)

            predicted_reads_intron[intron_name] = torch.exp(
                self.intercept_intron[intron_name] + self.log_phi_zero[
                    intron_name] + log_library_sizes + gene_expression_term - speed_term) + torch.exp(
                self.intercept_intron[intron_name] + log_library_sizes + gene_expression_term + splicing_term)

        return predicted_log_reads_exon, predicted_reads_intron, phi


class CoverageLoss(nn.Module):

    def __init__(self, num_position_coverage: int = 100):
        super().__init__()
        locations = torch.linspace(start=1 / (2 * num_position_coverage),
                                   end=1 - 1 / (2 * num_position_coverage),
                                   steps=num_position_coverage).unsqueeze(0)
        self.register_buffer("locations", locations)

    def forward(self, phi, num_intronic_reads, coverage):
        loss_per_location = -torch.log(1 + phi.unsqueeze(1) - 2 * phi.unsqueeze(1) * self.locations)
        return torch.sum(loss_per_location * coverage * num_intronic_reads.unsqueeze(1))


class Pol2TotalLoss(nn.Module):
    def __init__(self, num_position_coverage: int = 100):
        super().__init__()
        self.loss_function_exon = nn.PoissonNLLLoss(log_input=True, full=True, reduction='sum')
        self.loss_function_intron = nn.PoissonNLLLoss(log_input=False, full=True, reduction='sum')
        self.loss_function_coverage = CoverageLoss(num_position_coverage=num_position_coverage)

    def forward(self, gene_data: GeneData,
                predicted_log_reads_exon: torch.Tensor,
                predicted_reads_intron: dict[str, torch.Tensor],
                phi: dict[str, torch.Tensor]):
        loss_exon = self.loss_function_exon(predicted_log_reads_exon, gene_data.exon_reads)

        loss_intron = sum(self.loss_function_intron(predicted_reads_intron[name], gene_data.intron_reads[name])
                          for name in gene_data.intron_names)

        loss_coverage = sum(self.loss_function_coverage(phi[name],
                                                        gene_data.intron_reads[name],
                                                        gene_data.coverage_density[name])
                            for name in gene_data.intron_names)

        total_loss = loss_exon + loss_intron + loss_coverage
        return total_loss


class CoverageLossOptimized(nn.Module):

    def __init__(self, num_position_coverage: int = 100):
        super().__init__()
        locations = torch.linspace(start=1 / (2 * num_position_coverage),
                                   end=1 - 1 / (2 * num_position_coverage),
                                   steps=num_position_coverage)
        self.register_buffer("locations", locations)

    def forward(self, phi, coverage):
        loss_per_location = -torch.log(1 + phi.unsqueeze(2) * (1 - 2 * self.locations))
        return torch.sum(loss_per_location * coverage)


class Pol2TotalLossOptimized(nn.Module):
    def __init__(self, num_position_coverage: int = 100):
        super().__init__()
        self.loss_function_exon = nn.PoissonNLLLoss(log_input=True, full=True, reduction='sum')
        self.loss_function_intron = nn.PoissonNLLLoss(log_input=False, full=True, reduction='sum')
        self.loss_function_coverage = CoverageLossOptimized(num_position_coverage=num_position_coverage)

    def forward(self,
                reads_exon: torch.Tensor,
                reads_introns: torch.Tensor,
                coverage: torch.Tensor,
                predicted_log_reads_exon: torch.Tensor,
                predicted_reads_intron: torch.Tensor,
                phi: torch.Tensor):
        loss_exon = self.loss_function_exon(predicted_log_reads_exon, reads_exon)

        loss_intron = self.loss_function_intron(predicted_reads_intron, reads_introns)

        loss_coverage = self.loss_function_coverage(phi, coverage)

        total_loss = loss_exon + loss_intron + loss_coverage
        return total_loss

from dataclasses import dataclass
from typing import Optional

import torch
import torch.nn as nn


@dataclass
class GeneData:
    exon_reads: torch.tensor
    intron_names: list[str]
    intron_reads: dict[str, torch.tensor]
    coverage_density: dict[str, torch.tensor]


@dataclass
class GeneDataWithSolution:
    gene_data: GeneData
    mu_1: float
    mu_2: float
    nu_1: dict[str, float]
    nu_2: dict[str, float]
    phi_1: dict[str, float]
    phi_2: dict[str, float]


class Pol2Model(nn.Module):

    def __init__(self,
                 num_features: int,
                 intron_names: list[str],
                 exon_intercept_init: Optional[float] = None,
                 intron_intercepts_init: Optional[dict[str, float]] = None):
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

    def forward(self, X, log_library_sizes):
        gene_expression_term = X @ self.alpha
        predicted_log_reads_exon = self.intercept_exon + log_library_sizes + gene_expression_term

        predicted_reads_intron = {}
        phi = {}
        for intron_name in self.intron_names:
            speed_term = X @ self.beta[intron_name]
            splicing_term = X @ self.gamma[intron_name]

            phi[intron_name] = torch.sigmoid(self.log_phi_zero[intron_name] - speed_term - splicing_term)

            predicted_reads_intron[intron_name] = torch.exp(
                self.intercept_intron[intron_name] + self.log_phi_zero[
                    intron_name] + log_library_sizes + gene_expression_term - speed_term) + torch.exp(
                self.intercept_intron[intron_name] + log_library_sizes + gene_expression_term + splicing_term)

        return predicted_log_reads_exon, predicted_reads_intron, phi


class CoverageLoss(nn.Module):

    def __init__(self, num_position_coverage: int = 100):
        super().__init__()
        self.locations = torch.linspace(start=1 / (2 * num_position_coverage),
                                        end=1 - 1 / (2 * num_position_coverage),
                                        steps=num_position_coverage).unsqueeze(0)

    def forward(self, phi, num_intronic_reads, coverage):
        loss_per_location = -torch.log(1 + phi.unsqueeze(1) - 2 * phi.unsqueeze(1) * self.locations)
        return torch.mean(loss_per_location * coverage * num_intronic_reads.unsqueeze(1))

from dataclasses import dataclass
from typing import Optional
import math

import pandas as pd
import torch
import torch.nn as nn

EXP_INPUT_CLAMP = 40


def safe_exp(x: torch.Tensor, output_threshold: float = 1e9) -> torch.Tensor:
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

    def to(self, device):
        self.exon_reads = self.exon_reads.to(device)
        self.intron_reads = self.intron_reads.to(device)
        self.coverage = self.coverage.to(device)
        return self


@dataclass
class DatasetMetadata:
    design_matrix: torch.Tensor
    library_sizes: torch.Tensor
    log_library_sizes: torch.Tensor
    feature_names: list[str]
    num_coverage_bins: int = 100

    def to(self, device):
        self.design_matrix = self.design_matrix.to(device)
        self.library_sizes = self.library_sizes.to(device)
        self.log_library_sizes = self.log_library_sizes.to(device)
        return self


@dataclass
class ParameterMask:
    logical_mask: torch.Tensor
    fixed_value_mask: Optional[torch.Tensor] = None

    def apply(self, param: torch.Tensor) -> torch.Tensor:
        param_masked = param * self.logical_mask
        if self.fixed_value_mask is not None:
            param_masked = param_masked + self.fixed_value_mask
        return param_masked


class Pol2Model(nn.Module):

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

        self.mask_alpha: Optional[ParameterMask] = None
        self.mask_beta: Optional[ParameterMask] = None
        self.mask_gamma: Optional[ParameterMask] = None

    def initialize_intercepts(self, gene_data: GeneData, library_sizes: torch.Tensor) -> None:
        with torch.no_grad():
            intercept_exon_scalar = torch.log(gene_data.exon_reads.mean() / library_sizes.mean())
            self.intercept_exon.data[:] = intercept_exon_scalar  # preserve shape [1]

            intercept_intron_vector = torch.log(gene_data.intron_reads.mean(axis=0) / library_sizes.mean())
            self.intercept_intron.data.copy_(intercept_intron_vector)

    def set_parameter_mask(self,
                           param_name: str,
                           feature_name: str,
                           intron_name: Optional[str] = None,
                           value=0.0) -> None:
        if param_name == 'alpha':
            logical_mask = torch.ones_like(self.alpha)
            feature_index = self.feature_names.index(feature_name)
            logical_mask[feature_index] = 0.0
            self.mask_alpha = ParameterMask(logical_mask=logical_mask)
            if value != 0:
                self.mask_alpha.fixed_value_mask = torch.zeros_like(self.alpha)
                self.mask_alpha.fixed_value_mask[feature_index] = value
        elif param_name in ('beta', 'gamma'):
            logical_mask = torch.ones_like(self.beta)  # beta and gamma have the same shape
            feature_index = self.feature_names.index(feature_name)
            intron_index = self.intron_names.index(intron_name)
            logical_mask[feature_index, intron_index] = 0.0
            fixed_value_mask = None
            if value != 0:
                fixed_value_mask = torch.zeros_like(self.beta)
                fixed_value_mask[feature_index, intron_index] = value
            parameter_mask = ParameterMask(logical_mask=logical_mask,
                                           fixed_value_mask=fixed_value_mask)
            if param_name == 'beta':
                self.mask_beta = parameter_mask
            elif param_name == 'gamma':
                self.mask_gamma = parameter_mask

        else:
            raise RuntimeError(f"Unexpected parameter name: {param_name}")

    def forward(self, design_matrix: torch.Tensor, log_library_sizes: torch.Tensor):
        alpha = self.alpha if self.mask_alpha is None else self.mask_alpha.apply(self.alpha)
        beta = self.beta if self.mask_beta is None else self.mask_beta.apply(self.beta)
        gamma = self.gamma if self.mask_gamma is None else self.mask_gamma.apply(self.gamma)

        gene_expression_term = design_matrix @ alpha
        predicted_log_reads_exon = (self.intercept_exon + log_library_sizes + gene_expression_term)

        speed_term = design_matrix @ beta
        splicing_term = design_matrix @ gamma

        phi = torch.sigmoid(self.log_phi_zero - speed_term - splicing_term)

        intron_gene_expression_term = self.intercept_intron + log_library_sizes.unsqueeze(
            1) + gene_expression_term.unsqueeze(1)
        reads_intronic_polymerases = safe_exp(intron_gene_expression_term + self.log_phi_zero - speed_term)
        reads_unspliced_transcripts = safe_exp(intron_gene_expression_term + splicing_term)
        predicted_reads_intron = reads_intronic_polymerases + reads_unspliced_transcripts
        
        # Check for NaNs and Infs
        assert torch.all(torch.isfinite(predicted_log_reads_exon)), f"NaN or Inf in predicted_log_reads_exon: {predicted_log_reads_exon=}"
        assert torch.all(torch.isfinite(predicted_reads_intron)), "NaN or Inf in predicted_reads_intron"
        assert torch.all(torch.isfinite(phi)), "NaN or Inf in phi"
        
        # Check for suspiciously large values
        assert torch.all(predicted_log_reads_exon < 100), "Suspiciously large log exon prediction"
        assert torch.all(predicted_reads_intron < 1e12), "Suspiciously large intron prediction"
        
        # Check for negative predictions
        assert torch.all(predicted_reads_intron >= 0), "Negative predicted intron reads"
        assert torch.all(safe_exp(predicted_log_reads_exon) >= 0), "Negative predicted exon reads (after exp)"

        

        return safe_exp(predicted_log_reads_exon), predicted_reads_intron, phi

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
                    for intron_index, intron_name in enumerate(self.intron_names):
                        parameter_data.append({'parameter_type': param_name,
                                               'intron_name': intron_name,
                                               'feature_name': feature_name,
                                               'value': param_value[feature_index, intron_index].item()})
            elif param_name in ('intercept_intron', 'log_phi_zero'):
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

    def forward(self, phi, coverage):
        loss_per_location = -torch.log(1 + phi.unsqueeze(2) * self.location_term)
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
                phi: torch.Tensor):
        loss_exon = self.loss_function_exon(predicted_reads_exon,
                                            reads_exon)
        loss_intron = self.loss_function_intron(predicted_reads_intron, reads_introns)
        loss_coverage = self.loss_function_coverage(phi, coverage)

        total_loss = loss_exon + loss_intron + loss_coverage
        return total_loss

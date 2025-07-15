from dataclasses import dataclass

import pandas as pd
import torch
import torch.nn as nn

EXP_INPUT_CLAMP = 40


@dataclass
class GeneData:
    gene_name: str
    intron_names: list[str]
    exon_reads: torch.Tensor
    intron_reads: torch.Tensor
    coverage: torch.Tensor

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

    def initialize_intercepts(self, gene_data: GeneData, library_sizes: torch.Tensor) -> None:
        with torch.no_grad():
            self.intercept_exon.data = torch.log(gene_data.exon_reads.mean() / library_sizes.mean())
            self.intercept_intron.data = torch.log(gene_data.intron_reads.mean(axis=0) / library_sizes.mean())

    def forward(self, design_matrix: torch.Tensor, log_library_sizes: torch.Tensor):
        alpha = self.alpha
        # if self.mask_alpha is not None:
        #     pass
        # alpha = self.alpha.clone()
        # alpha[self.mask_alpha.feature_index] = self.mask_alpha.value

        gene_expression_term = design_matrix @ alpha
        predicted_log_reads_exon = (self.intercept_exon + log_library_sizes + gene_expression_term)

        beta = self.beta
        gamma = self.gamma
        # TODO: handle masking for LRT

        speed_term = design_matrix @ beta
        splicing_term = design_matrix @ gamma

        phi = torch.sigmoid(self.log_phi_zero - speed_term - splicing_term)  # (num_samples, num_introns)

        intron_gene_expression_term = self.intercept_intron + log_library_sizes.unsqueeze(
            1) + gene_expression_term.unsqueeze(1)
        reads_intronic_polymerases = torch.exp(
            (intron_gene_expression_term + self.log_phi_zero - speed_term).clamp(max=EXP_INPUT_CLAMP,
                                                                                 min=-EXP_INPUT_CLAMP))
        reads_unspliced_transcripts = torch.exp(
            (intron_gene_expression_term + splicing_term).clamp(max=EXP_INPUT_CLAMP, min=-EXP_INPUT_CLAMP))
        predicted_reads_intron = reads_intronic_polymerases + reads_unspliced_transcripts

        if torch.isnan(predicted_log_reads_exon).any() or torch.isnan(predicted_reads_intron).any():
            current_params = dict(self.named_parameters())
            print(f"{current_params=}")
            assert False, "NaNs in predicted outputs"
        if torch.isinf(predicted_reads_intron).any():
            current_params = dict(self.named_parameters())
            print(f"{current_params=}")
            assert False, "Infs in predicted_reads_intron"
        return predicted_log_reads_exon, predicted_reads_intron, phi

    def get_param_df(self):
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
        self.loss_function_exon = nn.PoissonNLLLoss(log_input=True, full=True, reduction='sum')
        self.loss_function_intron = nn.PoissonNLLLoss(log_input=False, full=True, reduction='sum')
        self.loss_function_coverage = CoverageLoss(num_position_coverage=num_position_coverage)

    def forward(self,
                reads_exon: torch.Tensor,
                reads_introns: torch.Tensor,
                coverage: torch.Tensor,
                predicted_log_reads_exon: torch.Tensor,
                predicted_reads_intron: torch.Tensor,
                phi: torch.Tensor):
        loss_exon = self.loss_function_exon(predicted_log_reads_exon.clamp(max=EXP_INPUT_CLAMP, min=-EXP_INPUT_CLAMP),
                                            reads_exon)
        loss_intron = self.loss_function_intron(predicted_reads_intron, reads_introns)
        loss_coverage = self.loss_function_coverage(phi, coverage)

        total_loss = loss_exon + loss_intron + loss_coverage
        return total_loss

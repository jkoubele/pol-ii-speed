import torch
import torch.nn as nn
import torch.optim as optim
from tqdm import tqdm
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# %%

folder_aggregated_counts = Path('/cellfile/datapublic/jkoubele/data_pol_ii/drosophila_mutants/aggregated_counts')
folder_sample_annotation = Path('/cellfile/datapublic/jkoubele/data_pol_ii/drosophila_mutants/sample_annotation')
folder_intron_coverage = Path('/cellfile/datapublic/jkoubele/data_pol_ii/drosophila_mutants/intron_coverage')

df_genes = pd.read_csv(folder_aggregated_counts / 'exonic_reads.tsv', sep='\t').set_index('gene_id')
df_introns = pd.read_csv(folder_aggregated_counts / 'intron_counts.tsv', sep='\t').set_index('name')
df_library_size = pd.read_csv(folder_aggregated_counts / 'library_size_factors.tsv', sep='\t').set_index('sample_name')
df_design_matrix = pd.read_csv(folder_sample_annotation / 'design_matrix.tsv', sep='\t').set_index('sample_name')

assert (df_design_matrix.index == df_introns.columns).all()
assert (df_design_matrix.index == df_genes.columns).all()
assert (df_design_matrix.index == df_library_size.index).all()

sample_names: list[str] = df_design_matrix.index.to_list()

# %%
df_introns['gene_id'] = [intron_name.split('_')[0] for intron_name in df_introns.index]
intron_count_by_gene = {}
for gene_id in tqdm(df_introns['gene_id'].unique(), desc='Preparing intron read counts'):
    intron_count_by_gene[gene_id] = df_introns[df_introns['gene_id'] == gene_id].drop(columns=['gene_id'])

# %%
coverage_df_by_sample = {}
for sub_folder in folder_intron_coverage.iterdir():
    coverage_df_by_sample[sub_folder.name] = pd.read_parquet(
        sub_folder / 'rescaled_intron_coverages.parquet').set_index('intron_name')

# %%
X = torch.tensor(df_design_matrix.to_numpy(), dtype=torch.float32)
log_library_sizes = torch.log(torch.tensor(df_library_size['library_size_factor'].to_numpy(), dtype=torch.float32))


# %%

class Pol2Model(nn.Module):

    def __init__(self,
                 num_features: int,
                 intron_names: list[str],
                 exon_intercept_init: float,
                 intron_intercepts_init: dict[str, float]):
        super().__init__()
        self.intron_names = intron_names

        self.alpha = nn.Parameter(torch.zeros(num_features))
        # self.alpha = nn.Parameter(torch.randn(num_features) / 10)         
        self.intercept_exon = nn.Parameter(torch.tensor(exon_intercept_init))

        self.beta = nn.ParameterDict({
            intron: nn.Parameter(torch.randn(num_features) / 10)
            for intron in intron_names})
        self.gamma = nn.ParameterDict({
            intron: nn.Parameter(torch.randn(num_features) / 10)
            for intron in intron_names})

        self.intercept_intron = nn.ParameterDict({
            intron_name: nn.Parameter(torch.tensor(intercept_init))
            for intron_name, intercept_init in intron_intercepts_init.items()})
        self.log_phi_zero = nn.ParameterDict({
            intron: nn.Parameter(torch.tensor(0.0))
            for intron in intron_names})

    def forward(self, X, log_library_sizes):
        predicted_log_reads_exon = self.intercept_exon + log_library_sizes + X @ self.alpha

        predicted_reads_intron = {}
        for intron_name in self.intron_names:
            predicted_reads_intron[intron_name] = torch.exp(
                self.intercept_intron[intron_name] + log_library_sizes + X @ self.alpha + self.log_phi_zero[
                    intron_name] - X @ self.beta[intron_name]) + torch.exp(
                self.intercept_intron[intron_name] + log_library_sizes + X @ self.alpha + X @ self.gamma[intron_name])

        phi = {}
        for intron_name in self.intron_names:
            # Use softmax-like trick for numerical stability in division         
            speed_term = self.log_phi_zero[intron_name] - X @ self.beta[intron_name]
            splicing_term = X @ self.gamma[intron_name]
            max_term = torch.maximum(speed_term, splicing_term)
            phi[intron_name] = torch.exp(speed_term - max_term) / (
                    torch.exp(speed_term - max_term) + torch.exp(splicing_term - max_term))
        return predicted_log_reads_exon, predicted_reads_intron, phi


# %%
gene_name = 'FBgn0024352'
exon_read_count = df_genes.loc[gene_name]
intron_read_counts = intron_count_by_gene[gene_name]

intron_names: list[str] = intron_read_counts.index.to_list()

# TODO: remove genes / introns with zero counts 


# %%
data_exon_count = torch.tensor(df_genes.loc[gene_name].to_numpy(), dtype=torch.float32)
data_intron_count = {intron_name: torch.tensor(row.to_numpy(), dtype=torch.float32) for
                     intron_name, row in intron_read_counts.iterrows()}
data_coverage = {}
for intron_name in intron_names:
    data_coverage_list = []
    for sample_name in sample_names:
        density = coverage_df_by_sample[sample_name].loc[intron_name].to_numpy()
        density /= np.sum(density)
        data_coverage_list.append(density)
    data_coverage[intron_name] = torch.tensor(np.array(data_coverage_list))

# %%
model = Pol2Model(num_features=X.shape[1],
                  intron_names=intron_names,
                  exon_intercept_init=np.log(exon_read_count.mean()),
                  intron_intercepts_init={intron_name: np.log(row.mean()) for intron_name, row in
                                          intron_read_counts.iterrows()})


class CoverageLoss(nn.Module):

    def __init__(self, num_position_coverage: int = 100):
        super().__init__()
        self.locations = torch.linspace(start=1 / (2 * num_position_coverage),
                                        end=1 - 1 / (2 * num_position_coverage),
                                        steps=num_position_coverage).unsqueeze(0)

    def forward(self, phi, num_intronic_reads, coverage):
        loss_per_location = -torch.log(1 + phi.unsqueeze(1) - 2 * phi.unsqueeze(1) * self.locations)
        return torch.mean(loss_per_location * coverage * num_intronic_reads.unsqueeze(1))


loss_function_exon = nn.PoissonNLLLoss(log_input=True)
loss_function_intron = nn.PoissonNLLLoss(log_input=False)
loss_function_coverage = CoverageLoss(num_position_coverage=100)

optimizer = optim.Adam(model.parameters(), lr=0.01)
from time import time

t0 = time()
for epoch in tqdm(range(100)):
    optimizer.zero_grad()

    predicted_log_reads_exon, predicted_reads_intron, phi = model(X, log_library_sizes)  # Forward pass

    loss_exon = loss_function_exon(predicted_log_reads_exon, data_exon_count)
    loss_intron = torch.tensor(0.0)
    loss_coverage = torch.tensor(0.0)
    for intron_name in intron_names:
        loss_intron += loss_function_intron(predicted_reads_intron[intron_name], data_intron_count[intron_name])
        loss_coverage += loss_function_coverage(phi[intron_name], data_intron_count[intron_name],
                                                data_coverage[intron_name])
    loss = loss_exon + loss_intron + loss_coverage

    loss.backward()  # Backprop
    optimizer.step()
t1 = time()

print(f"{loss=}")
print(f"{loss_exon=}")
print(f"{loss_intron=}")
print(f"{loss_coverage=}")
print(f"{phi}")
print(f"{model.alpha=}")
print("----BETA----")
for key, value in model.beta.items():
    print(key, value)
print("----GAMMA----")
for key, value in model.gamma.items():
    print(key, value)

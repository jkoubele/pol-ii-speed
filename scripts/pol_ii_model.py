import torch
import torch.nn as nn
import torch.optim as optim
from tqdm import tqdm
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from torch.utils.tensorboard import SummaryWriter
import json
from visualize_read_density import plot_density, estimate_phi

# %%

project_name = 'human_astrocytes'

folder_aggregated_counts = Path('/cellfile/datapublic/jkoubele/data_pol_ii') / project_name / 'aggregated_counts'
folder_sample_annotation = Path('/cellfile/datapublic/jkoubele/data_pol_ii') / project_name / 'sample_annotation'
folder_intron_coverage = Path('/cellfile/datapublic/jkoubele/data_pol_ii') / project_name / 'intron_coverage'

df_genes = pd.read_csv(folder_aggregated_counts / 'exonic_reads.tsv', sep='\t').set_index('gene_id')
df_introns = pd.read_csv(folder_aggregated_counts / 'intron_counts.tsv', sep='\t').set_index('name')
df_library_size = pd.read_csv(folder_aggregated_counts / 'library_size_factors.tsv', sep='\t').set_index('sample_name')
df_design_matrix = pd.read_csv(folder_sample_annotation / 'design_matrix.tsv', sep='\t').set_index('sample_name')

with open('/cellfile/datapublic/jkoubele/reference_genomes/GRCh38.p14/protein_coding_genes.json') as protein_coding_json_file:
    protein_coding_genes = set(json.load(protein_coding_json_file))

assert (df_design_matrix.index == df_introns.columns).all()
assert (df_design_matrix.index == df_genes.columns).all()
assert (df_design_matrix.index == df_library_size.index).all()

sample_names: list[str] = df_design_matrix.index.to_list()

# %%
df_introns['gene_id'] = [intron_name.split('_')[0] for intron_name in df_introns.index]
df_introns = df_introns[[gene in protein_coding_genes for gene in df_introns['gene_id']]]
intron_count_by_gene = {}
gene_to_introns = defaultdict(list) 
for intron, gene in zip(df_introns.index, df_introns['gene_id']):
    gene_to_introns[gene].append(intron)
    
for gene_id in tqdm(df_introns['gene_id'].unique(), desc='Preparing intron read counts'):
    intron_count_by_gene[gene_id] = df_introns.loc[gene_to_introns[gene_id]].drop(columns=['gene_id'])

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
        self.intercept_exon = nn.Parameter(torch.tensor(exon_intercept_init))

        self.beta = nn.ParameterDict({
            intron: nn.Parameter(torch.zeros(num_features))
            for intron in intron_names})
        self.gamma = nn.ParameterDict({
            intron: nn.Parameter(torch.zeros(num_features))
            for intron in intron_names})

        self.intercept_intron = nn.ParameterDict({
            intron_name: nn.Parameter(torch.tensor(intercept_init))
            for intron_name, intercept_init in intron_intercepts_init.items()})
        self.log_phi_zero = nn.ParameterDict({
            intron: nn.Parameter(torch.tensor(0.0))
            for intron in intron_names})

    def forward(self, X, log_library_sizes):
        gene_expression_term =  X @ self.alpha        
        predicted_log_reads_exon = self.intercept_exon + log_library_sizes + gene_expression_term

        predicted_reads_intron = {}
        phi = {}
        for intron_name in self.intron_names:
            speed_term = X @ self.beta[intron_name]
            splicing_term = X @ self.gamma[intron_name]            
            
            phi[intron_name] = torch.sigmoid(self.log_phi_zero[intron_name] - speed_term - splicing_term)
            
            predicted_reads_intron[intron_name] = torch.exp(
                self.intercept_intron[intron_name] + self.log_phi_zero[intron_name] + log_library_sizes + gene_expression_term - speed_term) + torch.exp(
                self.intercept_intron[intron_name] + log_library_sizes + gene_expression_term + splicing_term)  
                    
        return predicted_log_reads_exon, predicted_reads_intron, phi


# %%

for gene_name, intron_read_counts in tqdm(intron_count_by_gene.items()):
    if gene_name not in df_genes.index:
        continue
    exon_read_count = df_genes.loc[gene_name]    
    if exon_read_count.mean() < 50 or exon_read_count.mean() > 2000:
        continue
    intron_names: list[str] = intron_read_counts.index.to_list()
    # if len(intron_names) != 1:
    #     continue
    data_exon_count = torch.tensor(df_genes.loc[gene_name].to_numpy(), dtype=torch.float32)
    data_intron_count = {intron_name: torch.tensor(row.to_numpy(), dtype=torch.float32) for
                         intron_name, row in intron_read_counts.iterrows()}
    data_coverage = {}
    found_good_example = False
    for intron_name in intron_names:
        if (data_intron_count[intron_name]==0).sum() > 3:
            continue
        if data_intron_count[intron_name].mean() > 150:
            continue
        lfc_intron = np.log(data_intron_count[intron_name][0:3].mean() / data_intron_count[intron_name][3:7].mean())
        lfc_exon = np.log(exon_read_count.iloc[0:3].mean() / exon_read_count.iloc[3:7].mean())
        if lfc_intron < lfc_exon:
            break
        data_coverage_list = []
        has_na = False
        for sample_name in sample_names:
            density = coverage_df_by_sample[sample_name].loc[intron_name].to_numpy()
            density /= np.sum(density)
            if np.isnan(density).any():
                has_na = True                
            data_coverage_list.append(density)
        if has_na:
            break
        data_coverage[intron_name] = torch.tensor(np.array(data_coverage_list))
        
        coverage_group_0 = data_coverage[intron_name][X[:,0]==0].mean(dim=0).numpy()
        coverage_group_1 = data_coverage[intron_name][X[:,0]==1].mean(dim=0).numpy()
    
        phi_0 = estimate_phi(coverage_group_0)
        phi_1 = estimate_phi(coverage_group_1)
        print(phi_0)
        if ((phi_0 > 0.05 and phi_0 < 0.95) and (phi_1 > 0.05 and phi_1 < 0.95)) and (phi_0 -0.05 > phi_1): 
            found_good_example = True
            break
    if found_good_example:
        break

# TODO: remove genes / introns with zero counts 


# %%
gene_name = 'ENSG00000100330'
exon_read_count = df_genes.loc[gene_name]
intron_read_counts = intron_count_by_gene[gene_name]

intron_names: list[str] = intron_read_counts.index.to_list()

data_exon_count = torch.tensor(df_genes.loc[gene_name].to_numpy(), dtype=torch.float32)
data_intron_count = {intron_name: torch.tensor(row.to_numpy(), dtype=torch.float32) for
                     intron_name, row in intron_read_counts.iterrows()}
data_coverage = {}
intron_names = ['ENSG00000100330_2']
for intron_name in intron_names:   
    data_coverage_list = []
    for sample_name in sample_names:
        density = coverage_df_by_sample[sample_name].loc[intron_name].to_numpy()
        density /= np.sum(density)
        data_coverage_list.append(density)
    data_coverage[intron_name] = torch.tensor(np.array(data_coverage_list))
    
# %%

coverage_group_0 = data_coverage['ENSG00000100330_2'][X[:,0]==0].mean(dim=0).numpy()
coverage_group_1 = data_coverage['ENSG00000100330_2'][X[:,0]==1].mean(dim=0).numpy()

phi_0 = estimate_phi(coverage_group_0)


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


loss_function_exon = nn.PoissonNLLLoss(log_input=True, full=True)
loss_function_intron = nn.PoissonNLLLoss(log_input=False, full=True)
loss_function_coverage = CoverageLoss(num_position_coverage=100)

optimizer = optim.Adam(model.parameters(), lr=1e-3)
from time import time

t0 = time()
for epoch in tqdm(range(500)):
    optimizer.zero_grad()
    
    tf0 = time()
    predicted_log_reads_exon, predicted_reads_intron, phi = model(X, log_library_sizes)  # Forward pass
    tf1 = time()
    

    loss_exon = loss_function_exon(predicted_log_reads_exon, data_exon_count)
    loss_intron = torch.tensor(0.0)
    loss_coverage = torch.tensor(0.0)
    for intron_name in intron_names:
        loss_intron += loss_function_intron(predicted_reads_intron[intron_name], data_intron_count[intron_name])
        tc0 = time()
        loss_coverage += loss_function_coverage(phi[intron_name], data_intron_count[intron_name],
                                                data_coverage[intron_name])
        tc1 = time()
    loss = loss_exon + loss_intron + 100*loss_coverage
    tl0=time()
    loss.backward()  # Backprop
    tl1=time()
    to0 = time()
    optimizer.step()
    to1 = time()
t1 = time()
print('')
# print(f"{(t1-t0)*1000=}")
# print(f"Forward pass: {(tf1-tf0)*1000=}")
# print(f"Coverage loss: {(tc1-tc0)*1000=}")
# print(f"Backward pass: {(tl1-tl0)*1000=}")
# print(f"Optimizer time: {(to1-to0)*1000=}")

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

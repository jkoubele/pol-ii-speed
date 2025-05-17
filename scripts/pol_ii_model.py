import torch
import torch.nn as nn
import torch.optim as optim
from tqdm import tqdm
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

class Pol2Model(nn.Module):

    def __init__(self, num_features: int):
        super().__init__()

        self.alpha = nn.Parameter(torch.randn(num_features) /10) 
        self.beta = nn.Parameter(torch.randn(num_features) / 10 )  
        self.gamma = nn.Parameter(torch.randn(num_features) / 10)
        
        # self.alpha = nn.Parameter(torch.zeros(num_features))
        # self.beta = nn.Parameter(torch.zeros(num_features))
        # self.gamma = nn.Parameter(torch.zeros(num_features))

        self.intercept_exon = nn.Parameter(torch.tensor(-11.0))
        self.intercept_intron = nn.Parameter(torch.tensor(-13.0))
        self.phi_zero = nn.Parameter(torch.tensor(0.0))

    def forward(self, Z, log_library_sizes):
        predicted_log_reads_exon = self.intercept_exon + log_library_sizes + Z @ self.alpha
        predicted_reads_intron = torch.exp(
            self.intercept_intron + log_library_sizes + Z @ self.alpha + self.phi_zero - Z @ self.beta) + torch.exp(
            self.intercept_intron + log_library_sizes + Z @ self.alpha + Z @ self.gamma)

        # Use softmax-like trick for numerical stability in division         
        speed_term = self.phi_zero - Z @ self.beta
        splicing_term = Z @ self.gamma
        max_term = torch.maximum(speed_term, splicing_term)
        phi = torch.exp(speed_term - max_term) / (
                torch.exp(speed_term - max_term) + torch.exp(splicing_term - max_term))
        return predicted_log_reads_exon, predicted_reads_intron, phi


input_folder = Path('/home/jakub/Desktop/data_pol_ii/human_astrocytes/aggregated_counts')
df_genes = pd.read_csv(input_folder / 'exonic_reads.tsv', sep='\t').set_index('gene_id')
df_introns = pd.read_csv(input_folder / 'intron_counts.tsv', sep='\t').set_index('name')

annotation = pd.DataFrame(data = {'young': [1,1,1,0,0,0,0], 'old': [0,0,0,1,1,1,1]},
                          index=['GSM1901309', 'GSM1901310', 'GSM1901311', 'GSM1901316', 'GSM1901318', 'GSM1901319', 'GSM1901320'])

intron_name = 'ENSG00000009335_1'
gene_name = intron_name.split('_')[0]

model = Pol2Model(2)
num_position_coverage = 100
optimizer = optim.Adam(model.parameters(), lr=0.01)

Z = torch.tensor(annotation.to_numpy(), dtype=torch.float32)

intron_coverage_folder = Path('/home/jakub/Desktop/data_pol_ii/human_astrocytes/intron_coverage')
coverage_df_by_sample = {}
for sub_folder in intron_coverage_folder.iterdir():    
    coverage_df_by_sample[sub_folder.name] = pd.read_parquet(sub_folder / 'rescaled_intron_coverages.parquet').set_index('intron_name')
    
#%%
log_library_sizes = torch.tensor(7 * [16], dtype=torch.float32)

output = model(Z, log_library_sizes)

data_exon_count = torch.tensor(df_genes.loc[gene_name].to_numpy(), dtype=torch.float32)
data_intron_count = torch.tensor(df_introns.loc[intron_name].to_numpy(), dtype=torch.float32)

density = torch.linspace(start=10, end=5, steps=num_position_coverage)
density /= density.sum()

data_coverage = density.repeat(7, 1)
data_coverage_list = []
for sample_name in  df_genes.loc[gene_name].index:
    x = coverage_df_by_sample[sample_name].loc[intron_name].to_numpy()
    x /= np.sum(x)
    data_coverage_list.append(x)
    
data_coverage = torch.tensor(data_coverage_list)

#%%


class CoverageLoss(nn.Module):

    def __init__(self, num_position_coverage: int = 10):
        super().__init__()
        self.locations = torch.linspace(start=1 / (2 * num_position_coverage),
                                        end=1 - 1 / (2 * num_position_coverage),
                                        steps=num_position_coverage).unsqueeze(0)

    def forward(self, phi, num_intronic_reads, coverage):
        loss_per_location = -torch.log(1 + phi.unsqueeze(1) - 2 * phi.unsqueeze(1) * self.locations)
        return torch.mean(loss_per_location * coverage * data_intron_count.unsqueeze(1))


loss_function_exon = nn.PoissonNLLLoss(log_input=True)
loss_function_intron = nn.PoissonNLLLoss(log_input=False)
loss_function_coverage = CoverageLoss(num_position_coverage)

for epoch in tqdm(range(10000)):
    optimizer.zero_grad()

    predicted_log_reads_exon, predicted_reads_intron, phi = model(Z, log_library_sizes)  # Forward pass
    loss_exon = loss_function_exon(predicted_log_reads_exon, data_exon_count)
    loss_intron = loss_function_intron(predicted_reads_intron, data_intron_count)
    loss_coverage = loss_function_coverage(phi, data_intron_count, data_coverage)
    loss = loss_exon + loss_intron + loss_coverage

    loss.backward()  # Backprop
    optimizer.step()

print(f"{loss=}")
print(f"{loss_exon=}")
print(f"{loss_intron=}")
print(f"{loss_coverage=}")
print(f"{model.alpha[1] - model.alpha[0]=}")
print(f"{model.beta[1] - model.beta[0]=}")
print(f"{model.gamma[1] - model.gamma[0]=}")
print(f"{phi}")

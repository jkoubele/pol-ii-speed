import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.tensorboard import SummaryWriter

import pickle
from tqdm import tqdm
from pathlib import Path
from pol_ii_model import Pol2Model, CoverageLoss, GeneData, GeneDataWithSolution
import math

project_path = Path('/home/jakub/Desktop/dev-pol-ii-analysis/data/drosophila_mutants')
with open(project_path / 'train_data' / 'data_train_with_solutions.pkl', 'rb') as file:
    data_train_with_solutions = pickle.load(file)

train_data_with_solutions = data_train_with_solutions['gene_data_with_solutions']
design_matrix = data_train_with_solutions['design_matrix'].float()
log_library_size = torch.log(data_train_with_solutions['library_size']).float()

gene_data_with_solutions = train_data_with_solutions[0]
gene_data = gene_data_with_solutions.gene_data

model = Pol2Model(num_features=design_matrix.shape[1],
                  intron_names=gene_data.intron_names,
                  exon_intercept_init=float(torch.log(gene_data.exon_reads.mean())),
                  intron_intercepts_init={intron_name: float(torch.log(intron_reads.mean()))
                                          for intron_name, intron_reads in gene_data.intron_reads.items()})

loss_function_exon = nn.PoissonNLLLoss(log_input=True, full=True)
loss_function_intron = nn.PoissonNLLLoss(log_input=False, full=True)
loss_function_coverage = CoverageLoss()

# optimizer = optim.Adam(model.parameters(), lr=2 * 1e-2)
optimizer = optim.LBFGS(model.parameters(),
                        lr=1.0,
                        max_iter=100,
                        tolerance_change=1e-9,
                        tolerance_grad=1e-7,
                        history_size=100,
                        line_search_fn=None)

optimal_intercept_exon = math.log(gene_data_with_solutions.mu_1)
optimal_alpha = math.log(gene_data_with_solutions.mu_2 / gene_data_with_solutions.mu_1)

writer = SummaryWriter()
previous_loss = None


def closure():
    optimizer.zero_grad()
    predicted_log_reads_exon, predicted_reads_intron, phi = model(design_matrix, log_library_size)

    loss_exon = loss_function_exon(predicted_log_reads_exon, gene_data.exon_reads)
    loss_intron = torch.zeros(1)
    loss_coverage = torch.zeros(1)

    for intron_name in gene_data.intron_names:
        loss_intron += loss_function_intron(
            predicted_reads_intron[intron_name],
            gene_data.intron_reads[intron_name]
        )
        loss_coverage += loss_function_coverage(
            phi[intron_name],
            gene_data.intron_reads[intron_name],
            gene_data.coverage_density[intron_name]
        )

    loss = loss_exon + loss_intron + loss_coverage
    # loss = loss_exon

    loss.backward()
    return loss


intron_of_interest = gene_data.intron_names[7]
intron_of_interest = 'FBgn0000008_3'

for epoch in tqdm(range(50)):

    loss = optimizer.step(closure)

    # Evaluate loss again (optional, for logging):
    with torch.no_grad():
        predicted_log_reads_exon, predicted_reads_intron, phi = model(design_matrix, log_library_size)

        loss_exon = loss_function_exon(predicted_log_reads_exon, gene_data.exon_reads)
        loss_intron = torch.zeros(1)
        loss_coverage = torch.zeros(1)

        for intron_name in gene_data.intron_names:
            loss_intron += loss_function_intron(
                predicted_reads_intron[intron_name],
                gene_data.intron_reads[intron_name]
            )
            loss_coverage += loss_function_coverage(
                phi[intron_name],
                gene_data.intron_reads[intron_name],
                gene_data.coverage_density[intron_name]
            )

        total_loss = loss_exon + loss_intron + loss_coverage

        writer.add_scalar('Total loss', total_loss, epoch)
        writer.add_scalar('Loss exons', loss_exon, epoch)
        writer.add_scalar('Loss introns', loss_intron, epoch)
        writer.add_scalar('Loss coverage', loss_coverage, epoch)

        writer.add_scalar('Intercept exon', model.intercept_exon, epoch)
        writer.add_scalar('Intercept first intron', model.intercept_intron[intron_of_interest], epoch)
        writer.add_scalar('Alpha', model.alpha, epoch)
        writer.add_scalar('Beta (first intron)', model.beta[intron_of_interest], epoch)
        writer.add_scalar('Gamma (first intron)', model.gamma[intron_of_interest], epoch)
        writer.add_scalar('Log(phi_0) (first intron)', model.log_phi_zero[intron_of_interest], epoch)
        writer.add_scalar('Phi - group 1 (first intron)', phi[intron_of_interest][0], epoch)
        writer.add_scalar('Phi - group 2 (first intron)', phi[intron_of_interest][-1], epoch)

        # for name, param in model.named_parameters():
        #     if param.grad is not None:
        #         grad_norm = param.grad.norm().item()
        #         print(f"{name}: gradient norm = {grad_norm:.3e}")

        if previous_loss is not None:
            loss_difference = previous_loss - total_loss
            writer.add_scalar('Loss difference', loss_difference, epoch)
        previous_loss = total_loss

# %%
import torch
from torch.func import functional_call, hessian

# Extract parameters and buffers from the model
params = dict(model.named_parameters())
buffers = dict(model.named_buffers())


def loss_fn(params):
    # Perform a functional call with the current parameters and buffers
    outputs = functional_call(model, (params, buffers), (design_matrix, log_library_size))
    predicted_log_reads_exon, predicted_reads_intron, phi = outputs

    loss_exon = loss_function_exon(predicted_log_reads_exon, gene_data.exon_reads)
    loss_intron = torch.zeros(1)
    loss_coverage = torch.zeros(1)
    for intron_name in gene_data.intron_names:
        loss_intron += loss_function_intron(
            predicted_reads_intron[intron_name],
            gene_data.intron_reads[intron_name]
        )
        loss_coverage += loss_function_coverage(
            phi[intron_name],
            gene_data.intron_reads[intron_name],
            gene_data.coverage_density[intron_name]
        )
    total_loss = loss_exon + loss_intron + loss_coverage
    return total_loss


# Compute the Hessian
H = hessian(loss_fn)(params)

param_names = []
param_sizes = []
for name, p in params.items():
    param_names.append(name)
    param_sizes.append(p.numel())

index_map = {}
start = 0
for name, size in zip(param_names, param_sizes):
    index_map[name] = (start, start + size)
    start += size
total_params = start

H_matrix = torch.zeros((total_params, total_params))

for name1, block_row in H.items():
    idx1_start, idx1_end = index_map[name1]
    for name2, block in block_row.items():
        idx2_start, idx2_end = index_map[name2]
        H_matrix[idx1_start:idx1_end, idx2_start:idx2_end] = block.reshape(
            idx1_end - idx1_start, idx2_end - idx2_start
        )

H_inv = torch.linalg.pinv(H_matrix)

standard_errors = torch.sqrt(torch.diagonal(H_inv))
flat_params = torch.cat([p.flatten() for p in params.values()])
z_scores = flat_params / standard_errors
from scipy.stats import norm

p_values = 2 * (1 - norm.cdf(torch.abs(z_scores).detach().numpy()))

import pandas as pd

param_labels = []
for name, p in params.items():
    for idx in range(p.numel()):
        param_labels.append(f"{name}[{idx}]")

df = pd.DataFrame({
    'Parameter': param_labels,
    'Estimate': flat_params.detach().numpy(),
    'StdError': standard_errors.detach().numpy(),
    'Z': z_scores.detach().numpy(),
    'p-value': p_values
})

# %%

import numpy as np
import matplotlib.pyplot as plt

group_mask_1 = torch.tensor(6 * [True] + 6 * [False])
group_mask_2 = torch.tensor(6 * [False] + 6 * [True])

# Debug exon parameters

exon_reads = gene_data.exon_reads

exon_reads_1 = exon_reads[group_mask_1]
exon_reads_2 = exon_reads[group_mask_2]

optimal_intercept_exon = torch.log(exon_reads_1.mean()).item()
optimal_alpha = torch.log(exon_reads_2.mean() / exon_reads_1.mean()).item()

fitted_intercept_exon = model.intercept_exon.item()
fitted_alpha = model.alpha.item()

print('')
print(f"{optimal_intercept_exon=}")
print(f"{fitted_intercept_exon=}")
print(50 * '-')
print(f"{optimal_alpha=}")
print(f"{fitted_alpha=}")

# Beta and gamma
mu_1 = gene_data_with_solutions.mu_1
mu_2 = gene_data_with_solutions.mu_2

nu_1 = gene_data_with_solutions.nu_1[intron_of_interest]
nu_2 = gene_data_with_solutions.nu_2[intron_of_interest]

phi_1 = gene_data_with_solutions.phi_1[intron_of_interest]
phi_2 = gene_data_with_solutions.phi_2[intron_of_interest]

optimal_beta = math.log(mu_2 / mu_1) - math.log(nu_2 / nu_1) - math.log(phi_2 / phi_1)
fitted_beta = model.beta[intron_of_interest].item()

optimal_gamma = -math.log(mu_2 / mu_1) + math.log(nu_2 / nu_1) + math.log((1 - phi_2) / (1 - phi_1))
fitted_gamma = model.gamma[intron_of_interest].item()

print(50 * '-')
print(f"{optimal_beta=}")
print(f"{fitted_beta=}")

print(50 * '-')
print(f"{optimal_gamma=}")
print(f"{fitted_gamma=}")

optimal_intercept_intron = math.log((1 - phi_1) * nu_1)
fitted_intercept_intron = model.intercept_intron[intron_of_interest].item()

print(50 * '-')
print(f"{optimal_intercept_intron=}")
print(f"{fitted_intercept_intron=}")

optimal_log_phi_zero = math.log(phi_1 / (1 - phi_1))
fitted_log_phi_zero = model.log_phi_zero[intron_of_interest].item()

print(50 * '-')
print(f"{optimal_log_phi_zero=}")
print(f"{fitted_log_phi_zero=}")

predicted_log_reads_exon_1 = predicted_log_reads_exon[0]
predicted_log_reads_exon_2 = predicted_log_reads_exon[-1]

# Debug phi
phi_1_fitted = phi[intron_of_interest][0].item()
phi_2_fitted = phi[intron_of_interest][-1].item()

phi_1_optimal = gene_data_with_solutions.phi_1[intron_of_interest]
phi_2_optimal = gene_data_with_solutions.phi_2[intron_of_interest]

print(50 * '-')
print(f"{phi_1_optimal=}")
print(f"{phi_1_fitted=}")
print(50 * '-')
print(f"{phi_2_optimal=}")
print(f"{phi_2_fitted=}")

phi_1 = phi[intron_of_interest][group_mask_1]
phi_2 = phi[intron_of_interest][group_mask_2]

intron_reads_first_intron = gene_data.intron_reads[intron_of_interest]
coverage_density_first_intron = gene_data.coverage_density[intron_of_interest]

intron_reads_first_intron_1 = intron_reads_first_intron[group_mask_1]
intron_reads_first_intron_2 = intron_reads_first_intron[group_mask_2]

coverage_density_first_intron_1 = coverage_density_first_intron[group_mask_1]
coverage_density_first_intron_2 = coverage_density_first_intron[group_mask_2]

all_phi_experimental = np.linspace(0.01, 0.99, 99)
losses_experimental = []
for phi_experimental in all_phi_experimental:
    losses_experimental.append(loss_function_coverage(torch.tensor(6 * [phi_experimental]),
                                                      intron_reads_first_intron_1,
                                                      coverage_density_first_intron_1))

best_phi_1_experimental = all_phi_experimental[np.argmin(losses_experimental)]

plt.plot(all_phi_experimental, losses_experimental)
plt.show()

# %%
from read_location_model import plot_density

intron_reads_1 = gene_data.intron_reads[intron_of_interest][group_mask_1]
intron_reads_2 = gene_data.intron_reads[intron_of_interest][group_mask_2]

densities_1 = gene_data.coverage_density[intron_of_interest][group_mask_1]
densities_2 = gene_data.coverage_density[intron_of_interest][group_mask_2]

aggregate_density_1 = (densities_1 * intron_reads_1.unsqueeze(1)).sum(dim=0)
aggregate_density_1 /= aggregate_density_1.sum()

aggregate_density_2 = (densities_2 * intron_reads_2.unsqueeze(1)).sum(dim=0)
aggregate_density_2 /= aggregate_density_2.sum()

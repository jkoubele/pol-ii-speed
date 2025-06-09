import torch
import torch.nn as nn
import torch.optim as optim
from torch.func import functional_call, hessian
from torch.utils.tensorboard import SummaryWriter
import pandas as pd

import pickle
from tqdm import tqdm, trange
from pathlib import Path
from pol_ii_model import Pol2Model, CoverageLoss, GeneData, Pol2TotalLoss
import math

project_path = Path('/home/jakub/Desktop/dev-pol-ii-analysis/data/drosophila_mutants')
with open(project_path / 'train_data' / 'data_train_with_solutions.pkl', 'rb') as file:
    data_train_with_solutions = pickle.load(file)

# gene_data, design_matrix and library_sizes are input to model training
gene_data = data_train_with_solutions['gene_data_with_solutions'][0].gene_data
design_matrix = pd.DataFrame(data={'genotype_rp2': 6 * [0] + 6 * [1],
                                   'dummy': 3 * [1] + 6 * [0] + 3 * [1]})
library_sizes = torch.ones(len(design_matrix)).float()



def train_model(gene_data: GeneData,
                design_matrix: pd.DataFrame,
                library_sizes: torch.tensor,
                pol_2_total_loss: Pol2TotalLoss,
                max_epochs=100,
                use_cuda=True) -> Pol2Model:
    device = torch.device("cuda" if torch.cuda.is_available() and use_cuda else "cpu")

    X = torch.tensor(design_matrix.to_numpy()).float().to(device)
    log_library_size = torch.log(library_sizes).to(device)
    gene_data = gene_data.to(device)
    model = Pol2Model(num_features=X.shape[1],
                      intron_names=gene_data.intron_names,
                      exon_intercept_init=float(torch.log((gene_data.exon_reads / library_sizes).mean())),
                      intron_intercepts_init={intron_name: float(torch.log((intron_reads / library_sizes).mean()))
                                              for intron_name, intron_reads in gene_data.intron_reads.items()})
    model = model.to(device)

    optimizer = optim.LBFGS(model.parameters(),
                            lr=1.0,
                            max_iter=100,
                            tolerance_change=1e-9,
                            tolerance_grad=1e-7,
                            history_size=100,
                            line_search_fn=None)

    def closure():
        optimizer.zero_grad()
        predicted_log_reads_exon, predicted_reads_intron, phi = model(X, log_library_size)
        loss = pol_2_total_loss(gene_data, predicted_log_reads_exon, predicted_reads_intron, phi)
        loss.backward()
        return loss

    for epoch in trange(max_epochs):
        loss = optimizer.step(closure)
    return model


# if __name__ == "__main__":
pol_2_total_loss = Pol2TotalLoss()
model = train_model(gene_data, design_matrix, library_sizes, pol_2_total_loss)

# device = next(model.parameters()).device  # detect model device
# X = torch.tensor(design_matrix.to_numpy()).float().to(device)
# log_library_size = torch.log(library_sizes).to(device)
#
# # Extract model parameters and buffers
# params = dict(model.named_parameters())
# buffers = dict(model.named_buffers())
#
# # Define the loss function (same as training)
# loss_function_exon = nn.PoissonNLLLoss(log_input=True, full=True)
# loss_function_intron = nn.PoissonNLLLoss(log_input=False, full=True)
# loss_function_coverage = CoverageLoss().to(device)
#
#
# def loss_fn(params):
#     outputs = functional_call(model, (params, buffers), (X, log_library_size))
#     predicted_log_reads_exon, predicted_reads_intron, phi = outputs
#
#     loss_exon = loss_function_exon(predicted_log_reads_exon, gene_data.exon_reads.to(device))
#     loss_intron = torch.zeros(1, device=device)
#     loss_coverage = torch.zeros(1, device=device)
#     for intron_name in gene_data.intron_names:
#         loss_intron += loss_function_intron(
#             predicted_reads_intron[intron_name],
#             gene_data.intron_reads[intron_name].to(device)
#         )
#         loss_coverage += loss_function_coverage(
#             phi[intron_name],
#             gene_data.intron_reads[intron_name].to(device),
#             gene_data.coverage_density[intron_name].to(device)
#         )
#     return loss_exon + loss_intron + loss_coverage
#
#
# # Compute the Hessian
# H = hessian(loss_fn)(params)
#
# # Flatten the Hessian into a 2D matrix
# param_names = []
# param_sizes = []
# for name, p in params.items():
#     param_names.append(name)
#     param_sizes.append(p.numel())
#
# index_map = {}
# start = 0
# for name, size in zip(param_names, param_sizes):
#     index_map[name] = (start, start + size)
#     start += size
# total_params = start
#
# H_matrix = torch.zeros((total_params, total_params), device=device)
#
# for name1, block_row in H.items():
#     idx1_start, idx1_end = index_map[name1]
#     for name2, block in block_row.items():
#         idx2_start, idx2_end = index_map[name2]
#         H_matrix[idx1_start:idx1_end, idx2_start:idx2_end] = block.reshape(
#             idx1_end - idx1_start, idx2_end - idx2_start
#         )

import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.optim as optim
from scipy import stats
from torch.func import functional_call, hessian
from tqdm import trange

from pol_ii_model import Pol2Model, GeneData, Pol2TotalLoss


def train_model(gene_data: GeneData,
                X: torch.Tensor,
                log_library_size: torch.Tensor,
                pol_2_total_loss: Pol2TotalLoss,
                device: str,
                max_epochs=100
                ) -> Pol2Model:
    X = torch.tensor(design_matrix.to_numpy()).float().to(device)
    log_library_size = torch.log(library_sizes).to(device)
    model = Pol2Model(num_features=X.shape[1],
                      intron_names=gene_data.intron_names,
                      exon_intercept_init=float(torch.log((gene_data.exon_reads / library_sizes).mean())),
                      intron_intercepts_init={intron_name: float(torch.log((intron_reads / library_sizes).mean()))
                                              for intron_name, intron_reads in gene_data.intron_reads.items()})
    model = model.to(device)
    pol_2_total_loss = pol_2_total_loss.to(device)
    gene_data = gene_data.to(device)

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


project_path = Path('/home/jakub/Desktop/dev-pol-ii-analysis/data/drosophila_mutants')
with open(project_path / 'train_data' / 'data_train_with_solutions.pkl', 'rb') as file:
    data_train_with_solutions = pickle.load(file)

# gene_data, design_matrix and library_sizes are input to model training
gene_data = data_train_with_solutions['gene_data_with_solutions'][0].gene_data

# design_matrix = pd.DataFrame(data={'genotype_rp2': 6 * [0] + 6 * [1],
#                                    'dummy': 3 * [1] + 6 * [0] + 3 * [1]})

design_matrix = pd.DataFrame(data={'genotype_rp2': 6 * [0] + 6 * [1]})

library_sizes = torch.ones(len(design_matrix)).float()
use_cuda = True

# Above is the input to the 'main' per-gene function
feature_names = design_matrix.columns.to_list()
device = torch.device("cuda" if torch.cuda.is_available() and use_cuda else "cpu")

X = torch.tensor(design_matrix.to_numpy()).float().to(device)
log_library_size = torch.log(library_sizes).to(device)
pol_2_total_loss = Pol2TotalLoss().to(device)

model = train_model(gene_data, X, log_library_size, pol_2_total_loss, device)

model_parameters = dict(model.named_parameters())


def loss_by_model_parameters(model_parameters):
    outputs = functional_call(model, model_parameters, (X, log_library_size))
    predicted_log_reads_exon, predicted_reads_intron, phi = outputs
    return pol_2_total_loss(gene_data, predicted_log_reads_exon, predicted_reads_intron, phi)


# Compute the Hessian
H = hessian(loss_by_model_parameters)(model_parameters)

# Flatten the Hessian into a 2D matrix
parameter_names = []
parameter_sizes = []
for name, p in model_parameters.items():
    parameter_names.append(name)
    parameter_sizes.append(p.numel())

index_map = {}
start = 0
for name, size in zip(parameter_names, parameter_sizes):
    index_map[name] = (start, start + size)
    start += size
total_params = start

H_matrix = torch.zeros((total_params, total_params), device=device)

for name1, block_row in H.items():
    idx1_start, idx1_end = index_map[name1]
    for name2, block in block_row.items():
        idx2_start, idx2_end = index_map[name2]
        H_matrix[idx1_start:idx1_end, idx2_start:idx2_end] = block.reshape(
            idx1_end - idx1_start, idx2_end - idx2_start
        )

H_inv = torch.linalg.pinv(H_matrix).detach()

# Standard errors (sqrt of diagonal elements)
standard_errors = torch.sqrt(torch.diag(H_inv)).cpu().numpy()

# %%
parameter_data: list[dict] = []
for name, size in zip(parameter_names, parameter_sizes):
    if '.' in name:
        parameter_type = name.split('.', 1)[0]
        intron_name = name.split('.', 1)[1]
    else:
        parameter_type = name.split('.', 1)[0]
        intron_name = None

    parameter_value_list = [model_parameters[name].item()] if model_parameters[name].numel() == 1 \
        else model_parameters[name].tolist()
    for feature_index in range(size):
        feature_name = None
        if parameter_type in ('alpha', 'beta', 'gamma'):
            feature_name = feature_names[feature_index]
        parameter_data.append({'parameter_type': parameter_type,
                               'intron_name': intron_name,
                               'feature_name': feature_name,
                               'value': parameter_value_list[feature_index]})

df_param = pd.DataFrame(data=parameter_data)
df_param['SE'] = standard_errors
df_param['z_value'] = df_param['value'] / standard_errors
df_param['p_value'] = 2 * (1 - stats.norm.cdf(np.abs(df_param['z_value'])))

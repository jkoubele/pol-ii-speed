import numpy as np
import pandas as pd
import cvxpy as cp
import torch
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.special import expit
import pickle

from pol_ii_speed_modeling.pol_ii_model import (
    GeneData, DatasetMetadata, Pol2Model, Pol2TotalLoss
)


from pol_ii_speed_modeling.load_dataset import load_dataset_from_results_folder

with open("gene_data_list.pkl", "rb") as f:
    gene_data_list = pickle.load(f)
    
with open("dataset_metadata.pkl", "rb") as f:
    dataset_metadata = pickle.load(f)
    
# dataset_metadata, gene_data_list = load_dataset_from_results_folder(
#     results_folder=Path('/cellfile/datapublic/jkoubele/data_pol_ii/mouse_myocardium/results/'),
#     log_output_folder=Path('/cellfile/datapublic/jkoubele/data_pol_ii/mouse_myocardium/results/logs'),
#     gene_names_file_name='test.csv')
#%%

gene_data = gene_data_list[0]
design_matrix = dataset_metadata.design_matrix.numpy()

n, p = design_matrix.shape

exon_reads = gene_data.exon_reads.numpy()  # shape (n,)
intron_reads = gene_data.intron_reads[:,2].numpy() # shape (n,)
alpha = cp.Variable(p) # shape (p,)
intercept_exon = cp.Variable()

predicted_log_reads_exon = intercept_exon + design_matrix @ alpha
exon_loss = cp.sum(cp.exp(predicted_log_reads_exon) - cp.multiply(exon_reads, predicted_log_reads_exon))


beta = cp.Variable(p) # shape (p,)
gamma = cp.Variable(p) # shape (p,)
intercept_intron = cp.Variable()
theta = cp.Variable()

custom_value = cp.sum(2 * alpha) 



objective = cp.Minimize(exon_loss)
problem = cp.Problem(objective)
ret = problem.solve()

ret_stats = problem.solver_stats

print(f"{problem.status=}")

reads_1 = exon_reads[(design_matrix==0).flatten()]
reads_2 = exon_reads[(design_matrix==1).flatten()]


# ret = problem.solve()





# coverage = gene_data.coverage[1,0,:].numpy()
# num_position_coverage = len(coverage)
# locations = np.linspace(start=1 / (2 * num_position_coverage),
#                                    stop=1 - 1 / (2 * num_position_coverage),
#                                    num=num_position_coverage)
# location_term = 1 - 2 * locations



# logit_linspace = np.linspace(-10, 10, 100)
# pi_linspace = np.linspace(0.01, 0.99, 100)


# lossess_by_logit = [np.sum(-np.log(1 + expit(logit) * location_term) * coverage) for logit in logit_linspace]
# lossess_by_pi = [np.sum(-np.log(1 + pi * location_term) * coverage) for pi in pi_linspace]
# # total_loss = np.sum(coverage * loss_per_location)0

# plt.plot(coverage)
# plt.title('coverage')
# plt.show()

# plt.plot(pi_linspace, lossess_by_pi)
# plt.title('lossess_by_pi')
# plt.show()

# plt.plot(logit_linspace, lossess_by_logit)
# plt.title('Loss by logit')
# plt.show()
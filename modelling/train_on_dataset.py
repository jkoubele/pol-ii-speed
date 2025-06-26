import pickle
from pathlib import Path
import torch

import pickle
from pathlib import Path

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import scipy
import torch
import torch.nn as nn
import torch.optim as optim
from scipy import stats
from torch.func import functional_call, hessian
from tqdm import trange, tqdm

from pol_ii_model import Pol2Model, GeneData, Pol2TotalLoss
from read_location_model import estimate_phi

from train import fit_analytical_solution, get_param_df, compute_hessian_matrix, add_wald_test_results, train_model



project_path = Path('/cellfile/datapublic/jkoubele/data_pol_ii/drosophila_mutants/')
with open(project_path / 'train_data' / 'data_train_2_groups.pkl', 'rb') as file:
    dataset = pickle.load(file)
 
#%%
data_train = dataset['gene_data']
design_matrix = dataset['design_matrix']
library_sizes = dataset['library_size']
library_sizes = library_sizes.to(dtype=torch.float32)

feature_names = ['rp2']

use_cuda = True
device = torch.device("cuda" if torch.cuda.is_available() and use_cuda else "cpu")

pol_2_total_loss = Pol2TotalLoss().to(device)
log_library_size = torch.log(library_sizes).to(device)

results = []
for gene_data in tqdm(data_train):   
    model_analytical = fit_analytical_solution(gene_data, design_matrix, library_sizes, device)
    df_param_analytical = get_param_df(model_analytical, feature_names)

    hessian_matrix_analytical = compute_hessian_matrix(model=model_analytical,
                                                       pol_2_total_loss=pol_2_total_loss,
                                                       X=design_matrix,
                                                       log_library_size=log_library_size,
                                                       gene_data=gene_data,
                                                       device=device)

    df_param_analytical = add_wald_test_results(df_param_analytical, hessian_matrix_analytical)
    df_param_analytical['gene'] = gene_data.intron_names[0].split('_')[0]
    
    model_numerical = fit_n(gene_data, design_matrix, library_sizes, device)
    df_param_analytical = get_param_df(model_analytical, feature_names)

    hessian_matrix_analytical = compute_hessian_matrix(model=model_analytical,
                                                       pol_2_total_loss=pol_2_total_loss,
                                                       X=design_matrix,
                                                       log_library_size=log_library_size,
                                                       gene_data=gene_data,
                                                       device=device)

    df_param_analytical = add_wald_test_results(df_param_analytical, hessian_matrix_analytical)
    df_param_analytical['gene'] = gene_data.intron_names[0].split('_')[0]
    
    
    results.append(df_param_analytical)
        
   
    
df_all = pd.concat(results, ignore_index=True)

#%%

b = df_all[df_all['parameter_type']=='beta']
g = df_all[df_all['parameter_type']=='gamma']
b = b[b['identifiable']]
b['FDR'] = stats.false_discovery_control(b['p_value_wald'])
# b_significant = b[b['p_value_wald'] < 0.05]
b_significant = b[b['FDR'] < 0.05]

b_up = b_significant[b_significant['value']>0]
b_down = b_significant[b_significant['value']<0]
  
print(f"num insignificant: {len(b) - len(b_significant)}")  
print(f'num up: {len(b_up)}')
print(f'num down: {len(b_down)}')
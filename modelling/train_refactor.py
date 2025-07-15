from pathlib import Path
import pickle
import torch
from torch import nn
from torch import optim
import pandas as pd
import numpy as np

from tqdm import trange
from dataclasses import dataclass
from pol_ii_model_refactor import GeneData, DatasetMetadata, Pol2TotalLoss, Pol2Model
from torch.func import functional_call, hessian
from time import time


def train_model(gene_data: GeneData,
                dataset_metadata: DatasetMetadata,
                pol_2_total_loss: Pol2TotalLoss,
                device: str,
                max_epochs=100,
                max_patience=5,
                loss_change_tolerance=1e-7
                ) -> tuple[Pol2Model, float]:
    model = Pol2Model(feature_names=feature_names, intron_names=gene_data.intron_names).to(device)
    model.initialize_intercepts(gene_data, library_sizes)

    optimizer = optim.LBFGS(model.parameters(),
                            lr=1.0,
                            max_iter=50,
                            tolerance_change=1e-09,
                            tolerance_grad=1e-07,
                            history_size=100,
                            line_search_fn=None)

    def closure():
        optimizer.zero_grad()
        predicted_log_reads_exon, predicted_reads_intron, phi = model(dataset_metadata.design_matrix,
                                                                      dataset_metadata.log_library_sizes)
        loss = pol_2_total_loss(reads_exon=gene_data.exon_reads,
                                reads_introns=gene_data.intron_reads,
                                coverage=gene_data.coverage,
                                predicted_log_reads_exon=predicted_log_reads_exon,
                                predicted_reads_intron=predicted_reads_intron,
                                phi=phi)
        loss.backward()
        return loss

    previous_loss = None
    patience_counter = 0

    for epoch in trange(max_epochs):
        loss = optimizer.step(closure).item()

        if previous_loss is not None:
            relative_loss_change = abs(loss - previous_loss) / (abs(previous_loss) + 1e-10)
            if relative_loss_change < loss_change_tolerance:
                patience_counter += 1
                if patience_counter >= max_patience:
                    break
            else:
                patience_counter = 0

        previous_loss = loss

    return model, loss


# def get_param_names_and_sizes(model_parameters: dict[str, nn.Parameter]
#                               ) -> tuple[list[str], list[int]]:
#     parameter_names: list[str] = []
#     parameter_sizes: list[int] = []
#     for name, parameter in model_parameters.items():
#         parameter_names.append(name)
#         parameter_sizes.append(parameter.numel())
#     return parameter_names, parameter_sizes


def flatten_hessian_dict(hessian_dict, model_parameters):
    param_names = list(model_parameters.keys())
    param_sizes = {name: model_parameters[name].numel() for name in param_names}
    param_offsets = {name: sum(param_sizes[n] for n in param_names[:i])
                     for i, name in enumerate(param_names)}
    total_size = sum(param_sizes.values())
    H = torch.zeros((total_size, total_size))

    for i_name in param_names:
        for j_name in param_names:
            block = hessian_dict[i_name][j_name].detach()
            block = block.reshape(param_sizes[i_name], param_sizes[j_name])
            H[param_offsets[i_name]:param_offsets[i_name] + param_sizes[i_name],
            param_offsets[j_name]:param_offsets[j_name] + param_sizes[j_name]] = block

    return H


if __name__ == "__main__":
    # Load training data
    input_folder = Path('/home/jakub/Desktop/dev-pol-ii-analysis/data/drosophila_mutants/refactored_train_data/')
    with open(input_folder / 'train_data.pkl', 'rb') as file:
        train_data = pickle.load(file)

    gene_data_dict = train_data[0]
    gene_data = GeneData(**gene_data_dict)

    # Load or create dataset metadata    
    design_matrix_df = pd.DataFrame(data={'genotype_rp2': 6 * [0] + 6 * [1],
                                          'dummy': 3 * [1] + 6 * [0] + 3 * [1]})
    # design_matrix_df = pd.DataFrame(data={'genotype_rp2': 6 * [0] + 6 * [1]})    

    feature_names = design_matrix_df.columns.to_list()
    design_matrix = torch.tensor(design_matrix_df.to_numpy()).float()

    library_sizes = torch.ones(len(design_matrix))

    dataset_metadata = DatasetMetadata(design_matrix=design_matrix,
                                       library_sizes=library_sizes,
                                       log_library_sizes=torch.log(library_sizes),
                                       feature_names=feature_names)

    # both gene_data and dataset_metadata are loaded now
    device = 'cpu'
    gene_data = gene_data.to(device)
    dataset_metadata = dataset_metadata.to(device)

    pol_2_total_loss = Pol2TotalLoss().to(device)

    t0 = time()
    model, loss = train_model(gene_data=gene_data,
                              dataset_metadata=dataset_metadata,
                              pol_2_total_loss=pol_2_total_loss,
                              device=device)
    t1 = time()
    time_train = 1000 * (t1 - t0)
    print(f"{time_train=}")

    model_parameters = dict(model.named_parameters())
    # parameter_names, parameter_sizes = get_param_names_and_sizes(model_parameters)

    t2 = time()


    def loss_by_model_parameters(model_parameters):
        outputs = functional_call(model, model_parameters, (dataset_metadata.design_matrix,
                                                            dataset_metadata.log_library_sizes))
        predicted_log_reads_exon, predicted_reads_intron, phi = outputs
        return pol_2_total_loss(reads_exon=gene_data.exon_reads,
                                reads_introns=gene_data.intron_reads,
                                coverage=gene_data.coverage,
                                predicted_log_reads_exon=predicted_log_reads_exon,
                                predicted_reads_intron=predicted_reads_intron,
                                phi=phi)


    hessian_dict = hessian(loss_by_model_parameters)(model_parameters)
    t3 = time()
    time_hessian_computation = 1000 * (t3 - t2)
    print(f"{time_hessian_computation=}")

    H = flatten_hessian_dict(hessian_dict, model_parameters)

    # parameter_data: list[dict] = []
    # for name, size in zip(parameter_names, parameter_sizes):
    #     if '.' in name:
    #         parameter_type = name.split('.', 1)[0]
    #         intron_name = name.split('.', 1)[1]
    #     else:
    #         parameter_type = name.split('.', 1)[0]
    #         intron_name = None

    #     parameter_value_list = [model_parameters[name].item()] if model_parameters[name].numel() == 1 \
    #         else model_parameters[name].tolist()
    #     for feature_index in range(size):
    #         feature_name = None
    #         if parameter_type in ('alpha', 'beta', 'gamma'):
    #             feature_name = feature_names[feature_index]
    #         parameter_data.append({'parameter_type': parameter_type,
    #                                'intron_name': intron_name,
    #                                'feature_name': feature_name,
    #                                'value': parameter_value_list[feature_index]})

    # df_param = pd.DataFrame(data=parameter_data)

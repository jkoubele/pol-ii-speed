import numpy as np
import pandas as pd
import cvxpy as cp
import torch
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.special import expit
import pickle
from dataclasses import dataclass
from queue import PriorityQueue

from cvxpy.constraints import Constraint

from pol_ii_speed_modeling.pol_ii_model import (
    GeneData, DatasetMetadata, Pol2Model, Pol2TotalLoss
)


from pol_ii_speed_modeling.load_dataset import load_dataset_from_results_folder

with open("gene_data_list.pkl", "rb") as file:
    gene_data_list = pickle.load(file)
    
with open("dataset_metadata.pkl", "rb") as file:
    dataset_metadata = pickle.load(file)
    
# dataset_metadata, gene_data_list = load_dataset_from_results_folder(
#     results_folder=Path('/cellfile/datapublic/jkoubele/data_pol_ii/mouse_myocardium/results/'),
#     log_output_folder=Path('/cellfile/datapublic/jkoubele/data_pol_ii/mouse_myocardium/results/logs'),
#     gene_names_file_name='test.csv')
#%%
POISSON_LOSS_EPSILON = 1e-12
MAX_GAP_TOLERANCE = 0.25

gene_data = gene_data_list[0]
design_matrix = dataset_metadata.design_matrix.numpy()

n, p = design_matrix.shape
constraints = []

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

#%% 1 intron only, ignore exon term


# node_constraints = []

# example_row_0 = np.array([0.0], dtype=np.float32)
# example_row_1 = np.array([1.0], dtype=np.float32)

# node_constraints += [intercept_intron + theta - example_row_1 @ beta <= 1.0]
# node_constraints += [intercept_intron + example_row_1 @ gamma <= 1.0]
# # node_constraints = []



# loss = 0
# problem_constraints = []
# problem_variables = {'u': [], 'a': [], 'b': []}


# for i, row in enumerate(design_matrix):
#     a = intercept_intron + theta - row @ beta
#     b = intercept_intron + row @ gamma
#     # construct LP to get bounds on A and B   
#     min_a = cp.Problem(cp.Minimize(a), node_constraints).solve()
#     max_a = cp.Problem(cp.Maximize(a), node_constraints).solve()
    
#     min_b = cp.Problem(cp.Minimize(b), node_constraints).solve()
#     max_b = cp.Problem(cp.Maximize(b), node_constraints).solve()
    
#     # These min/max may not be tight for more complex constrains --> optimize v directly 
#     # in a convex program instead of the LPs above
#     min_v = np.exp(min_a) + np.exp(min_b)
#     max_v = np.exp(max_a) + np.exp(max_b)
    
#     n = intron_reads[i]
    
#     if max_v <= n:
#         loss += max_v - n * np.log(max_v + POISSON_LOSS_EPSILON)
#         print(f"{max_v} <= {n}")
#     else:
#         u = cp.Variable()
#         problem_constraints.append(u >= cp.exp(a) + cp.exp(b))
        
#         problem_variables['u'] += [u]
#         problem_variables['a'] += [a]
#         problem_variables['b'] += [b]
#         loss += u - n * cp.log(u + POISSON_LOSS_EPSILON)
#         pass
    
    
# problem = cp.Problem(cp.Minimize(loss), node_constraints + problem_constraints)
# node_loss = problem.solve() # replace by a certified lower bound (dual objective) later on
# print(f"{node_loss=}")

#%% BnB

beta = cp.Variable(p) # shape (p,)
gamma = cp.Variable(p) # shape (p,)
intercept_intron = cp.Variable()
theta = cp.Variable()

logit_transcribing = [intercept_intron + theta - row @ beta for row in design_matrix]
logit_unspliced = [intercept_intron + row @ gamma for row in design_matrix]


best_feasible_loss = sum([intron_reads.mean() - n * np.log(intron_reads.mean()) for n in intron_reads])

def poisson_loss(prediction: np.array, target: np.array) -> np.array:
    return np.array([pred - target_value * np.log(pred + POISSON_LOSS_EPSILON)
                     for pred, target_value in zip(prediction, target)])
    
@dataclass
class Node:
    parent_lower_bound: float
    node_constraints: list[Constraint]
    
    def __lt__(self, other):
        return self.parent_lower_bound < other.parent_lower_bound
    
root_node = Node(parent_lower_bound=-np.inf,
                 node_constraints=[])
stack = PriorityQueue()
stack.put(root_node)

gaps = []
i = 0
while not stack.empty():
    node = stack.get()
    print(f"{stack.qsize()=}")
    print(50*'--')
    
    i += 1    
    
    loss = cp.Constant(0.0)
    problem_constraints: list[Constraint] = []
    lifted_variables = [cp.Variable() for _ in range(len(design_matrix))]
    
    node_is_infeasible = False
    for a, b, n, u in zip(logit_transcribing, logit_unspliced, intron_reads, lifted_variables):
        
        # Solve small LPs to get bounds on logits. Optionally solve full convex program
        # for predicted intron reads directly, which may be tighter that LPs for more complex constraints.
        
        problem_min_a = cp.Problem(cp.Minimize(a), node.node_constraints)
        min_a = problem_min_a.solve(solver=cp.GLPK)
          
        
        problem_max_a = cp.Problem(cp.Maximize(a), node.node_constraints)
        max_a = problem_max_a.solve(solver=cp.GLPK)
        
        
        problem_min_b = cp.Problem(cp.Minimize(b), node.node_constraints)
        min_b = problem_min_b.solve(solver=cp.GLPK)
             
        
        problem_max_b = cp.Problem(cp.Maximize(b), node.node_constraints)
        max_b = problem_max_b.solve(solver=cp.GLPK)
        
        
        if any([problem_min_a.status == cp.INFEASIBLE, 
                problem_max_a.status == cp.INFEASIBLE,
                problem_min_b.status == cp.INFEASIBLE,
                problem_max_b.status == cp.INFEASIBLE]):
            node_is_infeasible = True
            break        
        
        # min_b = cp.Problem(cp.Minimize(b), node.node_constraints).solve()
        # max_b = cp.Problem(cp.Maximize(b), node.node_constraints).solve()        
        
        min_v = np.exp(min_a) + np.exp(min_b)
        max_v = np.exp(max_a) + np.exp(max_b)  
        
        if max_v <= 0:
            print(f"{max_v}")
            assert False
        
        loss += u - n * cp.log(u + POISSON_LOSS_EPSILON)
        
        if max_v <= n:
            print(f"{max_v=}")
            problem_constraints += [u <= max_v]
        else:
            problem_constraints += [u >= cp.exp(a) + cp.exp(b)]   
            
    print(f"{node_is_infeasible=}")
    if node_is_infeasible:
        continue    
            
    problem = cp.Problem(cp.Minimize(loss), node.node_constraints + problem_constraints)    
    loss_relaxed = problem.solve() # replace by a certified lower bound (dual objective) later on
    print(f"{loss_relaxed=}")
    print(f"{problem.status=}")
    
    if problem.status == cp.INFEASIBLE:
        continue
    
    if loss_relaxed >= best_feasible_loss:
        continue
    
    # Compute node feasible loss
    A = intercept_intron + theta - design_matrix @ beta
    B = intercept_intron + design_matrix @ gamma
    V = cp.exp(A) + cp.exp(B)
    print(f"{V.value=}")
    node_feasible_loss= sum(poisson_loss(V.value, intron_reads))
    
    if node_feasible_loss + MAX_GAP_TOLERANCE< best_feasible_loss:
        print(f"{best_feasible_loss=}")
        best_feasible_loss = loss_relaxed
        print(f"{node_feasible_loss=}")
        # Store also the best parameters
        
    if loss_relaxed + MAX_GAP_TOLERANCE < best_feasible_loss:
        sample_loss_relaxed = poisson_loss(np.array([u.value for u in lifted_variables]), intron_reads)
        sample_loss_original = poisson_loss(V.value, intron_reads)
        sample_gap =  sample_loss_original - sample_loss_relaxed
        
        print(f"{sample_gap.max()}")
        gaps.append(sample_gap.max())
        
        largest_gap_sample = int(sample_gap.argmax())
        x_for_bound = design_matrix[largest_gap_sample]
        bound_reads = (V.value[largest_gap_sample] + intron_reads[largest_gap_sample]) / 2
        print(f"{bound_reads=}")
        # Expand the current node
        cut_boundary_1 = intercept_intron + theta - x_for_bound @ beta - np.log(bound_reads / 2)
        cut_boundary_2 = intercept_intron + x_for_bound @ gamma -  np.log(bound_reads / 2)
        
        cut_1 = [intercept_intron + theta - x_for_bound @ beta <=  np.log(bound_reads / 2),
                 intercept_intron + x_for_bound @ gamma <=  np.log(bound_reads / 2)]
        
        cut_2 = [intercept_intron + theta - x_for_bound @ beta <=  np.log(bound_reads / 2),
                 intercept_intron + x_for_bound @ gamma >=  np.log(bound_reads / 2)]
        
        cut_3 = [intercept_intron + theta - x_for_bound @ beta >=  np.log(bound_reads / 2),
                 intercept_intron + x_for_bound @ gamma <=  np.log(bound_reads / 2)]
        
        cut_4 = [intercept_intron + theta - x_for_bound @ beta >=  np.log(bound_reads / 2),
                 intercept_intron + x_for_bound @ gamma >=  np.log(bound_reads / 2)]
        
        for cut in [cut_1, cut_2,cut_3,cut_4]:
            stack.put(Node(parent_lower_bound=loss_relaxed,
                             node_constraints=node.node_constraints + cut))
        
        
        # for inequality_signs in [(1,1), (1,-1), (-1,1), (-1,-1)]:
        #     new_constraints = [cut_boundary_1 * inequality_signs[0] <= 0,
        #                        cut_boundary_2 * inequality_signs[1] <= 0]
        #     child = Node(parent_lower_bound=loss_relaxed,
        #                  node_constraints=node.node_constraints + new_constraints)
        #     stack.put(child)
        #     print(f"{child.node_constraints=}")
    
    
    
        
        


# A = intercept_intron + theta - design_matrix @ beta
# B = intercept_intron + design_matrix @ gamma


#%%
# predicted_reads_intron = cp.exp(transcribing_introns_term) + cp.exp(unspliced_introns_term)
# intron_loss_original = cp.sum(predicted_reads_intron - cp.multiply(intron_reads, cp.log(predicted_reads_intron)))


# t_vec = np.full(n, 0.5)
# transcribing_introns_surrogate = cp.Variable(n, nonneg=True)
# unspliced_introns_surrogate = cp.Variable(n, nonneg=True)
# predicted_reads_intron_surrogate = cp.Variable(n, nonneg=True)

# constraints += [
#     transcribing_introns_surrogate >= cp.exp(transcribing_introns_term),
#     unspliced_introns_surrogate >= cp.exp(unspliced_introns_term),
#     predicted_reads_intron_surrogate == transcribing_introns_surrogate + unspliced_introns_surrogate,
# ]



# lse_bound = cp.log(predicted_reads_intron_surrogate)

# intron_loss_relaxed = cp.sum(predicted_reads_intron_surrogate - cp.multiply(intron_reads, lse_bound))

# custom_value = cp.sum(2 * alpha) 



# objective = cp.Minimize(exon_loss + intron_loss_relaxed)

# problem = cp.Problem(objective, constraints)
# ret = problem.solve()

# ret_stats = problem.solver_stats

# print(f"{problem.status=}")

# reads_1 = exon_reads[(design_matrix==0).flatten()]
# reads_2 = exon_reads[(design_matrix==1).flatten()]


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
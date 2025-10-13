import numpy as np
import pandas as pd
import cvxpy as cp
import torch
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.special import expit
import pickle
from dataclasses import dataclass, field
from queue import PriorityQueue
from typing import Optional, NamedTuple
from scipy.spatial import ConvexHull

from cvxpy.constraints import Constraint

from pol_ii_speed_modeling.pol_ii_model import (
    GeneData, DatasetMetadata, Pol2Model, Pol2TotalLoss
)

from pol_ii_speed_modeling.load_dataset import load_dataset_from_results_folder

with open("gene_data_list.pkl", "rb") as file:
    gene_data_list = pickle.load(file)

with open("dataset_metadata.pkl", "rb") as file:
    dataset_metadata = pickle.load(file)



# %%
MAX_GAP_TOLERANCE = 0.01
num_position_coverage = 100
locations = np.linspace(1 / (2 * num_position_coverage), 1 - 1 / (2 * num_position_coverage), num_position_coverage)


def get_loss_by_pi(pi: float, coverage: np.array) -> float:
    loss_by_location = -np.log(1 + pi - 2 * pi * locations)
    return loss_by_location @ coverage


def get_loss_by_logit(logit: float, coverage: np.array) -> float:
    return get_loss_by_pi(pi=expit(logit), coverage=coverage)

gene_data = gene_data_list[1]
design_matrix = dataset_metadata.design_matrix.numpy()
coverage = gene_data.coverage[2, 3, :].numpy()

loss_left_asymp = get_loss_by_pi(1e-10, coverage)
loss_right_asymp = get_loss_by_pi(1 - 1e-10, coverage)

num_logits_sampled = 100
logit_max_abs_value = 5

# TO-DO: Choose logits_range such that logits_range[0] ≈ loss_left_asymp and logits_range[-1] ≈ loss_right_asymp

logits_range = np.linspace(-logit_max_abs_value,
                           logit_max_abs_value,
                           num_logits_sampled)

loss_by_logits = np.array([get_loss_by_logit(logit, coverage) for logit in logits_range])


plt.plot(logits_range, loss_by_logits)
plt.title("Loss by logit (initial plot)")
plt.xlabel("Logit")
plt.ylabel("Loss")
plt.show()



#%%


class LogitBounds(NamedTuple):
    min_logit: float
    max_logit: float      
    
def create_cvx_envelope(xs: np.ndarray, ys: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    pts = np.column_stack([xs, ys])
    hull = ConvexHull(pts)
    ring = hull.vertices  # CCW order around the hull
    ring_pts = pts[ring]

    # indices of leftmost and rightmost vertices on the ring
    i_left = np.argmin(ring_pts[:, 0])
    i_right = np.argmax(ring_pts[:, 0])

    if i_left <= i_right:
        lower_chain = ring_pts[i_left:i_right + 1]
    else:
        # wrap-around slice
        lower_chain = np.vstack([ring_pts[i_left:], ring_pts[:i_right + 1]])

    # ensure strictly increasing x and remove any duplicates/flat collinearities
    order = np.argsort(lower_chain[:, 0])
    x_chain = lower_chain[order, 0]
    y_chain = lower_chain[order, 1]
    return x_chain, y_chain
    

class GlobalBoundAndCvxEnvelope(NamedTuple):
    global_lb: float
    create_envelope: bool  
    envelope_x: Optional[np.ndarray]
    envelope_y: Optional[np.ndarray]
    
def get_global_lb_and_cvx_envelope(logit_lb: float,
                                   logit_ub: float,
                                   coverage: np.ndarray,
                                   logits_range: np.ndarray,
                                   loss_by_logits: np.ndarray) -> GlobalBoundAndCvxEnvelope:
    
    loss_by_logits_argmin = loss_by_logits.argmin()
    logit_optimal = logits_range[loss_by_logits_argmin]
    loss_min = loss_by_logits[loss_by_logits_argmin]
    
    create_envelope = False
    envelope_lb: Optional[float] = None
    envelope_ub: Optional[float] = None
    global_lb = loss_min

    if logit_lb == -np.inf:
        if logit_ub <= logit_optimal:
            # Only global bound of logit_ub
            global_lb = get_loss_by_logit(logit_ub, coverage)
            
        elif logit_optimal < logit_ub < np.inf:
            # Global bound of logit_optimal + envelope on (logit_optimal, logit_ub)        
            create_envelope = True
            envelope_lb = logit_optimal
            envelope_ub = logit_ub        
            
        elif logit_ub == np.inf:
            pass # Only global bound of logit_optimal            
            
        else:
            raise RuntimeError()
    
    
    elif -np.inf < logit_lb:
        if logit_ub < np.inf:
            # Envelope between (logit_lb, logit_ub). Global bound as safeguard.
            create_envelope = True
            envelope_lb = logit_lb
            envelope_ub = logit_ub    
            
        elif logit_ub == np.inf:
            if logit_lb < logit_optimal:
                # Envelope on (logit_lb, logit_optimal) + global bound of logit_optimal
                create_envelope = True
                envelope_lb = logit_lb
                envelope_ub = logit_optimal    
                
                
            elif logit_lb >= logit_optimal:
                # Global bound of logit_lb
                global_lb = get_loss_by_logit(logit_lb, coverage)
                
            else:
                raise RuntimeError()
        else:
            raise RuntimeError()
    
    else:
        raise RuntimeError()
    
    
    if create_envelope:
        # TO-DO: handle corner cases
        envelope_lb_logit_index = logits_range.searchsorted(envelope_lb, side='left')
        envelope_ub_logit_index = max(0, logits_range.searchsorted(envelope_ub, side='left') - 1)
        
        if envelope_ub_logit_index - envelope_lb_logit_index  < 1:
            # envelope has less than 2 points, replace it by a global bound
            global_lb = min(get_loss_by_logit(envelope_lb, coverage), 
                            get_loss_by_logit(envelope_ub, coverage))   
            create_envelope = False
        else:        
            xs = logits_range[envelope_lb_logit_index:envelope_ub_logit_index]
            ys = loss_by_logits[envelope_lb_logit_index:envelope_ub_logit_index]
            envelope_x, envelope_y = create_cvx_envelope(xs, ys)
            # TO-DO: move envelope down by safety margin (calculated from curvature bound)
            
    global_bound_and_cvx_envelope =  GlobalBoundAndCvxEnvelope(global_lb=global_lb,
                                     create_envelope=create_envelope,
                                     envelope_x=envelope_x if create_envelope else None,
                                     envelope_y=envelope_y if create_envelope else None)
            
    return global_bound_and_cvx_envelope
            
theta = cp.Variable(name='theta')

logit = 1 * theta # Will replace for formula with LFC and design matrix

@dataclass
class Node:
    parent_lower_bound: float
    node_constraints: list[Constraint]
    constraints_string: Optional[list[str]] = field(default_factory=list)
    
    def __lt__(self, other):
        return self.parent_lower_bound < other.parent_lower_bound
    
root_node = Node(parent_lower_bound=-np.inf,
                 node_constraints=[])
stack = PriorityQueue()
stack.put(root_node)
    

# Initiall (heuristic) solution to get feasible loss
best_theta  = 0.0
theta.value = best_theta # Replace by some dict of all params
best_feasible_loss = get_loss_by_logit(logit = logit.value, coverage=coverage)


while not stack.empty():
    node = stack.get()
    
    print(100*'-')
    print(node)
    
    if best_feasible_loss <= node.parent_lower_bound:
        print("Prunning", f"{best_feasible_loss}", f"{node.parent_lower_bound}")
        continue
    
    logit_bounds_list: list[LogitBounds] = []
    
    node_is_infeasible = False        
    
    
    
    # Iterate over samples/strata
    problem_min_logit = cp.Problem(cp.Minimize(logit), node.node_constraints)
    min_logit = problem_min_logit.solve(solver=cp.GLPK)   

    problem_max_logit = cp.Problem(cp.Maximize(logit), node.node_constraints)
    max_logit = problem_max_logit.solve(solver=cp.GLPK)     
    
    if problem_min_logit.status == cp.INFEASIBLE or problem_max_logit.status == cp.INFEASIBLE:          
        node_is_infeasible = True
        
    if node_is_infeasible:
        raise NotImplementedError() # handle infeasible nodes - continue the while loop
        
    logit_bounds = LogitBounds(min_logit=min_logit,
                               max_logit=max_logit)
    
    global_bound_and_cvx_envelope = get_global_lb_and_cvx_envelope(logit_lb=logit_bounds.min_logit,
                                                                   logit_ub=logit_bounds.max_logit,
                                                                   coverage=coverage,
                                                                   logits_range=logits_range,
                                                                   loss_by_logits=loss_by_logits)
    
    # Construct the convex program and solve it    
    t = cp.Variable()

    program_constraints: list[Constraint] = [t >= global_bound_and_cvx_envelope.global_lb]
    
    print(f"{global_bound_and_cvx_envelope=}")

    if global_bound_and_cvx_envelope.create_envelope:
        envelope_x = global_bound_and_cvx_envelope.envelope_x
        envelope_y = global_bound_and_cvx_envelope.envelope_y        
        
        
        for i in range(len(envelope_x)-1):
            x_1 = envelope_x[i]
            x_2 = envelope_x[i+1]
            y_1 = envelope_y[i]
            y_2 = envelope_y[i+1]
            program_constraints += [t >= y_1 + (logit - x_1) * (y_2 - y_1) / (x_2 - x_1)]

    prob = cp.Problem(cp.Minimize(t), program_constraints + node.node_constraints)
    relaxed_loss = prob.solve()
    
    if theta.value is None:
        theta.value = 0.0

    original_loss = get_loss_by_logit(logit.value, coverage)
    
    plt.plot(logits_range, loss_by_logits, label='Original loss')
    if global_bound_and_cvx_envelope.create_envelope:
        plt.plot(global_bound_and_cvx_envelope.envelope_x, 
                 global_bound_and_cvx_envelope.envelope_y,
                 label='Convex envelope')    
    plt.plot(logits_range, global_bound_and_cvx_envelope.global_lb * np.ones_like(logits_range), 
             label='Global lower bound')
    plt.title("Loss by logit")
    plt.xlabel("Logit")
    plt.ylabel("Loss")
    plt.legend()
    plt.show()
    
    print(f"{theta.value=}")
    print(f"{relaxed_loss=}")
    print(f"{original_loss=}")
    
    if relaxed_loss >= best_feasible_loss:
        print("Prunning the node (not expanding it)")
        continue # Prunning the node (not expanding it)
    
    if original_loss < best_feasible_loss:
        print("Setting the best new feasible loss")
        best_feasible_loss = original_loss
        best_theta = theta.value
        
    if relaxed_loss + MAX_GAP_TOLERANCE < best_feasible_loss:
        print(f'Expanding at theta={float(theta.value)} ')
        # Heuristic to expand the node
        cut_1 = [theta <= float(theta.value)]
        cut_2 = [theta >= float(theta.value)]
        
        cut_1_string = [f"theta <= {float(theta.value)}"]
        cut_2_string = [f"theta >= {float(theta.value)}"]
        
        for cut, cut_string in zip((cut_1, cut_2), (cut_1_string, cut_2_string)):
            child = Node(parent_lower_bound=relaxed_loss,
                       node_constraints=node.node_constraints + cut,
                       constraints_string=node.constraints_string + cut_string) 
            
            stack.put(child)
    
        
    
           
    

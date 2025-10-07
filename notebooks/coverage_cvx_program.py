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
    
    
def sigmoid_second_derivative(s: float) -> float:
    return (1-2*s)*(1-s)*s

#%%
num_position_coverage = 100
locations = np.linspace(1/(2*num_position_coverage), 1-1/(2*num_position_coverage), num_position_coverage)

def loss_by_pi(pi: float, coverage: np.array) -> float:
    loss_by_location = -np.log(1 + pi - 2*pi*locations)
    return loss_by_location @ coverage

# def lower_convex_envelope_via_hull(xs, ys):
    

    # Optional: clean with the same monotone-chain pruning to remove collinear points
    # xe, ye = lower_convex_envelope(x_chain, y_chain)
    # return xe, ye

gene_data = gene_data_list[1]
design_matrix = dataset_metadata.design_matrix.numpy()
coverage = gene_data.coverage[2,3,:].numpy()

num_logits_sampled = 200
logit_max_abs_value = 3
logits_range = np.linspace(-logit_max_abs_value,
                          logit_max_abs_value,
                          num_logits_sampled)

loss_by_logits = np.array([loss_by_pi(pi, coverage) for pi in expit(logits_range)])

loss_left_asymp = loss_by_pi(1e-10, coverage)
loss_right_asymp = loss_by_pi(1-1e-10, coverage)

plt.plot(logits_range, loss_by_logits)
plt.title("Loss by logit")
plt.xlabel("Logit")
plt.ylabel("Loss")
plt.show()


# points = np.stack((logits_range, loss_by_logits), axis=1)
# hull = ConvexHull(points)


# Obtain from LP
logit_lb = -2
logit_ub = 1

if logit_lb == -np.inf and logit_ub == np.inf:
    # Add just a global lower bound on the loss
    pass

else:
    envelope_lb = logit_lb
    envelope_ub =  logit_ub
    
    if logit_lb == -np.inf:        
        pass
    





xs = logits_range
ys = loss_by_logits

def get_convex_envelope(xs: np.ndarray, xy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    pts = np.column_stack([xs, ys])
    hull = ConvexHull(pts)
    ring = hull.vertices  # CCW order around the hull
    ring_pts = pts[ring]
    
    # indices of leftmost and rightmost vertices on the ring
    i_left = np.argmin(ring_pts[:,0])
    i_right = np.argmax(ring_pts[:,0])
    
    if i_left <= i_right:
        lower_chain = ring_pts[i_left:i_right+1]
    else:
        # wrap-around slice
        lower_chain = np.vstack([ring_pts[i_left:], ring_pts[:i_right+1]])
    
    # ensure strictly increasing x and remove any duplicates/flat collinearities
    order = np.argsort(lower_chain[:,0])
    x_chain = lower_chain[order,0]
    y_chain = lower_chain[order,1]
    return x_chain, y_chain

x_chain, y_chain = get_convex_envelope(logits_range, loss_by_logits)
plt.plot(logits_range, loss_by_logits, label='Original loss')
plt.plot(x_chain, y_chain, label='Convex envelope')
plt.title("Loss by logit")
plt.xlabel("Logit")
plt.ylabel("Loss")
plt.legend()
plt.show()




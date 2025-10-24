import pickle
from dataclasses import dataclass, field
from queue import PriorityQueue
from typing import Optional, NamedTuple

import cvxpy as cp
import numpy as np
from cvxpy.constraints import Constraint

from coverage_cvx_envelope import IntronCoverage, get_loss_by_logit, fit_coverage

MAX_GAP_TOLERANCE = 0.3


class IntronAndSample(NamedTuple):
    intron_name: str
    sample_id: int


class CoverageData(NamedTuple):
    intron_names: list[str]
    feature_vectors: dict[int, np.ndarray]
    intron_coverages: dict[IntronAndSample, IntronCoverage]


def load_coverage_data(collapse_samples_by_strata=True) -> CoverageData:
    with open("gene_data_list.pkl", "rb") as file:
        gene_data_list = pickle.load(file)

    with open("dataset_metadata.pkl", "rb") as file:
        dataset_metadata = pickle.load(file)

    gene_data = gene_data_list[1]
    design_matrix = dataset_metadata.design_matrix.numpy()
    intron_names = ['ENSMUSG00000073131_3', 'ENSMUSG00000073131_4', 'ENSMUSG00000073131_5', 'ENSMUSG00000073131_6']

    intron_coverages: dict[IntronAndSample, IntronCoverage] = {}

    if collapse_samples_by_strata:
        design_matrix_as_tuple_list = [tuple(row) for row in dataset_metadata.design_matrix.numpy()]
        strata_tuple_to_mask = {strata_tuple: [row == strata_tuple for row in design_matrix_as_tuple_list] for
                                strata_tuple in set(design_matrix_as_tuple_list)}

        reduced_design_matrix = np.asarray(list(strata_tuple_to_mask.keys()))

        num_introns_and_stratum = len(intron_names) * len(reduced_design_matrix)
        print(f"{num_introns_and_stratum=}")

        max_envelope_gap = 0.25 * MAX_GAP_TOLERANCE / num_introns_and_stratum

        for intron_name in intron_names:
            for strata_id, strata_row_as_tuple in enumerate(strata_tuple_to_mask):
                coverage_vector = gene_data.coverage.numpy()[strata_tuple_to_mask[strata_row_as_tuple],
                                  gene_data.intron_names.index(intron_name), :].sum(axis=0)
                intron_coverage = IntronCoverage(coverage=coverage_vector, max_envelope_gap=max_envelope_gap)
                intron_coverages[IntronAndSample(intron_name, strata_id)] = intron_coverage

        feature_vectors = {sample_id: row for sample_id, row in enumerate(reduced_design_matrix)}

    else:
        num_introns_and_samples = len(intron_names) * len(design_matrix)
        max_envelope_gap = 0.25 * MAX_GAP_TOLERANCE / num_introns_and_samples
        for intron_name in intron_names:
            for sample_id, _ in enumerate(design_matrix):
                coverage_vector = gene_data.coverage[sample_id, gene_data.intron_names.index(intron_name), :].numpy()
                intron_coverage = IntronCoverage(coverage=coverage_vector, max_envelope_gap=max_envelope_gap)
                intron_coverages[IntronAndSample(intron_name, sample_id)] = intron_coverage

        feature_vectors = {sample_id: row for sample_id, row in enumerate(design_matrix)}

    intron_coverages = {key: intron_coverage for key, intron_coverage in intron_coverages.items() if
                        intron_coverage.coverage.sum() > 0}
    coverage_data = CoverageData(
        intron_names=intron_names,
        feature_vectors=feature_vectors,
        intron_coverages=intron_coverages
    )
    return coverage_data


coverage_data = load_coverage_data()

coverage_group_0 = np.sum(
    [coverage_data.intron_coverages[IntronAndSample(coverage_data.intron_names[0], sample_id)].coverage for
     sample_id, feature_vector in
     coverage_data.feature_vectors.items() if feature_vector == 0], axis=0)

coverage_group_1 = np.sum(
    [coverage_data.intron_coverages[IntronAndSample(coverage_data.intron_names[0], sample_id)].coverage for
     sample_id, feature_vector in
     coverage_data.feature_vectors.items() if feature_vector == 1], axis=0)

optim_results_group_0 = fit_coverage(coverage_group_0)
optim_results_group_1 = fit_coverage(coverage_group_1)

optim_loss = optim_results_group_0.primal_objective + optim_results_group_1.primal_objective

# %%

thetas = {intron_name: cp.Variable(name=f"theta_{intron_name}") for intron_name in coverage_data.intron_names}
lfc_parameter = cp.Variable(coverage_data.feature_vectors[0].shape[0], name='LFC')

logits: dict[IntronAndSample, cp.Expression] = {}
logit_strings: dict[IntronAndSample, str] = {}

for intron_name in coverage_data.intron_names:
    for sample_id, feature_vector in coverage_data.feature_vectors.items():
        logits[IntronAndSample(intron_name, sample_id)] = thetas[intron_name] + lfc_parameter @ feature_vector
        logit_strings[IntronAndSample(intron_name,
                                      sample_id)] = f"{thetas[intron_name].name()} + {lfc_parameter.name()} @ {feature_vector}" if not np.all(
            feature_vector == 0) else f"{thetas[intron_name].name()}"

    # Thetas fitted to aggregated coverages for each intron
default_thetas: dict[str, float] = {}
for intron_name in coverage_data.intron_names:
    aggregated_intron_coverage = np.sum(
        [coverage_data.intron_coverages[IntronAndSample(intron_name, sample_id)].coverage for sample_id in
         coverage_data.feature_vectors.keys()], axis=0)
    optim_results = fit_coverage(aggregated_intron_coverage)
    default_thetas[intron_name] = np.clip(optim_results.logit_optimal, -5.0, 5.0).item()

default_lfc_parameter = np.zeros(lfc_parameter.shape)

best_thetas = dict(default_thetas)
best_lfc_parameter = default_lfc_parameter.copy()

best_feasible_loss = 0
for key, intron_coverage in coverage_data.intron_coverages.items():
    logit = best_thetas[key.intron_name] + np.dot(coverage_data.feature_vectors[key.sample_id], best_lfc_parameter)
    best_feasible_loss += get_loss_by_logit(logit, intron_coverage.coverage)


# %%


class LogitBounds(NamedTuple):
    min_logit: float
    max_logit: float


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

num_nodes_searched = 0
while not stack.empty():
    node = stack.get()
    num_nodes_searched += 1

    if best_feasible_loss <= node.parent_lower_bound:
        continue

    logit_bounds_by_intron_and_sample: dict[IntronAndSample, LogitBounds] = {}
    node_is_infeasible = False

    for key, logit in logits.items():
        problem_min_logit = cp.Problem(cp.Minimize(logit), node.node_constraints)
        min_logit = problem_min_logit.solve(solver=cp.GLPK)
        if problem_min_logit.status == cp.INFEASIBLE:
            node_is_infeasible = True
            break

        problem_max_logit = cp.Problem(cp.Maximize(logit), node.node_constraints)
        max_logit = problem_max_logit.solve(solver=cp.GLPK)

        if problem_max_logit.status == cp.INFEASIBLE:
            node_is_infeasible = True
            break

        logit_bounds_by_intron_and_sample[key] = LogitBounds(min_logit, max_logit)

    if node_is_infeasible:
        continue

    # Get lower bounds and envelopes, use them to construct cvx program
    epigraph_variables: dict[IntronAndSample, cp.Variable] = {}
    program_constraints: list[Constraint] = []
    for intron_and_sample, intron_coverage in coverage_data.intron_coverages.items():
        logit_bounds = logit_bounds_by_intron_and_sample[intron_and_sample]
        epigraph_lower_bound, envelope_points = intron_coverage.get_bound_and_envelope(logit_bounds.min_logit,
                                                                                       logit_bounds.max_logit)
        epigraph_variable = cp.Variable()
        epigraph_variables[intron_and_sample] = epigraph_variable
        program_constraints += [epigraph_variable >= epigraph_lower_bound]
        if envelope_points is not None:
            logit = logits[intron_and_sample]
            for i in range(len(envelope_points.envelope_x) - 1):
                x_1 = envelope_points.envelope_x[i]
                x_2 = envelope_points.envelope_x[i + 1]
                y_1 = envelope_points.envelope_y[i]
                y_2 = envelope_points.envelope_y[i + 1]
                program_constraints += [epigraph_variable >= y_1 + (logit - x_1) * (y_2 - y_1) / (x_2 - x_1)]

    objective = sum(epigraph_variables.values())

    problem = cp.Problem(cp.Minimize(objective), program_constraints + node.node_constraints)
    relaxed_loss = problem.solve()

    if problem.status == cp.INFEASIBLE:
        continue

    if relaxed_loss >= best_feasible_loss:
        continue

    original_loss = 0.0
    largest_gap = -np.inf
    intron_and_sample_with_largest_gap = None
    for intron_and_sample, epigraph_variable in epigraph_variables.items():
        logit = logits[intron_and_sample]
        # For unconstrained logits, model parameters are not assigned, we assign them default values
        if thetas[intron_and_sample.intron_name].value is None:
            thetas[intron_and_sample.intron_name].value = default_thetas[intron_and_sample.intron_name]
        if lfc_parameter.value is None:
            lfc_parameter.value = default_lfc_parameter

        intron_coverage = coverage_data.intron_coverages[intron_and_sample]
        original_loss_in_intron_and_sample = get_loss_by_logit(logit.value, intron_coverage.coverage)
        original_loss += original_loss_in_intron_and_sample
        gap = original_loss_in_intron_and_sample - epigraph_variable.value
        print(f"{intron_and_sample} has gap {gap=}")
        if gap > largest_gap:
            largest_gap = gap
            intron_and_sample_with_largest_gap = intron_and_sample

    if original_loss < best_feasible_loss:
        print(50 * '*')
        print(f"Updating best feasible loss {best_feasible_loss} -> {original_loss}")
        print(50 * '*')
        best_feasible_loss = original_loss
        best_thetas = {intron_name: theta.value for intron_name, theta in thetas.items()}
        best_lfc_parameter = lfc_parameter.value

    if relaxed_loss + MAX_GAP_TOLERANCE < best_feasible_loss:
        print(f"Expanding node {node=}")
        cut_1 = logits[intron_and_sample_with_largest_gap] >= logits[intron_and_sample_with_largest_gap].value
        cut_2 = logits[intron_and_sample_with_largest_gap] <= logits[intron_and_sample_with_largest_gap].value

        cut_1_string = f"{logit_strings[intron_and_sample_with_largest_gap]} >= {logits[intron_and_sample_with_largest_gap].value}"
        cut_2_string = f"{logit_strings[intron_and_sample_with_largest_gap]} <= {logits[intron_and_sample_with_largest_gap].value}"

        child_1 = Node(parent_lower_bound=relaxed_loss,
                       node_constraints=node.node_constraints + [cut_1],
                       constraints_string=node.constraints_string + [cut_1_string])
        child_2 = Node(parent_lower_bound=relaxed_loss,
                       node_constraints=node.node_constraints + [cut_2],
                       constraints_string=node.constraints_string + [cut_2_string])

        stack.put(child_1)
        stack.put(child_2)

    if len(node.node_constraints) > 20:
        print("DEBUG BREAK")
        break

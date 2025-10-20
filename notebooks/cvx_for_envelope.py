import numpy as np
import cvxpy as cp
import pickle
from scipy.special import expit
from scipy.special import logit as logit_function

NUM_POSITION_COVERAGE = 100
LOCATIONS = np.linspace(1 / (2 * NUM_POSITION_COVERAGE),
                        1 - 1 / (2 * NUM_POSITION_COVERAGE),
                        NUM_POSITION_COVERAGE)

with open('coverage.pkl', 'rb') as f:
    coverage_mid = pickle.load(f)


def get_loss_by_pi(pi: float, coverage: np.array) -> float:
    loss_by_location = -np.log(1 + pi - 2 * pi * LOCATIONS)
    return loss_by_location @ coverage


def get_loss_by_logit(logit: float, coverage: np.array) -> float:
    return get_loss_by_pi(pi=expit(logit), coverage=coverage)


def get_loss_derivative_by_logit(logit: float, coverage: np.array) -> float:
    s = expit(logit)
    derivative_by_location = - (1 - 2 * LOCATIONS) * s * (1 - s) / (1 + s - 2 * s * LOCATIONS)
    return coverage @ derivative_by_location


coverage_left = np.zeros_like(coverage_mid)
coverage_left[0:10] = 3.0
coverage_uniform = np.ones_like(coverage_mid)

coverage = coverage_mid
pi = cp.Variable()
arg = 1 + pi - 2 * pi * LOCATIONS
loss = -cp.log(arg) @ coverage

pi_lower_bound = pi >= 0
pi_upper_bound = pi <= 1
problem = cp.Problem(cp.Minimize(loss), [pi_lower_bound, pi_upper_bound])
minimum_value = problem.solve(solver=cp.SCS)
dual_objective = problem.solver_stats.extra_stats['info']['dobj']
pi_value = float(pi.value)

tolerance_dual = 1e-6
tolerance_pi_value = 1e-4

lower_bound_active = False
upper_bound_active = False

if abs(pi_lower_bound.dual_value) > tolerance_dual:
    lower_bound_active = True
    if not pi_value <= tolerance_pi_value:
        raise RuntimeError("Lower bound active, but pi not close to 0 ({pi_value=})")

if abs(pi_upper_bound.dual_value) > tolerance_dual:
    upper_bound_active = True
    if not pi_value >= 1 - tolerance_pi_value:
        raise RuntimeError("Upper bound active, but pi not close to  1 ({pi_value=})")

print(f"{pi_value=}")
print(f"{lower_bound_active=}")
print(f"{pi_lower_bound.dual_value=}")
print(f"{upper_bound_active=}")
print(f"{pi_upper_bound.dual_value=}")


class IntronCoverage:

    def __init__(self, coverage: np.ndarray, locations: np.ndarray) -> None:
        pass

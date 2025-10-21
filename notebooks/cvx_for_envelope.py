import numpy as np
import cvxpy as cp
import pickle
import sympy as sp
from scipy.special import expit
from scipy.special import logit as logit_function
from typing import NamedTuple, Optional
import math

NUM_POSITION_COVERAGE = 100
LOCATIONS = np.linspace(1 / (2 * NUM_POSITION_COVERAGE),
                        1 - 1 / (2 * NUM_POSITION_COVERAGE),
                        NUM_POSITION_COVERAGE)

CURVATURE_BOUNDS = np.asarray([0.06493384, 0.06401682, 0.06308885, 0.06214971, 0.06119918,
                               0.06023702, 0.059263, 0.05827685, 0.05727834, 0.05626719,
                               0.05524313, 0.05420588, 0.05315516, 0.05209067, 0.05101209,
                               0.04991912, 0.04881142, 0.04768866, 0.04655048, 0.04539653,
                               0.04422643, 0.04303979, 0.04183621, 0.04061527, 0.03937656,
                               0.03811962, 0.03684399, 0.03554919, 0.03423472, 0.03290008,
                               0.03154472, 0.03016808, 0.0287696, 0.02734867, 0.02590466,
                               0.02443693, 0.02294479, 0.02142754, 0.01988444, 0.01831472,
                               0.01671759, 0.01509219, 0.01343766, 0.01175308, 0.01003749,
                               0.00828989, 0.00650921, 0.00469437, 0.00284419, 0.00095747,
                               0.00096709, 0.00293083, 0.00493516, 0.00698159, 0.00907169,
                               0.01120713, 0.01338968, 0.0156212, 0.01790367, 0.02023919,
                               0.02262999, 0.02507845, 0.02758708, 0.03015858, 0.03279583,
                               0.03550191, 0.03828011, 0.04113397, 0.0440673, 0.0470842,
                               0.05018908, 0.05338675, 0.05668238, 0.06008162, 0.0635906,
                               0.06721604, 0.07096528, 0.07484641, 0.07886833, 0.0830409,
                               0.08737508, 0.09188308, 0.0965786, 0.10147708, 0.10659599,
                               0.11195526, 0.11757777, 0.12348997, 0.12972272, 0.13631237,
                               0.14330225, 0.15074465, 0.15870368, 0.16725939, 0.176514,
                               0.18660168, 0.19770503, 0.21008465, 0.22413856, 0.24054176])

with open('coverage.pkl', 'rb') as f:
    coverage_mid = pickle.load(f)


def compute_curvature_bounds() -> np.ndarray:
    s, r = sp.symbols('s r', real=True)
    a = 1 - 2 * r

    curvature = a * s * (1 - s) * (a * s ** 2 + 2 * s - 1) / (1 + a * s) ** 2
    curvature_diff_wrt_s = sp.diff(curvature, s)
    diff_num, diff_den = sp.together(curvature_diff_wrt_s).as_numer_denom()

    # For s in (0,1) and r in (0,1), diff_den > 0 so we find critical points using diff_num only
    diff_num = sp.expand(diff_num)

    num_position_coverage = 100
    read_locations = np.linspace(1 / (2 * num_position_coverage), 1 - 1 / (2 * num_position_coverage),
                                 num_position_coverage)

    def get_curvature_bound(r_value: float) -> float:
        num_digits_precision = 50
        diff_num_substituted = sp.poly(diff_num.subs(r, r_value))
        root_expressions = sp.real_roots(diff_num_substituted)
        roots = [root.evalf(num_digits_precision) for root in root_expressions if 0 < root < 1]

        curvature_lambda = sp.lambdify(s, curvature.subs(r, r_value))
        curvature_bound = max([abs(float(curvature_lambda(root))) for root in roots], default=0.0)
        return curvature_bound

    curvature_bounds = np.asarray([get_curvature_bound(location) for location in read_locations])
    return curvature_bounds


def get_loss_by_pi(pi: float, coverage: np.array) -> float:
    loss_by_location = -np.log(1 + pi - 2 * pi * LOCATIONS)
    return loss_by_location @ coverage


def get_loss_by_logit(logit: float, coverage: np.array) -> float:
    return get_loss_by_pi(pi=expit(logit), coverage=coverage)


def get_loss_derivative_by_logit(logit: float, coverage: np.array) -> float:
    s = expit(logit)
    derivative_by_location = - (1 - 2 * LOCATIONS) * s * (1 - s) / (1 + s - 2 * s * LOCATIONS)
    return coverage @ derivative_by_location


class CoverageOptimResults(NamedTuple):
    primal_objective: float
    dual_objective: float
    pi_optimal: float
    logit_optimal: float
    pi_constraints_active: bool


def fit_coverage(coverage: np.ndarray) -> CoverageOptimResults:
    pi = cp.Variable()
    arg = 1 + pi - 2 * pi * LOCATIONS
    loss = -cp.log(arg) @ coverage

    pi_lower_bound = pi >= 0
    pi_upper_bound = pi <= 1
    problem = cp.Problem(cp.Minimize(loss), [pi_lower_bound, pi_upper_bound])
    minimum_value = problem.solve(solver=cp.SCS)
    dual_objective = problem.solver_stats.extra_stats['info']['dobj']
    pi_value = float(np.clip(pi.value, 0.0, 1.0))  # Clipping for values out of constraints, to prevent numeric issues

    tolerance_dual = 1e-6
    tolerance_primal = 1e-2

    lower_bound_active = False
    upper_bound_active = False

    if abs(pi_lower_bound.dual_value) > tolerance_dual and pi_value <= tolerance_primal:
        lower_bound_active = True
        pi_value = 0.0

    if abs(pi_upper_bound.dual_value) > tolerance_dual and pi_value >= 1 - tolerance_primal:
        upper_bound_active = True
        pi_value = 1.0

    optim_results = CoverageOptimResults(
        primal_objective=minimum_value,
        dual_objective=dual_objective,
        pi_optimal=pi_value,
        logit_optimal=logit_function(pi_value),
        pi_constraints_active=lower_bound_active or upper_bound_active)

    return optim_results


# curvature_bounds = compute_curvature_bounds()

# %%
coverage_left = np.zeros_like(coverage_mid)
coverage_left[0:10] = 3.0
coverage_uniform = np.ones_like(coverage_mid)

coverage = coverage_uniform
ret = fit_coverage(coverage / 100)
print(ret)


# %%

class EnvelopePoints(NamedTuple):
    envelope_x: np.ndarray
    envelope_y: np.ndarray


class IntronCoverage:

    def __init__(self, coverage: np.ndarray, max_envelope_gap=0.1) -> None:
        self.coverage = coverage
        self.optim_results = fit_coverage(coverage)

        loss_curvature_bound = coverage @ CURVATURE_BOUNDS

        self.max_sampling_distance = math.sqrt(max_envelope_gap / 8 / loss_curvature_bound)

    def get_bound_and_envelope(logit_lb: float, logit_ub: float) -> tuple[float, Optional[EnvelopePoints]]:
        # TODO
        pass


intron_coverage = IntronCoverage(coverage)

# %%

# coverage = coverage_mid
# pi = cp.Variable()
# arg = 1 + pi - 2 * pi * LOCATIONS
# loss = -cp.log(arg) @ coverage

# pi_lower_bound = pi >= 0
# pi_upper_bound = pi <= 1
# problem = cp.Problem(cp.Minimize(loss), [pi_lower_bound, pi_upper_bound])
# minimum_value = problem.solve(solver=cp.SCS)
# dual_objective = problem.solver_stats.extra_stats['info']['dobj']
# pi_value = float(pi.value)

# tolerance_dual = 1e-6
# tolerance_pi_value = 1e-4

# lower_bound_active = False
# upper_bound_active = False

# if abs(pi_lower_bound.dual_value) > tolerance_dual:
#     lower_bound_active = True
#     if not pi_value <= tolerance_pi_value:
#         raise RuntimeError("Lower bound active, but pi not close to 0 ({pi_value=})")

# if abs(pi_upper_bound.dual_value) > tolerance_dual:
#     upper_bound_active = True
#     if not pi_value >= 1 - tolerance_pi_value:
#         raise RuntimeError("Upper bound active, but pi not close to  1 ({pi_value=})")

# print(f"{pi_value=}")
# print(f"{lower_bound_active=}")
# print(f"{pi_lower_bound.dual_value=}")
# print(f"{upper_bound_active=}")
# print(f"{pi_upper_bound.dual_value=}")

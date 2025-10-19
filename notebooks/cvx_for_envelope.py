import numpy as np
import cvxpy as cp
import pickle

with open('coverage.pkl','rb') as f:
    coverage_mid = pickle.load(f)
    
coverage_left = np.zeros_like(coverage_mid)
coverage_left[0:10] = 3.0
coverage_uniform = np.ones_like(coverage_mid)


num_position_coverage = 100
locations = np.linspace(1/(2*num_position_coverage),
                        1 - 1/(2*num_position_coverage),
                        num_position_coverage)

coverage = coverage_uniform
pi = cp.Variable()
arg = 1 + pi - 2 * pi * locations
loss = -cp.log(arg) @ coverage

pi_lower_bound = pi >= 0
pi_upper_bound = pi <= 1
problem = cp.Problem(cp.Minimize(loss), [pi_lower_bound, pi_upper_bound])
minimum_value = problem.solve(solver=cp.SCS)
pi_value  = float(pi.value)


tolerance_dual = 1e-6
tolerance_pi_value = 1e-4

lower_bound_active = False
upper_bound_active = False

if abs(pi_lower_bound.dual_value) > tolerance_dual:
    lower_bound_active = True
    if not pi_value <= tolerance_pi_value:
        raise RuntimeError("Lower bound active, but pi not close to 0.")
        
if abs(pi_upper_bound.dual_value) > tolerance_dual:
    upper_bound_active = True
    if not pi_value >= 1 - tolerance_pi_value:
        raise RuntimeError("Upper bound active, but pi not close  1.")
        
print(f"{pi_value=}")
print(f"{lower_bound_active=}")
print(f"{pi_lower_bound.dual_value=}")
print(f"{upper_bound_active=}")
print(f"{pi_upper_bound.dual_value=}")
    



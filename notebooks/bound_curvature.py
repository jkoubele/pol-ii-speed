import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

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

plt.plot(read_locations, curvature_bounds)
plt.title("Loss curvature by read position in intron")
plt.ylabel("Upper bound on loss abs. curvature")
plt.xlabel("Read position in intron (5' to 3')")
plt.savefig('curvature_bounds.png')
plt.show()

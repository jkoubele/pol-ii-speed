import numpy as np
import matplotlib.pyplot as plt
from functools import partial
from scipy.optimize import minimize_scalar, OptimizeResult


def negative_ll(density: np.ndarray, phi: float) -> float:
    num_locations = density.size
    locations = np.linspace(start=1 / (2 * num_locations),
                            stop=1 - 1 / (2 * num_locations),
                            num=num_locations)
    loss_per_location = -np.log(1 + phi - 2 * phi * locations)
    return loss_per_location.T @ density


def estimate_phi(density: np.ndarray) -> float:
    objective_function = partial(negative_ll, density)
    return minimize_scalar(objective_function,
                           bounds=(0, 1),
                           method='bounded').x


def plot_density(density: np.ndarray, title='') -> None:
    phi = estimate_phi(density)

    num_locations = density.size
    locations = np.linspace(start=1 / (2 * num_locations),
                            stop=1 - 1 / (2 * num_locations),
                            num=num_locations)

    estimated_density = 1 + phi - 2 * phi * locations
    plt.fill_between(locations,
                     estimated_density,
                     alpha=0.4,
                     label='Estimated density')

    plt.bar(x=locations,
            height=density * estimated_density.max() / density.max(),
            width=1 / num_locations,
            label='Read locations')

    if title:
        title += '\n'

    plt.title(title + fr"$\varphi$={round(phi, 3)}")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    density_1 = [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                 0.0000, 0.0000, 0.0000, 0.0000, 0.0026, 0.0026, 0.0026, 0.0000, 0.0175,
                 0.0201, 0.0276, 0.0276, 0.0000, 0.0026, 0.0026, 0.0113, 0.0201, 0.0139,
                 0.0227, 0.0227, 0.0088, 0.0152, 0.0178, 0.0178, 0.0075, 0.0250, 0.0113,
                 0.0088, 0.0250, 0.0325, 0.0163, 0.0150, 0.0000, 0.0150, 0.0000, 0.0082,
                 0.0227, 0.0603, 0.0742, 0.0761, 0.0659, 0.0332, 0.0257, 0.0216, 0.0216,
                 0.0077, 0.0154, 0.0128, 0.0051, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                 0.0000, 0.0000, 0.0000, 0.0000, 0.0088, 0.0088, 0.0139, 0.0077, 0.0190,
                 0.0139, 0.0126, 0.0100, 0.0026, 0.0100, 0.0150, 0.0000, 0.0075, 0.0075,
                 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                 0.0000]

    density_2 = [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                 0.0000, 0.0000, 0.0000, 0.0031, 0.0253, 0.0335, 0.0335, 0.0130, 0.0319,
                 0.0297, 0.0333, 0.0599, 0.0459, 0.0506, 0.0494, 0.0380, 0.0283, 0.0445,
                 0.0391, 0.0314, 0.0247, 0.0114, 0.0063, 0.0031, 0.0031, 0.0000, 0.0000,
                 0.0031, 0.0530, 0.0680, 0.0661, 0.0099, 0.0099, 0.0135, 0.0135, 0.0135,
                 0.0135, 0.0135, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                 0.0000, 0.0139, 0.0139, 0.0000, 0.0139, 0.0139, 0.0139, 0.0139, 0.0000,
                 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                 0.0000]

    density = np.array(density_1)
    plot_density(density)

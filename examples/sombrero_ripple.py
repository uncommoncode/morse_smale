"""NOTE(emmett): It appears that this manifold with weak local minimum causes difficulties with partitioning."""
import numpy as np
import util

def gen_manifold(x_sample, y_sample):
    x_sample *= np.pi**3
    y_sample *= np.pi**3
    v = np.sqrt(x_sample**2 + y_sample**2)
    # Add a small amount to prevent divide by zero.
    epsilon = 1.0e-6
    return np.sin(v) / (epsilon + v)

if __name__ == "__main__":
    util.run_main(gen_manifold, default_n_samples=2000)

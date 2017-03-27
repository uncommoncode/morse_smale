import numpy as np
import util

def gen_manifold(n_samples):
    x_space = np.linspace(0, 1, int(np.sqrt(n_samples)))
    y_space = np.linspace(0, 1, int(np.sqrt(n_samples)))
    x_sample, y_sample = [np.array(a.flat) for a in np.meshgrid(x_space, y_space)]
    z_sample = 0.5 * (np.exp(-(x_sample - 0.25)**2/0.3**2) + np.exp(-(y_sample - 0.25)**2/0.3**2) +
                      np.exp(-(x_sample - 0.75)**2/0.1**2) + np.exp(-(y_sample - 0.75)**2/0.1**2))
    return x_sample, y_sample, z_sample

if __name__ == "__main__":
    util.plot_manifold_partitions(gen_manifold, n_samples=2000)
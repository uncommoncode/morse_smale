import numpy as np
import util

def gen_manifold(n_samples):
    x_space = np.linspace(-np.pi, np.pi, int(np.sqrt(n_samples)))
    y_space = np.linspace(-np.pi, np.pi, int(np.sqrt(n_samples)))
    x_sample, y_sample = [np.array(a.flat) for a in np.meshgrid(x_space, y_space)]
    z_sample = np.sin(x_sample) + np.cos(y_sample)
    return x_sample, y_sample, z_sample

if __name__ == "__main__":
    util.plot_manifold_partitions(gen_manifold, n_samples=300)

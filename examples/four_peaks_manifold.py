import numpy as np
import util

def gen_manifold(x_sample, y_sample):
    x_sample = x_sample * 0.5 + 0.5
    y_sample = y_sample * 0.5 + 0.5
    return 0.5 * (np.exp(-(x_sample - 0.25)**2/0.3**2) + np.exp(-(y_sample - 0.25)**2/0.3**2) +
                  np.exp(-(x_sample - 0.75)**2/0.1**2) + np.exp(-(y_sample - 0.75)**2/0.1**2))

if __name__ == "__main__":
    util.run_main(gen_manifold)

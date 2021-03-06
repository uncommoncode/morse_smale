import numpy as np
import util

def gen_manifold(x_sample, y_sample):
    x_sample *= np.pi
    y_sample *= np.pi
    return np.sin(x_sample) + np.cos(y_sample)

if __name__ == "__main__":
    util.run_main(gen_manifold)

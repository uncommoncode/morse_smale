import util

def gen_manifold(x_sample, y_sample):
    x_sample *= 1.5
    y_sample *= 1.5
    return x_sample**4 - 2.0 * x_sample**2 + y_sample**2

if __name__ == "__main__":
    util.run_main(gen_manifold)

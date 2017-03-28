"""NOTE(emmett): It appears that this manifold is not partitioned nicely."""
import util

def gen_manifold(x_sample, y_sample):
    return (x_sample + y_sample)**2

if __name__ == "__main__":
    util.run_main(gen_manifold)

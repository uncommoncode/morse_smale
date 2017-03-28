from __future__ import print_function
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import morse_smale
from sklearn.decomposition import PCA
import argparse

def _gen_samples(n_samples, random_grid):
    x_space = np.linspace(-1, 1, int(np.sqrt(n_samples)))
    y_space = np.linspace(-1, 1, int(np.sqrt(n_samples)))
    if random_grid:
        x_sample = np.random.choice(x_space, n_samples, replace=True)
        y_sample = np.random.choice(y_space, n_samples, replace=True)
    else:
        x_sample, y_sample = [np.array(a.flat) for a in np.meshgrid(x_space, y_space)]
    return x_sample, y_sample

def _plot_results(x_sample, y_sample, z_sample, result):
    # Plot results from Morse Smale partitioning
    print("Number of partitions found: {}".format(len(set(result.partitions))))

    fig = plt.figure()
    ax = fig.add_subplot(221, projection='3d')
    ax.title.set_text('Manifold Plot')
    ax.scatter(x_sample, y_sample, z_sample, c=z_sample)

    ax = fig.add_subplot(222, projection='3d')
    ax.title.set_text('Extrema Plot')
    ax.scatter(x_sample, y_sample, z_sample, c='gray', s=0.3)
    ax.scatter(x_sample[result.min_indices], y_sample[result.min_indices], z_sample[result.min_indices], c='red')
    ax.scatter(x_sample[result.max_indices], y_sample[result.max_indices], z_sample[result.max_indices], c='blue')

    ax = fig.add_subplot(223, projection='3d')
    ax.title.set_text('Partition Plot')
    ax.scatter(x_sample, y_sample, z_sample, c=result.partitions)

    plt.show()

def plot_manifold_partitions(manifold_generator, n_samples, random_grid):
    x_sample, y_sample = _gen_samples(n_samples, random_grid)
    z_sample = manifold_generator(x_sample, y_sample)

    # Create our morse smale partition.
    sample_points = np.stack([x_sample, y_sample])
    # Create our morse smale partition.
    sample_points = np.stack([x_sample, y_sample])
    result = morse_smale.nn_merged_partition(sample_points=sample_points, sample_values=z_sample, k_neighbors=15)
    _plot_results(x_sample, y_sample, z_sample, result)

def create_random_rotation(output_dimensions, input_dimensions):
    """Create a random rotation matrix that projects from R^input_dimensions to R^output_dimensions."""
    # Apply SVD to a random matrix and extract values from U, which has an orthonormal basis. This could also be done to
    # V.
    if output_dimensions < input_dimensions:
        raise ValueError("Having fewer output dimensions than input dimensions is not allowed.")
    random_projection = np.random.random((output_dimensions, output_dimensions))
    U, s, V = np.linalg.svd(np.dot(random_projection, random_projection.T))
    # Return only the number of input dimensions we need.
    return U[:, :input_dimensions]

def plot_pca_manifold_partitions(manifold_generator, n_samples, n_projected_dimensions, random_grid):
    x_sample, y_sample = _gen_samples(n_samples, random_grid)
    z_sample = manifold_generator(x_sample, y_sample)

    # Create our morse smale partition.
    sample_points = np.stack([x_sample, y_sample])

    # Create a random projection for us into R^projected_dimensions.
    projection = create_random_rotation(output_dimensions=n_projected_dimensions, input_dimensions=2)
    projected_sample_points = np.dot(projection, sample_points)

    result = morse_smale.nn_merged_partition(sample_points=projected_sample_points, sample_values=z_sample,
                                             k_neighbors=15)

    pca = PCA(n_components=2)
    pca_points = pca.fit_transform(projected_sample_points.T)

    print("Explained PCA variances: {}".format(pca.explained_variance_))

    pca_x_sample = pca_points[:, 0]
    pca_y_sample = pca_points[:, 1]

    _plot_results(pca_x_sample, pca_y_sample, z_sample, result)

def run_main(manifold_generator, default_n_samples=300):
    parser = argparse.ArgumentParser()
    parser.add_argument("--random_grid", action='store_true', help="Use random samples for manifold.")
    parser.add_argument("--n_samples", default=default_n_samples, type=int,
                        help="Number of samples to generate for the manifold.")
    parser.add_argument("--n_projected_dimensions", default=None,
                        help="Number of dimensions to project points from R^2 to R^n_projected_dimensions. Internally "
                             "projects points onto a random orthornormal basis which acts as a rotation. If not defined "
                             "no projection is applied.")
    args = parser.parse_args()
    if args.n_projected_dimensions is None:
        plot_manifold_partitions(manifold_generator, n_samples=args.n_samples, random_grid=args.random_grid)
    else:
        plot_pca_manifold_partitions(manifold_generator, n_samples=args.n_samples,
                                     n_projected_dimensions=int(args.n_projected_dimensions),
                                     random_grid=args.random_grid)

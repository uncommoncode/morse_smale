from __future__ import print_function
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import morse_smale

def plot_manifold_partitions(manifold_generator, n_samples):
    x_sample, y_sample, z_sample = manifold_generator(n_samples)

    # Create our morse smale partition.
    sample_points = np.stack([x_sample, y_sample])
    result = morse_smale.nn_merged_partition(sample_points=sample_points, sample_values=z_sample, k_neighbors=15)

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

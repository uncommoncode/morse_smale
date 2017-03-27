"""A python wrapper for the Morse Smale Complex library providing a topological exploration of manifolds.

References
----------
  [1] Samuel Gerber and Kristin Potter "The Morse-Smale Complex for Data Analysis, Journal of Statistical Software", 
  2012, vol. 50, no. 2, pp 1-22 
  
  [2] Samuel Gerber, Oliver Ruebel Peer-Timo Bremer, Valerio Pascucci, Ross Whitaker, "Morse-Smale Regression, Journal 
  of Computational and Graphical Statistics", 2012

  [3] Samuel Gerber, Peer-Timo Bremer, Valerio Pascucci, Ross Whitaker, "Visual Exploration of High Dimensional Scalar 
  Functions", IEEE Transactions on Visualization and Computer Graphics, vol. 16, no. 6, pp 1271-1280, Nov.-Dec. 2010.
"""
import _cpy_morse_smale as cms
import collections
import numpy as np

# These values are taken from the msr R package.
# The relative error bound for approximate nearest neighbor searching.
DEFAULT_ANN_EPSILON = 0.01

def _sample_points_to_data(sample_points):
    """Convert from sample points to the format expected by the C morse smale complex library."""
    return np.array(sample_points.T.flat)

def get_persistence_level(sample_values, p=0.2):
    """Return a persistence level given the percentage of the range to cover.
    
    Parameters
    ----------
    sample_values: numpy.array
         A numpy array with the scalar values for the manifold at each sample point of shape (n_samples).
    p: float
        A value between 0 and 1 describing the percentage to compute the persistance level for.
        
    Returns
    -------
    float
        The persistence_level to pass to the partition function.
    """
    return (sample_values.max() - sample_values.min()) * p

PartitionResult = collections.namedtuple("PartitionResult", ["partitions", "min_indices", "max_indices",
                                                             "persistences"])

def nn_partition(sample_points, sample_values, k_neighbors=None,  persistence_level=None, smooth_neighbors=False,
                 ann_epsilon=DEFAULT_ANN_EPSILON):
    """Partition the Morse-Smale complices and return the partitions, extrema, and persistance data.
    
    Parameters
    ----------
    sample_points: numpy.matrix
        A numpy matrix with the points specifying the location for each sample value of shape (n_dimensions, n_samples)
    sample_values: numpy.array
        A numpy array with the scalar values for the manifold at each sample point of shape (n_samples).
    k_neighbors: int
        The number of neighbors to use. Defaults to the number of dimensions.
    persistence_level: float
        Compute the Morse-Smale complex for a single persistence level. Extrema with persistence less than 
        persistance_level are ignored. Defaults to 20% of the range from the min value to the max value.
    smooth_neighbors: boolean
        If the data is very noise many extrema are introduced. If set to true the steepest ascent/descent is not 
        computed based on the raw function values but based on the function value obtained by averaging the function 
        values of the k neareast neighbors.
    ann_epsilon: float
        Specifies how close the approximate nearest neighbor approximation should be, i.e, the ratio of distance to 
        approximate neareast neighbor to true neareast neighbor is at most $1 + ann_epsilon$.
    
    Returns
    -------
    PartitionResult
        A named tuple containing the partitions ("crystals"), indices into the sample_points array for min and max
        extrema, and persistance information that hint at the stability of the result.
    """
    n_samples = sample_points.shape[1]
    n_dimensions = sample_points.shape[0]
    if n_samples != len(sample_values):
        raise ValueError("Dimensions for sample_values and sample_points dont match!")
    if k_neighbors is None:
        k_neighbors = n_dimensions
    if persistence_level is None:
        persistence_level = get_persistence_level(sample_values)
    result = cms.nearest_neighbor_morse_smale_partition(_sample_points_to_data(sample_points), sample_values,
                                                        n_dimensions, n_samples, persistence_level, k_neighbors,
                                                        smooth_neighbors, ann_epsilon)
    partitions, min_indices, max_indices, persistences = result
    return PartitionResult(partitions=partitions, min_indices=min_indices, max_indices=max_indices,
                           persistences=persistences)

def nn_merged_partition(sample_points, sample_values, k_neighbors=None, persistence_level=None, smooth_neighbors=False,
                        ann_epsilon=DEFAULT_ANN_EPSILON):
    """Partition the Morse-Smale complices merging them based on R^2 fit of linear models and return the partitions,
    extrema, and persistance data.
    
    Parameters
    ----------
    sample_points: numpy.matrix
        A numpy matrix with the points specifying the location for each sample value of shape (n_dimensions, n_samples)
    sample_values: numpy.array
        A numpy array with the scalar values for the manifold at each sample point of shape (n_samples).
    k_neighbors: int
        The number of neighbors to use. Defaults to the number of dimensions.
    persistence_level: float
        Compute the Morse-Smale complex for a single persistence level. Extrema with persistence less than 
        persistance_level are ignored. Defaults to 20% of the range from the min value to the max value.
    smooth_neighbors: boolean
        If the data is very noise many extrema are introduced. If set to true the steepest ascent/descent is not 
        computed based on the raw function values but based on the function value obtained by averaging the function 
        values of the k neareast neighbors.
    ann_epsilon: float
        Specifies how close the approximate nearest neighbor approximation should be, i.e, the ratio of distance to 
        approximate neareast neighbor to true neareast neighbor is at most $1 + ann_epsilon$.
    
    Returns
    -------
    PartitionResult
        A named tuple containing the partitions ("crystals"), indices into the sample_points array for min and max
        extrema, and persistance information that hint at the stability of the result.
    """
    n_samples = sample_points.shape[1]
    n_dimensions = sample_points.shape[0]
    if n_samples != len(sample_values):
        raise ValueError("Dimensions for sample_values and sample_points dont match!")
    if k_neighbors is None:
        k_neighbors = n_dimensions
    if persistence_level is None:
        persistence_level = get_persistence_level(sample_values)
    result = cms.nearest_neighbor_morse_smale_merged_partition(_sample_points_to_data(sample_points), sample_values,
                                                               n_dimensions, n_samples, persistence_level, k_neighbors,
                                                               smooth_neighbors, ann_epsilon)
    partitions, min_indices, max_indices, persistences = result
    return PartitionResult(partitions=partitions, min_indices=min_indices, max_indices=max_indices,
                           persistences=persistences)

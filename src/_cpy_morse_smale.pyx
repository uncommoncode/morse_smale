cimport cpy_morse_smale
cimport numpy as cnp
from libcpp.vector cimport vector
from libcpp cimport bool

import numpy as np

cdef cnp.ndarray[cnp.float64_t, ndim=1] vec_to_numpyd(vector[double] input):
  cdef cnp.ndarray[cnp.float64_t, ndim=1] data = np.zeros(input.size(), dtype=np.float64)
  for i in xrange(input.size()):
    data[i] = input[i]
  return data

cdef cnp.ndarray[cnp.int_t, ndim=1] vec_to_numpyi(vector[int] input):
  cdef cnp.ndarray[cnp.int_t, ndim=1] data = np.zeros(input.size(), dtype=np.int)
  for i in xrange(input.size()):
    data[i] = input[i]
  return data

def nearest_neighbor_morse_smale_partition(double[:] sample_points, double[:] sample_values, int n_dimensions,
                                           int n_samples, double persistence_level, int k_neighbors,
                                           bool smooth_neighbors, double ann_epsilon):
  result = cpy_morse_smale.nearest_neighbor_morse_smale_partition(&sample_points[0], &sample_values[0], n_dimensions, 
                                                                  n_samples, persistence_level, k_neighbors, 
                                                                  smooth_neighbors, ann_epsilon)
  return (vec_to_numpyi(result.partitions), vec_to_numpyi(result.min_indices),
          vec_to_numpyi(result.max_indices), vec_to_numpyd(result.persistences))

def nearest_neighbor_morse_smale_merged_partition(double[:] sample_points, double[:] sample_values, int n_dimensions,
                                            int n_samples, double persistence_level, int k_neighbors, 
                                            bool smooth_neighbors, double ann_epsilon):
  result = cpy_morse_smale.nearest_neighbor_morse_smale_merged_partition(&sample_points[0],
                &sample_values[0], n_dimensions, n_samples, persistence_level,
                k_neighbors, smooth_neighbors, ann_epsilon)
  return (vec_to_numpyi(result.partitions), vec_to_numpyi(result.min_indices),
          vec_to_numpyi(result.max_indices), vec_to_numpyd(result.persistences))

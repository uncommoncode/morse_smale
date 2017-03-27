from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "cpy_morse_smale.h":
  cdef struct PartitionResult:
    vector[int] partitions
    vector[int] min_indices
    vector[int] max_indices
    vector[double] persistences

  PartitionResult nearest_neighbor_morse_smale_partition(double *sample_points, double *sample_values, int n_dimensions,
                                                         int n_samples, double persistence_level, int k_neighbors,
                                                         bool smooth_neighbors, double ann_epsilon)
  PartitionResult nearest_neighbor_morse_smale_merged_partition(double *sample_points, double *sample_values,
                                                                int n_dimensions, int n_samples,
                                                                double persistence_level, int k_neighbors,
                                                                bool smooth_neighbors, double ann_epsilon)

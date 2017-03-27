#pragma once

#include <vector>

struct PartitionResult {
	// Partition assignment for each observation in x
	std::vector<int> partitions;
	// Indicies into scalar field of mimima for each partition.
	std::vector<int> min_indices;
	// Indicies into  scalar field of maxmima for each partition.
	std::vector<int> max_indices;
	// TODO(emmett): Still dont completly grok this
	std::vector<double> persistences;
};

PartitionResult nearest_neighbor_morse_smale_partition(double *sample_points,
											double *sample_values,
											int n_dimensions,
											int n_samples,
											double persistence_level,
											int k_neighbors, // number of neighbors to use
											bool smooth_neighbors, // used with noisy data
											double ann_epsilon // epsilon to tradeoff compute vs accuracy for nearest neighbor library
);

// merging based on R^2 of linear models
PartitionResult nearest_neighbor_morse_smale_merged_partition(double *sample_points,
													   double *sample_values,
													   int n_dimensions,
													   int n_samples,
													   double persistence_level,
													   int k_neighbors, // number of neighbors to use
													   bool smooth_neighbors, // used with noisy data
													   double ann_epsilon // epsilon to tradeoff compute vs accuracy for nearest neighbor library
);

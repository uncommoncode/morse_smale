#include "cpy_morse_smale.h"

// This file is based on the R bindings in "msr.c" provided by Samuel Gerber

#include "NNMSComplex.h"
#include "NNMSComplex2old.h"

#include "KernelDensity.h"
#include "GaussianKernel.h"

PartitionResult nearest_neighbor_morse_smale_partition(double *sample_points,
											double *sample_values,
											int n_dimensions,
											int n_samples,
											double persistence_level,
											int k_neighbors,
											bool smooth_neighbors,
											double ann_epsilon
) {
	PartitionResult result;
	
	using namespace FortranLinalg;
	
	// nearest neighbors
	//Steepest ascending KNNG(0,) and descending KNNG(1, ) neighbors for each point
	int knn = k_neighbors;
	// #columns
	int n = n_samples;
	// #rows
	int m = n_dimensions;
	//Data points
	// What is vector vs matrix?
	// #values
	double *x = sample_points;
	// #values
	double *y = sample_values;
	double pLevel = persistence_level;
	double eps = ann_epsilon;
	bool smooth = smooth_neighbors;
	
	if(knn > n){
		knn = n;
	}
	DenseMatrix<double> X(m, n, x);
	DenseVector<double> Y(n, y);
	
	NNMSComplex<double> msc(X, Y, knn, smooth, eps);
	msc.mergePersistence(pLevel);
	
	result.partitions.resize(n);
	DenseVector<int> c(n, result.partitions.data());
	msc.getPartitions(c);
	
	int nCrystals = msc.getNCrystals();
	
	result.min_indices.resize(nCrystals);
	DenseVector<int> vmins(nCrystals, result.min_indices.data());
	msc.getMin(vmins);
	
	result.max_indices.resize(nCrystals);
	DenseVector<int> vmaxs(nCrystals, result.max_indices.data());
	msc.getMax(vmaxs);
	
	int nExt = msc.getNAllExtrema();
	result.persistences.resize(nExt - 1);
	DenseVector<double> vps(nExt-1, result.persistences.data());
	msc.getPersistence(vps);
	
	msc.cleanup();
	
	return result;
}

PartitionResult nearest_neighbor_morse_smale_merged_partition(double *sample_points,
															  double *sample_values,
															  int n_dimensions,
															  int n_samples,
															  double persistence_level,
															  int k_neighbors,
															  bool smooth_neighbors,
															  double ann_epsilon
) {
	PartitionResult result;
	
	using namespace FortranLinalg;
	
	// nearest neighbors
	//Steepest ascending KNNG(0,) and descending KNNG(1, ) neighbors for each point
	int knn = k_neighbors;
	// #columns
	int n = n_samples;
	// #rows
	int m = n_dimensions;
	//Data points
	// What is vector vs matrix?
	// #values
	double *x = sample_points;
	// #values
	double *y = sample_values;
	double pLevel = persistence_level;
	double eps = ann_epsilon;
	bool smooth = smooth_neighbors;
	
	if(knn > n){
		knn = n;
	}
	DenseMatrix<double> X(m, n, x);
	DenseVector<double> Y(n, y);
	
	
	NNMSComplexR2<double> msc(X, Y, knn, -1, smooth, eps);
	
	result.partitions.resize(n);
	DenseVector<int> c(n, result.partitions.data());
	msc.getPartitions(c);
	
	int nCrystals = msc.getNCrystals();
	
	result.min_indices.resize(nCrystals);
	DenseVector<int> vmins(nCrystals, result.min_indices.data());
	msc.getMin(vmins);
	
	result.max_indices.resize(nCrystals);
	DenseVector<int> vmaxs(nCrystals, result.max_indices.data());
	msc.getMax(vmaxs);
	
	//TODO: No persistencies available in this version
	int nExt = msc.getNAllExtrema();
	result.persistences.resize(nExt - 1);
	DenseVector<double> vps(nExt-1, result.persistences.data());
	msc.getPersistence(vps);
	
	msc.cleanup();
	
	return result;
}

#ifndef MACQUEEN_H
#define MACQUEEN_H
#include "commonmacros.h"
#include <stddef.h>

/**
 * MacQueen k-means clustering algorithm.
 * points - data points
 * n - count of data points
 * d - dimension count
 * k - cluster count
 * iterations - iteration count to be done
 * clusterCenters  - cluster centers. Given as input and modified.
 * assignments - list of integers to which cluster point belongs to. Given as input and modified.
 */
void macqueen_clustering(double *points, size_t n, size_t d, size_t k, int iterations, double* clusterCenters, size_t *assignments, metricType metric);

#endif
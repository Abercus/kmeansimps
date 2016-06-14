#ifndef LLOYD_H
#define LLOYD_H
#include "commonmacros.h"
#include <stddef.h>


/**
 * Most standard k-means clustering method
 * points - data points
 * n - count of data points
 * d - dimension count
 * k - cluster count
 * iterations - iteration count to be done
 * clusterCenters  - cluster centers. Given as input and modified.
 * assignments  - list of integers to what cluster belongs to. Given as input and modified
 */
void lloyd_clustering(double *points, size_t n, size_t d, size_t k, unsigned int iterations, double *clusterCenters, size_t *assignments, metricType metric);
#endif //LLOYD_H
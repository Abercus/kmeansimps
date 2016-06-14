#ifndef ELKAN_H
#define ELKAN_H
#include "commonmacros.h"
#include <stddef.h>
/**
 * Elkan assignments method, which memorizes upper and lower bound and distances between cluster centers to decrease
 * point-center calculations count when it is not needed by using triangle inequality.
 * points - data points
 * n - count of data points
 * d - dimension count
 * k - cluster count
 * iterations - iteration count to be done
 * clusterCenters - cluster centers. Given as input and modified.
 * assignments - list of integers to which cluster point belongs to. Given as input and modified.
 */
void elkan_clustering(double *points, size_t n, size_t d, size_t  k, int iterations, double *clusterCenters, size_t *assignments, metricType metric);

#endif //ELKAN_H
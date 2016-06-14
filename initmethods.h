#ifndef BAKAC_INITMETHODS_H
#define BAKAC_INITMETHODS_H
#include "commonmacros.h"
#include <stddef.h>
/**
 * Function picks random k points as cluster centers. Called Forgy method.
 * [IN] data - data points, from which k points are picked from.
 * [IN] n - points count
 * [IN] d - dimension count
 * [IN] k - cluster count, k centers are returned
 * [IN] seed - integer value to be used as seed by pseudo RNG algorithm.
 * return - random n points from data as centers
 */
double* random_cluster_init(double* data, size_t n, size_t d, size_t k, unsigned int seed);


/**
 * Function picks first n points as cluster centers.
 * [IN] data - data points, from which k points are picked from.
 * [IN] n - points count
 * [IN] d - dimension count
 * [IN] k - cluster count, k centers are returned
 * return - first n points from data as centers.
 */
double* pick_n_first(double* data, size_t n, size_t d, size_t k);



/**
 * Function uses k-means++ initialization method to pick cluster centers.
 * [IN] data - data points, from which k points are picked from.
 * [IN] n - points count
 * [IN] k - cluster count, k centers are returned
 * return - k points from data as centers.
 */
double* kmeans_pp(double* data, size_t n, size_t d, size_t k, unsigned int seed, metricType m_metric);
 
 /**
 * Function uses furthest_first initialization method to pick cluster centers.
 * [IN] data - data points, from which k points are picked from.
 * [IN] n - points count
 * [IN] k - cluster count, k centers are returned
 * return - k points from data as centers.
 */
double* furthest_first(double* data, size_t n, size_t d, size_t k, unsigned int seed, metricType m_metric);
 
 

 /**
  * Function uses random partition initialization method to pick cluster centers.
  * [IN] data - data points, from which k centers are picked from.
  * [IN] n - points count
  * [IN] k - cluster count, k centers are returned
  * return - k points from data as centers.
  */
double* random_partition(double* data, size_t n, size_t d, size_t k, unsigned int seed);

 /**
  * Function uses uniform initialization method in order to pick cluster centers.
  * [IN] data - data points, from which k centers are picked from.
  * [IN] n - points count
  * [IN] k - cluster count, k centers are returned
  * return - k points from data as centers.
  */
double* uniform_init(double* data, size_t n, size_t d, size_t k, unsigned int seed);
  
#endif //BAKAC_INITMETHODS_H
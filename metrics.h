#ifndef METRICS_H
#define METRICS_H
#include <stddef.h>

/**
 * Euclidean distance metric
 * d - dimension count
 * X - data object 1
 * Y - data object 2
 * returns Euclidean distance between X and Y
 */
 
 
double euclid(size_t d, double* data1, double* data2);


/**
 * Taxicab / Manhattan distance metric
 * d - dimension count
 * X - data object 1
 * Y - data object 2
 * returns Manhattan distance between X and Y
 */
 
 double manhattan(size_t d, double* data1, double* data2);
 
 /**
 * Angular distance metric
 * d - dimension count
 * X - data object 1
 * Y - data object 2
 * returns angular similarity between X and Y. Range is [0, 1]
 */
 
 double angular(size_t d, double* data1, double* data2);
 
#endif
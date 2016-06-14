#ifndef BAKAC_AUXFUNCTIONS_H
#define BAKAC_AUXFUNCTIONS_H
#include "commonmacros.h"
#include <stddef.h>

#define LLOYD 1
#define ELKAN 2
#define MCQUEEN 3
#define HAMERLY 4
#define FIND_CLOSEST 5
#define HARTIGAN 6

#define RANDOM_INIT 1
#define KPP_INIT 2
#define FIRST_N 3
#define RANDOM_PARTITION 4
#define FURTHEST_FIRST 5

#define EUCLIDEAN_METRIC 1
#define MANHATTAN_METRIC 2


/**
 * Function to check if certain flag from terminal is true.
 * argc - argument count
 * argv - arguments
 * choice - chosen option
 * return: -1 if no argument found, otherwise location of the flag
 */
int get_option_location(int argc, char *argv[], char *choice);


/**
 * Function to read input data from files. Data points.
 * [IN] fileName - path of file to be openen
 * [OUT] n - datapoints count
 * [OUT] d - dimensions count
 * return value - 0 if failed to read data, otherwise points.
 */
double* get_points_from_file(const char *fileName, size_t *n, size_t *d);

/**
 * Function to print clustering.
 * [IN] clustering - which cluster does point belong to.
 * [IN] n - how many points
 */
void print_clustering(size_t *clustering, size_t n);

/**
 * Function to write data point assignments to a file. 
 * 1. posize_t data in input file is, 1. point assignment in output file
 * [IN] fileName - file name where assignments are written to.
 * [IN] clustering - assignments which are written to a file
 * [IN] n - assignments count
 * returns -1 if failed to write to a file,
 */
int write_assignments_to_file(char *fileName, size_t *clustering, size_t n);

/**
 * Function to write cluster center points to a file.
 * First line is cluster count, output format is cluster 
 * nr (corresponds to outputfile of assignments) and
 *  then all floating point numbers.
 * [IN] fileName - file name where centers are written to.
 * [IN] centers - centers which are to be written to a file
 * [IN] d - dimension count
 * [IN] k - centers count
 * return - 0 if managed to write to file
 */
int write_centers_to_file(char *fileName, double *centers, size_t d, size_t k);


/**
 * Function to get cluster centers from a file.
 * On the first line is the cluster count (k) and then
 * on each row are cluster coordinates. Dimensionality must be same as
 * input data points file.
 * We use this function to give our own centers.
 * [IN] fileName - name of a file where centers are read from
 * [IN] d - dimensionality of centers.
 * [OUT] k - cluster count
 * return - 0 if failed to load centers, otherwise points.
 */
double* get_centers_from_file(char *fileName, size_t d, size_t *k);


/**
 * Prints centers to screen.
 */ 
void print_centers(double* centers, size_t d, size_t k);


/**
 * Method to pick long random int.
 */
unsigned long long llrand();

/**
 * Method to test init methods. Only finds closest centers for each point and returns.
 * 
 */
void find_closest_centers(double *points, size_t n, size_t d, size_t k, double *clusterCenters, size_t *assignments, metricType metric);


/**
 * Method for choosing clustering algorithm.
 *
 */
int algorithm_choice(char* algorithm);


/**
 * Method for choosing initialization method.
 * Choice between k-means++ init, random point init and first n init.
 *
 */
int initmethod_choice(char* method);


/**
 * Method for choosing metric
 * Choice between Euclidean and Manhattan metrics.
 */
int metric_choice(char* metric);

/**
 * Function to calculate evalue 
 * TWGD - Total Within Group Distance
 *
 *
 */
double calculate_e_value(double *points, size_t n, size_t d, double *clusterCenters, size_t *assignments, metricType metric);

/**
 * Function to calculate squared error
 *
 *
 */
double calculate_sqrd_value(double *points, size_t n, size_t d, double *clusterCenters, size_t *assignments, metricType metric);


#endif



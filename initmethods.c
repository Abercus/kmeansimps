#include <assert.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include "initmethods.h"
#include "auxfunctions.h"
#include "commonmacros.h"


/**
Checks if an integer is already in an array.
1 if is, 0 otherwise
*/
static size_t in_array(size_t* array, size_t nr, size_t n) {
    size_t i;
    for (i = 0; i < n; ++i) {
        if (array[i] == nr)
            return 1;
    }
    return 0;
    
}

double* random_cluster_init(double* data, size_t n, size_t d, size_t k, unsigned int seed) {
    assert (k != 0 && k <= n);
    
    srand(seed);
    
    size_t *centerIds = malloc(sizeof(size_t) * n);
    
    double *centers = malloc(sizeof(double) * d * k);
    
    
    size_t picked = 0, randnr;
    
    // While haven't picked k centers.
    while (picked != k) {
        randnr = llrand() % n;
        if (!in_array(centerIds, randnr, picked)) {
            centerIds[picked++] = randnr;
        }
    }
    
    size_t i, j;
    
    for (i = 0; i < k; ++i) {
        for (j = 0; j < d; ++j) {
            centers[i*d+j] = data[centerIds[i] * d + j];
        }
    }
    
    ALFREE (centerIds);
    
    return centers;
}
 
double* pick_n_first(double* data, size_t n, size_t d, size_t k) {
    assert (k != 0 && k <= n);
    
    // Allocate memory for  cluster centers. K centers, D values per center
    
    double* centers = malloc(sizeof(double) * d * k);
    
    // Take first n points
    
    size_t* centerNrs = malloc(sizeof(size_t) * k);
    size_t i, j;
    for (i = 0; i < k; ++i) {
        centerNrs[i] = i;
    }
    
    for (i = 0; i < k; ++i) {
        for (j = 0; j < d; ++j) {
            centers[i*d+j] = data[centerNrs[i]*d+j];
        }
    }
    ALFREE(centerNrs);
    
    return centers;
    
}

double* kmeans_pp(double* data, size_t n, size_t d, size_t k, unsigned int seed, metricType m_metric) {
    assert (k != 0 && k <= n);
    
    srand(seed);
    
    double *centers = malloc(sizeof(double) * d * k);
    
    // Pick first point uniformly at random.
    size_t first_point = llrand() % n;
    
    size_t *centerIds = malloc(sizeof(size_t) * n);
    centerIds[0] = first_point;
    size_t picked = 1;
    
    size_t i, j;
    
    double *distances = malloc(sizeof(double)*n);
    for (i = 0; i < n; ++i) distances[i] = DBL_MAX;
    
    char *is_picked = calloc(n, sizeof(char)); // 0 if not picked as a center, 1 if is.
    
    // Choose next center ci, selecting ci = x' in X with probability proportional to the distance squared
    // When new center is picked calculate distance from each point to the newly picked center and memorize the lowest.
    while (picked < k) {
        // For each point find the distance to newly picked center.
        double total_distance = 0;
        for (i = 0; i < n; ++i) {
            double dist = 0;
            if (!is_picked[i]) {
                
                dist = m_metric(d, data + i*d, data + centerIds[picked-1]*d);
                total_distance += dist;
                if (dist < distances[i]) {
                    distances[i] = dist;
                }
                
            } else {
                distances[i] = 0;
            }
        }
        // Gets random float between 0 - total_distance, subtracts distances until <= 0. If <= 0 then found next center.
        double random_float = total_distance*((double)rand()/(double)RAND_MAX);
        for (i = 0; i < n; ++i) {
            if (is_picked[i]) continue;
            random_float -= distances[i];
            if (random_float <= 0) {
                centerIds[picked++] = i;
                is_picked[i] = 1;
                break;
            }
        }
    }
    // Add selected points as centers.
    for (i = 0; i < k; ++i) {
        for (j = 0; j < d; ++j) {
            centers[i*d+j] = data[centerIds[i] * d + j];
        }
    }
    
    ALFREE(centerIds);
    ALFREE(distances);
    ALFREE(is_picked);
    
    return centers;
}

double* furthest_first(double* data, size_t n, size_t d, size_t k, unsigned int seed, metricType m_metric) {
    assert (k != 0 && k <= n);
    
    srand(seed);
    
    double *centers = malloc(sizeof(double) * d * k);   
    
    // Pick first point uniformly at random. Only non-deterministic part in this one
    size_t first_point = llrand() % n;
    
    // Finds center with furthest distance from all current centers.
    
   
    size_t *centerIds = malloc(sizeof(size_t) * n);
    centerIds[0] = first_point;
    size_t picked = 1;
    
    size_t i, j;
    
    double *distances = malloc(sizeof(double)*n);
    for (i = 0; i < n; ++i) distances[i] = DBL_MAX;
    
    char *is_picked = calloc(n, sizeof(char)); // 0 if not picked as a center, 1 if is.
    
    // Choose next center ci, selecting ci = x' in X with probability ...
    // When new center is picked calculate distance from each point to that 
    
    double furthest_distance;
    size_t furthest;
    while (picked < k) {
        // For each point find the distance to newly picked center.
        for (i = 0; i < n; ++i) {
            double dist = 0;
            if (!is_picked[i]) {
                // Collect minimum distances.
                dist = m_metric(d, data + i*d, data + centerIds[picked-1]*d);
                if (dist < distances[i]) {
                    distances[i] = dist;
                }
            } else {
                distances[i] = 0;
            }
        }
        furthest = 0;
        furthest_distance = 0;
        // finding furthest from collect distances.
        for (i = 0; i < n; ++i) {
            if (distances[i] > furthest_distance) {
                furthest = i;
                furthest_distance = distances[i];
            }
        }
        is_picked[furthest] = 1;
        centerIds[picked++] = furthest;
        
    }
    
    for (i = 0; i < k; ++i) {
        for (j = 0; j < d; ++j) {
            centers[i*d+j] = data[centerIds[i] * d + j];
        }
    }
    
    ALFREE(centerIds);
    ALFREE(distances);
    ALFREE(is_picked);
    
    return centers;
}


double* random_partition(double* data, size_t n, size_t d, size_t k, unsigned int seed) {
    assert (k != 0 && k <= n);
    
    srand(seed);
    
    double *centers = calloc(d*k, sizeof(double));   
    size_t *centerIds = malloc(sizeof(size_t) * n);
    size_t *clusterSizes = calloc(k, sizeof(size_t));
    
    // Pick random center for each.
    size_t i, j, randnr;
    
    for (i = 0; i < n; ++i) {
        randnr = llrand() % k;
        centerIds[i] = randnr;
        clusterSizes[randnr]++;
    }
    // Some cluster can be left empty.. This is a quick fix.
    for (j = 0; j < k; ++j) {
        if (!clusterSizes[j]) {
            for (i = 0; i < n; ++i) {
                if (clusterSizes[centerIds[i]] > 1) {
                    clusterSizes[centerIds[i]]--;
                    clusterSizes[j]++;
                    centerIds[i] = j;
                }
            }
        }
    }
    
    // Accumulate points 
    for (i = 0; i < n; ++i) {
        for (j = 0; j < d; ++j) {
            centers[centerIds[i]*d+j] += data[i*d+j];
        }
    }
    // Divide by cluster sizes.
    for (i = 0; i < k; ++i) {
        for (j = 0; j < d; ++j) {
            centers[i*d+j] /= clusterSizes[i];
        }
    }
    
    ALFREE(centerIds)
    ALFREE(clusterSizes)
    return centers;
}


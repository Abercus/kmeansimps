#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "elkan.h"
#include "commonmacros.h"


void elkan_clustering(double *points, size_t n, size_t d, size_t k, int iterations, double *clusterCenters, size_t *assignments, metricType metric) {
    
 // allocate memory for lower, upper bounds and center to center distances
    double *lowerBound = calloc(n * k, sizeof(double));
    ALCHECK(lowerBound, "lowerBound");
    
    // centerToCenter is half of the distance, because we only need to use half.
    double *centerToCenter = malloc((k * k ) * sizeof(double));
    ALCHECK(centerToCenter, "centerToCenter");
    
    // n upper bound
    double *upperBound = malloc(n * sizeof(double));
    ALCHECK(upperBound, "upperBound");
    
    // Allocate memory to keep count of how many data points in each cluster
    size_t *clusterSizes = calloc(k, sizeof(size_t));
    ALCHECK(clusterSizes, "clusterSizes");
    
    // Allocate memory for cluster center locations
    // center minimums s(c)
    double *centerMinimums = malloc(k * sizeof(double));
    ALCHECK(centerMinimums, "centerMinimums");
    
    // For next cluster centers
    double *clusterCentersNext = malloc(k * d * sizeof(double));
    ALCHECK(clusterCentersNext, "clusterCentersNext");
    
    double *distMoved = malloc(k * sizeof(double));
    ALCHECK(distMoved, "distMoved");
    
    size_t i, ii, j, c1, dimNr;
    int reassignments = 1;
  
    double difference = 1000, CONVERGENCE = 0.0001;
    // Assign eachpoint to its closest initial center c(x) = argmin(c)d(x,c), using Lemma 1
    // To avoid redundant distance calculations
    
    int iterationNr = 0;
    
    size_t dataPointNr;
    for (dataPointNr = 0; dataPointNr < n; ++dataPointNr){
        size_t centerNr;
        for (centerNr = 0; centerNr < k; ++centerNr) {
            lowerBound[dataPointNr*k+centerNr] = 0;    
        }
        
        assignments[dataPointNr] = 0;
        upperBound[dataPointNr] = DBL_MAX;
        
        // Add to count
        clusterSizes[assignments[dataPointNr]]++;
        
    }                
    
    
    // While not converged and not enough iterations.
    while (iterationNr < iterations && difference > CONVERGENCE && reassignments > 0) {
        // Set change and reassignments count to 0.
        difference = 0;
        reassignments = 0;
    
        // Find center to center distances in the beginning of each iteration.
        for (i = 0; i < k; ++i) {
            centerMinimums[i] = DBL_MAX;
            for (j = 0; j < k; ++j) {
            centerToCenter[i*k+j] = 0;
            
                if (i != j) {
                    centerToCenter[i*k+j] = metric(d, clusterCenters + i*d, clusterCenters + j*d) / 2;
                    // We need to find minimum distance for each center.
                    if (centerToCenter[i*k+j] < centerMinimums[i]) {
                        centerMinimums[i] = centerToCenter[i*k+j];
                    }
                }
                
            }
        }
        
        
        // Iterate over all points
        for (i = 0; i < n; ++i) {
            // If upper bound for I is smaller than half of its center's distance to another center, then no change.
            if (upperBound[i] <= centerMinimums[assignments[i]]) {
                continue;
            }
            // We need to calculate upper bound once for each iteration.
            int needCalculation = 1;
            for (j = 0; j < k; ++j) {
                if (j == assignments[i]) continue;
                
                double zValue = MAX(lowerBound[i*k+j], centerToCenter[assignments[i]*k+j]);
                if (upperBound[i] <= zValue) {
                    continue;
                }
                
                // Update upper bound.
                if (needCalculation) {
                    upperBound[i] = 0;
                    for (dimNr = 0; dimNr < d; ++dimNr) upperBound[i] += pow(points[i*d+dimNr] - clusterCenters[assignments[i]*d+dimNr], 2);
                    upperBound[i] = sqrt(upperBound[i]);
                    needCalculation = 0;
                    if (upperBound[i] <= zValue) {
                        continue;
                    }

               }
               // Get distance.
               lowerBound[i*k+j] = metric(d, points + i*d, clusterCenters + j*d);
               
               // If lower bound is smaller than upper bound, then we need to change. Update upper bound as well
               if (lowerBound[i*k+j] < upperBound[i]) {
                   clusterSizes[assignments[i]]--;
                   clusterSizes[j]++;
                   assignments[i] = j;
                   reassignments++;
                   upperBound[i] = lowerBound[i*k+j];
                   
               }
            }
        }
        
        // Change accumulators to 0.
        for (i = 0; i < k*d; i++) {
            clusterCentersNext[i] = 0;
        }; 
        
        // For each center we find all the points which belong to it and add up their locations
        for (dataPointNr = 0; dataPointNr < n; ++dataPointNr) {
            size_t whichCluster = assignments[dataPointNr];
            for (dimNr = 0; dimNr < d; ++dimNr) {
                clusterCentersNext[whichCluster*d+dimNr] += points[dataPointNr*d+dimNr];
            }
            
        }
        // division by cluster size, find average in each dimension.
        for (c1 = 0; c1 < k; ++c1) {
            for (ii = 0; ii < d; ++ii) {
                clusterCentersNext[c1*d+ii] = clusterSizes[c1] ? (clusterCentersNext[c1*d+ii] / clusterSizes[c1]) : clusterCenters[c1*d+ii];
            }
            distMoved[c1] = metric(d, clusterCentersNext + c1*d, clusterCenters + c1*d);
            difference += distMoved[c1];
        }
            
        // Update upper bounds by distance moved by its corresponding center
        for (i = 0; i < n; ++i) {
            upperBound[i] += distMoved[assignments[i]];
            // Decrease lower bound by distance moved for the center and point pair.
            for (j = 0; j < k; ++j) lowerBound[i*k+j] -= distMoved[j];
        }
            
        // 7. replace each center c by m(c)
        memcpy(clusterCenters, clusterCentersNext, sizeof(double) * k * d);
        
        #if DEBUG == 0
        fprintf(stderr, "Iteration nr %d, difference: %lf, reassignments: %d\n", iterationNr, difference, reassignments);
        #endif
     iterationNr++;   
    }
            
    
    fprintf(stdout, "Elkan stopped after %d iterations\n", iterationNr);
    
    ALFREE(lowerBound);
    ALFREE(centerToCenter);
    ALFREE(upperBound);
    ALFREE(clusterSizes);
    ALFREE(centerMinimums);
    ALFREE(clusterCentersNext);
    ALFREE(distMoved);
}


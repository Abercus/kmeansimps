#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "hamerly.h"
#include "commonmacros.h"

void hamerly_clustering(double *points, size_t n, size_t d, size_t k, int iterations, double *clusterCenters, size_t *assignments, metricType metric) {
    
    //Allocate memory for:
    
    // n upper bounds
    double *upperBounds = malloc(sizeof(double)*n);
    ALCHECK(upperBounds, "upperBounds");
    
    // n lower bounds
    double *lowerBounds = malloc(sizeof(double)*n);
    ALCHECK(lowerBounds, "lowerBounds");

    // Distance to nearest other cluster for each cluster. 
    double *lowestClusterDistances = malloc(sizeof(double)*k);
    ALCHECK(lowestClusterDistances, "lowestClusterDistances");
    
    // to keep track how much distance did each cluster center move.
    double *centerDistMoved = malloc(sizeof(double)*k);
    ALCHECK(centerDistMoved, "centerDistMoved");
    
    // assignments count
    size_t  *clusterSizes = malloc(sizeof(size_t)*k);
    ALCHECK(clusterSizes, "clusterSizes");
    
    // For following iteration centers
    double *nextClusterCenters = malloc(sizeof(double)*k*d);
    ALCHECK(nextClusterCenters, "nextClusterCenters");
    
    size_t i, ii, j, dimensionNr, centerNr;
    
    // Initialize void bounds and assign each point to center 0.
    for (i = 0; i < n; ++i) {
        assignments[i] = 0;
        upperBounds[i] = DBL_MAX;
        lowerBounds[i] = 0;
    }
    
    clusterSizes[0] = n;
    for (centerNr = 1; centerNr < k; ++centerNr) clusterSizes[centerNr] = 0;
    size_t reassignments = 1;
    int iterationCount = 0;
    
    // While not converged and not enough iterations
    while (iterationCount < iterations && reassignments > 0) { 
    
        reassignments = 0;
        // Find minimum center to center distance between each two centers i != ii
        for (i = 0; i < k; i++) {
            double minDist = DBL_MAX;
            for (ii = 0; ii < k; ii++) {
                if (i != ii) {
                    // Calculate center to center distance
                    double currentDist = 0;
                    for (dimensionNr = 0; dimensionNr < d; dimensionNr++) currentDist += pow(clusterCenters[i*d+dimensionNr] - clusterCenters[ii*d+dimensionNr], 2);
                    currentDist = sqrt(currentDist) / 2;
                    if (currentDist < minDist) {minDist = currentDist; lowestClusterDistances[i] = currentDist;}
                }
            }
        }  
        
        // Next cluster centers reset
        for (i = 0; i < k*d; ++i) nextClusterCenters[i] = 0;
        // for all points
        for (i = 0; i < n; ++i) {
            // z = max(l(i), s(a(i)))
            double zValue = MAX(lowerBounds[i], lowestClusterDistances[assignments[i]]);
            // If upper bound is smaller than zVal then it is not possible that any other center is closer.
            if (upperBounds[i] <= zValue) continue;
            
             // Recalculate upper bound
            upperBounds[i] = metric(d, points + i*d, clusterCenters + assignments[i]*d);
            if (upperBounds[i] <= zValue) continue;
            
            // Find c(j) and c(j'), the two closest centers to x(i), as well as the distances to each
            size_t firstClosestCenter = 0;
            double distToFirstCenter = DBL_MAX;
            double distToSecondCenter = DBL_MAX;
            // calculate distance to each center, find two closest centers.
            
            // Iterate over all centers and find first and second closest distances.
            for (centerNr = 0; centerNr < k; centerNr++) {
                double currentDist = metric(d, points + i*d, clusterCenters + centerNr*d);
                // Find distace between the point and the center.
                if (currentDist < distToFirstCenter) {
                    distToSecondCenter = distToFirstCenter;
                    distToFirstCenter = currentDist;
                    firstClosestCenter = centerNr;
                } else if (currentDist < distToSecondCenter) {
                    distToSecondCenter = currentDist;
                }               
            
            
            }
            // If closest center is not currently assigned, then reassign. Update bounds.
            if (firstClosestCenter != assignments[i]) {
                upperBounds[i] = distToFirstCenter;
                clusterSizes[assignments[i]]--;
                assignments[i] = firstClosestCenter;
                clusterSizes[firstClosestCenter]++;
                reassignments++;
            }
            lowerBounds[i] = distToSecondCenter;
            
            
        }
        
        // Accumulate and calculate mean for each center.
        for (i = 0; i < n; i++) {
            for (j = 0; j < d; j++) {
            nextClusterCenters[assignments[i]*d+j] += points[i*d+j];
            }
        }
        
        // check for distance moved and move. Find max moved.
        double maxMoved = 0;
        for (centerNr = 0; centerNr < k; ++centerNr) {
            double currentDist = 0;
            if (clusterSizes[centerNr] > 0) {
                for (dimensionNr = 0; dimensionNr < d; dimensionNr++) {
                    nextClusterCenters[centerNr*d+dimensionNr] = clusterSizes[centerNr] ? 
                    (nextClusterCenters[centerNr*d+dimensionNr] / clusterSizes[centerNr]) : clusterCenters[centerNr*d+dimensionNr];
                }
                currentDist = metric(d, clusterCenters + centerNr*d, nextClusterCenters + centerNr*d);
            }
            centerDistMoved[centerNr] = currentDist;
            if (currentDist > maxMoved) maxMoved = currentDist;
            
        }
        // for all points update upper and lower distance bounds
        for (i = 0; i < n; ++i) {
            upperBounds[i] += centerDistMoved[assignments[i]];
            lowerBounds[i] -= maxMoved;
        }
        
        memcpy(clusterCenters, nextClusterCenters, sizeof(double) * k * d);
        
        #if DEBUG == 0
        printf("Hamerly - Iteration %d, maxmoved %lf, reassigned %d\n", iterationCount, maxMoved, reassignments);
        #endif
        iterationCount++;
    
    }
    
    printf("Hamerly stopped after %d iterations\n", iterationCount);
    
    ALFREE(upperBounds)
    ALFREE(lowerBounds)
    ALFREE(lowestClusterDistances)
    ALFREE(centerDistMoved)
    ALFREE(clusterSizes)
    ALFREE(nextClusterCenters)
}


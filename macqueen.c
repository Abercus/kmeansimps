#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include "macqueen.h"
#include "commonmacros.h"


void macqueen_clustering(double *points, size_t n, size_t d, size_t  k, int iterations, double *clusterCenters, size_t *assignments, metricType metric) {
    

        // Allocate memory for cluster sizes. How many objects in a cluster.
        size_t *clusterSizes = calloc(k, sizeof(size_t));
        ALCHECK(clusterSizes, "clusterSizes");
        
        
        
        int iteration_c = 0, changed = 1;
        size_t i, centerNr, dim;
        // First iteration is different. For each point find closest center and assign it to it.. Move centers and 
        // start usual macqueen. Each time point changes the cluster in belongs to, change both centers (previous and new).
        for (i = 0; i < n; ++i) {
            double lowest_distance = DBL_MAX;
            for (centerNr = 0; centerNr < k; ++centerNr) {
                // calculate distance
                double distance = metric(d, points + i*d, clusterCenters + centerNr*d);
                if (distance < lowest_distance) {
                    assignments[i] = centerNr;
                    lowest_distance = distance;
                }
            }
        }
        
        // Move centers to the means.
        for (i = 0; i < k*d; ++i) clusterCenters[i] = 0;
        
        for (i = 0; i < n; ++i) {
            clusterSizes[assignments[i]] += 1;
            for (dim = 0; dim < d; ++dim) {
                clusterCenters[assignments[i]*d+dim] += points[i*d+dim];
            }
        }
        
        // Divide each center's points sum by count of points in its cluster.
        for (centerNr = 0; centerNr < k; ++centerNr) {
            for (dim = 0; dim < d; ++dim) {
                if (clusterSizes[centerNr]) {
                    clusterCenters[centerNr*d+dim] = clusterCenters[centerNr*d+dim] / clusterSizes[centerNr];
                } else {
                    clusterCenters[centerNr*d+dim] = clusterCenters[centerNr*d+dim];
                    }
            }
        }
        
        // Every iteration check if point has a new closest center. If so then move point and update centers.
        while (iteration_c < iterations && changed > 0) {
            changed = 0;
            // For each point
            for (i = 0; i < n; ++i) {
                // find distance to current center and compare to other centers.

                double currentDist = metric(d, points + i*d, clusterCenters + assignments[i]*d);
                
                size_t nearestCenter = assignments[i];
                
                // find nearest center.
                for (centerNr = 0; centerNr < k; ++centerNr) {
                    if (centerNr != assignments[i]) {
                        
                        double newDist = metric(d, points + i*d, clusterCenters + centerNr*d);
                        if (newDist < currentDist) {
                            currentDist = newDist;
                            nearestCenter = centerNr;
                        }
                    }
                }
                // If this is true, then there is a new nearest center. Move centroids.
                if (nearestCenter != assignments[i]) {
                    // Assign point to the other cluster and move center.
                    for (dim = 0; dim < d; ++dim) {
                        // Subtracted from previous center.
                        // Check if any points assigned, if not, then don't move.
                        if (clusterSizes[assignments[i]] - 1 > 0) {
                            clusterCenters[assignments[i]*d+dim] *= clusterSizes[assignments[i]];
                            clusterCenters[assignments[i]*d+dim] -= points[i*d+dim];
                            clusterCenters[assignments[i]*d+dim] /= (clusterSizes[assignments[i]] - 1);
                        }
                        
                        // Add to new centroid.
                        clusterCenters[nearestCenter*d+dim] *= clusterSizes[nearestCenter];
                        clusterCenters[nearestCenter*d+dim] += points[i*d+dim];
                        clusterCenters[nearestCenter*d+dim] /= (clusterSizes[nearestCenter] + 1);
                        
                        // Change sizes and assignment
                        clusterSizes[assignments[i]] -= 1;
                        clusterSizes[nearestCenter] += 1;
                        assignments[i] = nearestCenter;
                        changed += 1;
                        
                    }
                }
                
                
            }

            
         // Set accumulators to 0
         for (centerNr = 0; centerNr < k; ++centerNr) {
             if (clusterSizes[centerNr]) {
                for (dim = 0; dim < d; ++dim) clusterCenters[centerNr*d+dim] = 0;	
             }
         }
          		
         // Add all up
         for (i = 0; i < n; ++i) {		
               for (dim = 0; dim < d; ++dim) {		
                   clusterCenters[assignments[i]*d+dim] += points[i*d+dim];		
                }		
          }		
             		
            // Divide each center's points sum by count of points in its cluster.		
          for (centerNr = 0; centerNr < k; ++centerNr) {		
                for (dim = 0; dim < d; ++dim) {		
                    if (clusterSizes[centerNr]) {		
                         clusterCenters[centerNr*d+dim] = clusterCenters[centerNr*d+dim] / clusterSizes[centerNr];		
                    } 
                 }		
             }            
            iteration_c++;
        }
        printf("MacQueen stopped after %d iterations\n", iteration_c);
            
            
        ALFREE(clusterSizes);
}

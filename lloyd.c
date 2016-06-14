#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include "lloyd.h"
#include "commonmacros.h"


void lloyd_clustering(double *points, size_t n, size_t d, size_t k, unsigned int iterations, double *clusterCenters, size_t *assignments, metricType metric) {       
        
        
        
        double CONVERGENCE = 0.0001;
        unsigned int iteration_c = 0;    
        
        // Allocate memory to keep count of how many data points in each cluster
        size_t *clusterSizes  = malloc(k * sizeof(size_t)); 
        ALCHECK(clusterSizes, "clusterSizes");          
    
        // Allocate memory for cluster center locations.
        double *clusterCentersNext = malloc(d * k * sizeof(double));
        ALCHECK(clusterCentersNext, "clusterCentersNext");        
        
        double sumPointsToCenters = DBL_MAX, prevSumPointsToCenters = 0, difference = DBL_MAX;
        
        // While less than given iteration count or converged.
        while (iteration_c < iterations && difference > CONVERGENCE) {                
            prevSumPointsToCenters = sumPointsToCenters;       
            sumPointsToCenters = 0;        
            
            // Clear center counters
            size_t dataPointNr;
            
            for (dataPointNr = 0; dataPointNr < k; ++dataPointNr) {                
                clusterSizes[dataPointNr] = 0;               
                size_t centerStart = dataPointNr*d, dimensionNr;
                for (dimensionNr = 0; dimensionNr < d; ++dimensionNr) {
                    clusterCentersNext[centerStart+dimensionNr] = 0;
                }
            }                                  
            
            // Step 1 - Find assignments for each cluster center.   
            // Iterate over all points
            for (dataPointNr = 0; dataPointNr < n; ++dataPointNr){              
                double minimumDistance = DBL_MAX;              
                size_t centerNr, dimensionNr;
                // Find closest center. 
                for (centerNr = 0; centerNr < k; ++centerNr) {                 
                    double currentDistance = 0; 
                    currentDistance = metric(d, points + dataPointNr*d, clusterCenters + centerNr*d);
                    // If current center is closer than previously found closest.
                    if (currentDistance < minimumDistance) {                      
                        minimumDistance = currentDistance;
                        assignments[dataPointNr] = centerNr;                    
                    }                                               
                }                
                // Add to count and to cluster.
               size_t assignedCluster = assignments[dataPointNr], clusterStart = assignedCluster*d, dataStart = dataPointNr*d;
               for (dimensionNr = 0; dimensionNr < d; ++dimensionNr) {
                    clusterCentersNext[clusterStart+dimensionNr] += points[dataStart+dimensionNr];
                }
                // Add to count
                clusterSizes[assignedCluster]++;
                
                // Add to distance (E value) to check for convergence.
                sumPointsToCenters += minimumDistance;
            }
            // Step 2 - Find new cluster centers based on new assignment
            size_t centerNr;
            for (centerNr = 0; centerNr < k; ++centerNr) {
                size_t clusterStart = centerNr*d, dimensionNr;
                // Fail safe. If no point in cluster, leave cluster at its own location.
                for (dimensionNr = 0; dimensionNr < d; ++dimensionNr) {
                    if (clusterSizes[centerNr]) {
                        clusterCenters[clusterStart+dimensionNr] = clusterCentersNext[clusterStart+dimensionNr] / clusterSizes[centerNr];
                    } else {
                        clusterCenters[clusterStart+dimensionNr] = clusterCenters[clusterStart+dimensionNr];
                    }
                }
            }
            
            // Check the E value difference comparing to last iteration
            difference = fabs(prevSumPointsToCenters - sumPointsToCenters);

            #if DEBUG == 0
            fprintf(stdout, "Iteration %d, E-value %lf\n", iteration_c, sqrt(sumPointsToCenters));
            #endif
            iteration_c++;
        }
    
    fprintf(stdout, "Lloyd stopped after %d iterations\n", iteration_c);
    
    // Free memory.
    ALFREE(clusterSizes);
    ALFREE(clusterCentersNext)
    
}
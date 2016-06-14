#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "hartiganwong.h"
#include "commonmacros.h"


void hartigan_clustering(double *points, size_t n, size_t d, size_t k, int iterations, double *clusterCenters, size_t *assignments, metricType metric) {

        // Allocate memory for:
        
        // How many points are in each luster
        size_t *clusterSizes = calloc(k, sizeof(size_t));
        ALCHECK(clusterSizes, "clusterSizes");
        
        
        // What is the first closest center for each point
        size_t *firstClosestCenters =  malloc(sizeof(size_t)*n);
        ALCHECK(firstClosestCenters, "firstClosestCenters");
        
        // What is the second closest center for each point
        size_t *secondClosestCenters = malloc(sizeof(size_t)*n);
        ALCHECK(secondClosestCenters, "secondClosestcCnters");
        
        
        // Which centers are in the liveset. 1 if live, 0 otherwise.
        int *isLiveSet = malloc(sizeof(int)*k);
        ALCHECK(isLiveSet, "isLiveSet");
        
        // What centers are made live in quicktransfer set.
        int *qtranLiveSet = malloc(sizeof(int)*k);
        ALCHECK(qtranLiveSet, "qtranLiveSet");
        
        // Swap memory for live sets.
        int *nextLiveSet = malloc(sizeof(int)*k);
        ALCHECK(nextLiveSet, "nextLiveSet");
        
        
        // What are R1 values for each point.
        double *R1_values = malloc(sizeof(double)*n);
        ALCHECK(R1_values, "R1_values");
        
        // To remember previous clusters.
        double *prevClusterCenters = malloc(sizeof(double) * k * d);
        ALCHECK(prevClusterCenters, "prevClusterCenters");
        
        memcpy(prevClusterCenters, clusterCenters, sizeof(double) * k * d);
        
        size_t i, centerNr, dim;

        int iteration_c = 0;

        
        // For each point, find closest and second closest centers.
        for (i = 0; i < n; ++i) {
            double fClosestDist = DBL_MAX;
            double sClosestDist = DBL_MAX;
            size_t fClosest = -1;
            size_t sClosest = -1;
            
            // Iterate over all centers
            for (centerNr = 0; centerNr < k; ++centerNr) {
                double distPointToCenter = metric(d, points + i*d, clusterCenters + centerNr*d);
               
                // Find distance from point to center.

                // Compare found distance to previously found distances.
                if (distPointToCenter < fClosestDist) {
                    // If newly found center is closer then previously found closest center, then we have a new closest center.
                    sClosestDist = fClosestDist;
                    sClosest = fClosest;
                    fClosestDist = distPointToCenter;
                    fClosest = centerNr;
                } else if (distPointToCenter < sClosestDist) {
                    // If newly found center is closer than 2nd previous center, then we got new 2nd previous center.
                    sClosestDist = distPointToCenter;
                    sClosest = centerNr;
                }
            }
            // After iterating over all set dists and assignments.
            firstClosestCenters[i] = fClosest;
            secondClosestCenters[i]= sClosest;
            assignments[i] = fClosest;
            clusterSizes[fClosest]++;
            
        }
        
        // Step 2. Update cluster centers to be averages of points contained within them.
        // All accumulators to 0.0
        for (i = 0; i < k*d; ++i) {
            clusterCenters[i] = 0.0;
        }
        
        // Find mean
        for (i = 0; i < n; ++i) {
            for (dim = 0; dim < d; ++dim) {
                clusterCenters[assignments[i]*d+dim] += points[i*d+dim];
            }
        }
        
        // Divide each center accumulator by count of points in center.
        for (centerNr = 0; centerNr < k; ++centerNr) {
            for (dim = 0; dim < d; ++dim) {
                // If at least 1 point in cluster.
                if (clusterSizes[centerNr]) {
                    clusterCenters[centerNr*d+dim] = clusterCenters[centerNr*d+dim] / clusterSizes[centerNr];
                } else {
                    //printf("Warning: Empty cluster when finding first point-to-center!\n");
                    clusterCenters[centerNr*d+dim] = prevClusterCenters[centerNr*d+dim];
                }
            }
        }
        
        // Precalculate R1 for each. Other R1 calculations are only when need to recalculate.
        
        // Step 3
        // All centers belong to the live set.
        for (i = 0; i < k; ++i) {
            isLiveSet[i] = 1;
        }
        
        // For n iterations do.
        while (iteration_c < iterations) {
            int changed = 0;
            
            // Step 4. OPTRA stage (Optimal transfer)
            
            // Clear next live set
            for (i = 0; i < k; ++i) nextLiveSet[i] = 0;
            
            
            // For each point
            for (i = 0; i < n; ++i) {
                size_t currentAssignment = assignments[i];
                // If point belongs to a cluster which is live then do 4a
                if (isLiveSet[currentAssignment]) {
                    // Compute minimum of R2 over all clusters except current.
                    double minR2 = DBL_MAX;
                    int L2 = -1;
                    for (centerNr = 0; centerNr < k; ++centerNr) {
                        if (currentAssignment == centerNr) { 
                            continue;
                        }
                        // Find distance from point to center.
                        double distancePointToCenter = metric(d, points + i*d, clusterCenters + centerNr*d);
                        distancePointToCenter *= distancePointToCenter;
                        // Calculate how much it would affect WCSS if we'd add this point to centerNr
                        double R2 = ((double) clusterSizes[centerNr] * distancePointToCenter) / (clusterSizes[centerNr] + 1);        
                        // If smaller than currently found.
                        if (R2 < minR2) {
                            minR2 = R2;
                            L2 = centerNr;
                        }
                        
                        
                    }
                    // Having found L2 with smallest R2. If this value is greater than or equal to NC(L1) * distance / (NC(L1) -1) then no reallocate
                    // Find R1
                    // Find distance of L1 to its cluster
                   
                    if (isLiveSet[currentAssignment]) {
                        // 4a. Current center is in the live set. Iterate over all centers.
                        
                        double distanceL1Center = metric(d, points + i*d, clusterCenters + currentAssignment*d);
                        distanceL1Center *= distanceL1Center;
                        // How much would WCSS change if we'd remove this center from previous center.
                        if (clusterSizes[currentAssignment] - 1 > 0) {
                            R1_values[i] = ((double) clusterSizes[currentAssignment] * distanceL1Center) / (clusterSizes[currentAssignment] - 1);
                        } else {
                            R1_values[i] = 0;
                        }
                    }
                    double R1 = R1_values[i];
                    // If WCSS would not decrease if switching cluster, then don't switch
                    if (minR2 >= R1) {
                        secondClosestCenters[i] = L2;
                    } else {
                        // Reallocation necessary, I is allocated to L2 and L1 is new IC2
                        secondClosestCenters[i] = firstClosestCenters[i];
                        firstClosestCenters[i] = L2;
                        
                        // Reallocation from one center to another.
                        // Take from previous cluster and update center
                        for (dim = 0; dim < d; ++dim) {
                            if (clusterSizes[currentAssignment] - 1 > 0) {
                                    clusterCenters[currentAssignment*d+dim] *= (double)clusterSizes[currentAssignment];
                                    clusterCenters[currentAssignment*d+dim] -= points[i*d+dim];
                                    clusterCenters[currentAssignment*d+dim] /= ((double)clusterSizes[currentAssignment] - 1);
                            }

                            // Add to new cluster and update center.
                            clusterCenters[L2*d+dim] *= (double)clusterSizes[L2];
                            clusterCenters[L2*d+dim] += points[i*d+dim];
                            clusterCenters[L2*d+dim] /= ((double)clusterSizes[L2] + 1);
                        }
                        
                        
                        // Change sizes, assignments and add centers in action to live set.
                        assignments[i] = L2;
                        
                        clusterSizes[currentAssignment] -= 1;
                        clusterSizes[L2] += 1;
                        
                        nextLiveSet[currentAssignment] = 1;
                        nextLiveSet[L2] = 1;
                        
                        changed += 1;
                        
                    }
                    
                } else {
                    // 4b Current center is not in the live set. Iterate only on over centers which are in live set.
                    double minR2 = DBL_MAX;
                    size_t L2 = -1;
                    for (centerNr = 0; centerNr < k; ++centerNr) {
                        // If center is not live then skip.
                        if (!isLiveSet[centerNr] || currentAssignment == centerNr) { 
                            continue;
                        }
                        double distancePointToCenter = metric(d, points + i*d, clusterCenters + centerNr*d);
                        distancePointToCenter *= distancePointToCenter;
                        // Calculate R2 for the center
                        double R2 = ((double) clusterSizes[centerNr] * distancePointToCenter) / (clusterSizes[centerNr] + 1);
                        
                        // If smaller than currently found.
                        if (R2 < minR2) {
                            minR2 = R2;
                            L2 = centerNr;
                        }
                        
                        
                    }
                    // Having found L2 with smallest R2. If this value is greater than or equal to NC(L1) * distance / (NC(L1) -1) then no reallocate
                    // Find R1
                    // Find distance of L1 to its cluster
                    
                    double R1 = R1_values[i];
                    // No reallocation necessary, L2 is new IC2
                    if (minR2 >= R1) {
                        secondClosestCenters[i] = L2;
                    } else {
                        // Reallocation necessary, I is allocated to L2 and L1 is new IC2
                        secondClosestCenters[i] = firstClosestCenters[i];
                        firstClosestCenters[i] = L2;
                        
                        // Reallocation
                        // Take from previous
                        for (dim = 0; dim < d; ++dim) {
                            if (clusterSizes[currentAssignment] - 1 > 0) {
                                    clusterCenters[currentAssignment*d+dim] *= (double)clusterSizes[currentAssignment];
                                    clusterCenters[currentAssignment*d+dim] -= points[i*d+dim];
                                    clusterCenters[currentAssignment*d+dim] /= ((double)clusterSizes[currentAssignment] - 1);
                            }

                            // Add to new.
                            clusterCenters[L2*d+dim] *= (double)clusterSizes[L2];
                            clusterCenters[L2*d+dim] += points[i*d+dim];
                            clusterCenters[L2*d+dim] /= ((double)clusterSizes[L2] + 1);
                        }
                        
                        
                        // Change assignments, add to cluster size.
                        assignments[i] = L2;
                        
                        clusterSizes[currentAssignment] -= 1;
                        clusterSizes[L2] += 1;
                        
                        nextLiveSet[currentAssignment] = 1;
                        nextLiveSet[L2] = 1;
                        
                        changed += 1;
                        
                    }
                }
            }
            
            // Step 5. Stop if live set is empty // Therefore no changes have occured 
            if (!changed) {
                break;
            }
            for (i = 0; i < k; ++i) {
                qtranLiveSet[i] = 0;
            }
        
            // Step 6. Quick-transfer (QTRAN) stage.
            do {
                if (changed == 2) {
                    break;
                }
                memcpy(isLiveSet, nextLiveSet, sizeof(int)*k);  
                for (i = 0; i < k; ++i) {
                    nextLiveSet[i] = 0;
                }
                changed = 0;
                for (i = 0; i < n; ++i) {
                    
                    size_t L1 = firstClosestCenters[i];
                    size_t L2 = secondClosestCenters[i];
                    
                    // No need to check if they haven't changed in previous n steps.
                    if (!isLiveSet[L1] && !isLiveSet[L2]) {
                        continue;
                    }
                    // Compute values R1 and R2. Need to calculate R1 only if it has changed.
                    if (isLiveSet[L1]) {
                        
                        double distanceL1Center = metric(d, points + i*d, clusterCenters + L1*d);
                        distanceL1Center *= distanceL1Center;

                        if (clusterSizes[L1] - 1 > 0) {
                            R1_values[i] = ((double) clusterSizes[L1] * distanceL1Center) / (clusterSizes[L1] - 1);
                        } else {
                            R1_values[i] = 0;
                        }
                    }             
                    double R1 = R1_values[i], distanceL2Center = 0;
                    
                    // Calculate distance between point and L2.
                    for (dim = 0; dim < d; ++dim) {
                        distanceL2Center += pow(points[i*d+dim] - clusterCenters[L2*d+dim],2);
                    }
                    
                    double R2 = ((double) clusterSizes[L2] * distanceL2Center) / (clusterSizes[L2] + 1);
                    
                    // If WCSS would decrease, then change clusters.
                    if (R1 < R2) {
                        continue;
                    } else {
                        // Take from previous
                        for (dim = 0; dim < d; ++dim) {
                            if (clusterSizes[L1] - 1 > 0) {
                                    clusterCenters[L1*d+dim] *= clusterSizes[L1];
                                    clusterCenters[L1*d+dim] -= points[i*d+dim];
                                    clusterCenters[L1*d+dim] /= (clusterSizes[L1] - 1);
                            }

                            // Add to new.
                            clusterCenters[L2*d+dim] *= clusterSizes[L2];
                            clusterCenters[L2*d+dim] += points[i*d+dim];
                            clusterCenters[L2*d+dim] /= (clusterSizes[L2] + 1);
                        }
                        
                        assignments[i] = L2;
                        
                        clusterSizes[L1] -= 1;
                        clusterSizes[L2] += 1;
                        
                        nextLiveSet[L1] = 1;
                        nextLiveSet[L2] = 1;
                        
                        if (!qtranLiveSet[L1]) {
                            qtranLiveSet[L1] = 1;
                            changed++;
                        }
                        if (!qtranLiveSet[L2]) {
                            qtranLiveSet[L2] = 1;
                            changed++;
                        }

                        firstClosestCenters[i] = L2;
                        secondClosestCenters[i] = L1;
                    }
                    
                }
            } while (changed > 0);
            
            memcpy(isLiveSet, qtranLiveSet, sizeof(int)*k);  
            iteration_c++;
        }
        
        printf("Hartigan-Wong stopped after %d iterations\n", iteration_c);
    
        ALFREE(clusterSizes);
        ALFREE(firstClosestCenters);
        ALFREE(secondClosestCenters);
        ALFREE(isLiveSet);
        ALFREE(nextLiveSet);
        ALFREE(R1_values);
        ALFREE(prevClusterCenters);
}
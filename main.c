#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "auxfunctions.h"
#include "clusteringalgorithms.h"
#include "initmethods.h"
#include "commonmacros.h"
#include "metrics.h"
#include <math.h>
#include <stddef.h>

int main(int argc, char *argv[]) {
    
    clock_t begin, end; 
    int time_spent;
    // Check if help flag is used. If so, then display help menu.
    if (argc < 2 || get_option_location(argc, argv, "-help") > 0 || get_option_location(argc, argv, "-h") > 0) {
        printf("Help menu. Argument and its value must be separated by a space:\n");
        printf("-h  - help menu\n");
        printf("-f  - input data points file name\n");
        printf("-o  - output file name\n");
        printf("-ci - centers input file name\n");
        printf("-a  - algorithm choice (lloyd, elkan, hamerly, macqueen, hartigan, closest)\n");
        printf("-i  - initialization method choice (kpp, forgy, partition, furthest)\n");
        printf("-m  - metric choice (euclidean, manhattan)\n");
        printf("-k  - clusters count, -ci flag has a higher priority\n");
        printf("-s  - random seed nr, otherwise uses current time as the seed. Used to confirm clustering results\n");
        printf("-n  - iteration count (default 100)\n");
        exit(0);
    }
    
    
    // Init default values.
    int algorithm = -1;
    int initmethod = -1;
    int metric = -1;
    int iterationCount = 100;
    
    // Which algorithm is used
    int algorithmChoiceFlagInd = get_option_location(argc, argv, "-a");
    if (algorithmChoiceFlagInd > 0) {
        // Check if enough arguments given
        if (algorithmChoiceFlagInd >= argc - 1) {
            printf("Give an algorithm as -a parameter!\n");
            exit(1);
        } else {
            // Check if correct algorithm given
            algorithm = algorithm_choice(argv[algorithmChoiceFlagInd + 1]);
            if (algorithm == -1) {
                printf("Incorrect parameter for -a (algorithm)\n");
                exit(1);
            }
        }
    } else {
        printf("No algorithm chosen! Will be using standard Lloyd clustering algorithm as default\n");
        algorithm = LLOYD;
    }
    
    
    // Which metric is used
    int metricChoiceFlagInd = get_option_location(argc, argv, "-m");
    if (metricChoiceFlagInd > 0) {
        // Check if enough arguments given
        if (metricChoiceFlagInd >= argc - 1) {
            printf("Give a metric as -m parameter!\n");
            exit(1);
        } else {
            // Check if correct metric given
            metric = metric_choice(argv[metricChoiceFlagInd+1]);
            if (metric == -1) {
                printf("Incorrect parameter for -m (metric)\n");
                exit(1);
            }
        }
    } else {
        printf("No metric chosen! Will be using Euclidean metric as default\n");
        metric = EUCLIDEAN_METRIC;
    }

    
    
    metricType m_metric;
    if (metric == EUCLIDEAN_METRIC) {
        m_metric = &euclid;
        printf("Using Euclidean metric\n");
    } else if (metric == MANHATTAN_METRIC) {
        m_metric = &manhattan;
        printf("Using Manhattan metric\n");
    } else {
        m_metric = &euclid;
    }
    

    begin = clock();
    
    // Allocate memory for file names.
    char *outputFilename = malloc(sizeof(char) * 100);
    ALCHECK(outputFilename, "outputFilename");
    char *outputCentersFilename = malloc(sizeof(char) * 110);
    ALCHECK(outputCentersFilename, "outputCentersFilename");
    outputFilename[0] = 0;
    outputCentersFilename[0] = 0;
    
    
    // Find outputfile flag. If no flag, then result assignments is output to screen.
    int outputfileFlagInd = get_option_location(argc, argv, "-o");
    
    if (outputfileFlagInd > 0) {
        if (outputfileFlagInd >= argc - 1) {
            printf("Give output filename as -o parameter!\n");
            exit(1);
        } else {
            strcpy(outputFilename, argv[outputfileFlagInd + 1]);
            
            strcpy(outputCentersFilename, outputFilename);
            strcat(outputCentersFilename, "_centers");
        }
    }
    
    char *inputPointsFilename = malloc(sizeof(char) * 100);
    ALCHECK(inputPointsFilename, "inputPointsFilename");
    
    inputPointsFilename[0] = 0;
    
    
    // Find inputfile flag and get it from terminal.
    int inputfileFlagInd = get_option_location(argc, argv, "-f");
    
    if (inputfileFlagInd > 0) {
        if (inputfileFlagInd >= argc - 1) {
            printf("Give input points filename as -f parameter!\n");
            exit(1);
        } else {
            strcpy(inputPointsFilename, argv[inputfileFlagInd + 1]);
        }

    } else {
        printf("Give input points filename as -f parameter!\n");
        exit(1);
    }
    
    // We need pointers here, because we are getting n and d from input files.
    size_t *n = malloc(sizeof(size_t));
    size_t *d = malloc(sizeof(size_t));
    ALCHECK(n, "integer n");
    ALCHECK(d, "integer d");
    
    
    // Getting datapoints from an input file.
    double *datapoints = 0;
    datapoints = get_points_from_file(inputPointsFilename, n, d);
    
    
    // Display n and d from objects file.
    if (datapoints) {
       printf("Data loaded successfully\n");
       printf("n: %zu, d: %zu\n", *n, *d);
             
    } else {
        printf("Failed to load data from objets file\n");
        exit(1);
    }

    end = clock();
    time_spent = (end - begin) * 1000 / CLOCKS_PER_SEC;
    
    printf("Data points loading time: %d seconds %d milliseconds\n", time_spent/1000, time_spent%1000);
    
    // if flag -ci then load cluster centers from given file, else look for flag -k.
    // If K given then use initialization method (k-means++ or random)
    // need pointer here, if k is from file.
    size_t *k = malloc(sizeof(size_t));
    ALCHECK(k, "integer k");
    *k = -1;
    
    int initMethodFlagInd = get_option_location(argc, argv, "-i");
    
    if (initMethodFlagInd > 0) {
        if (initMethodFlagInd >= argc - 1) {
            printf("Give an initialization method as -i parameter!\n");
            exit(1);
        } else {
            // Check if correct init method given
            initmethod = initmethod_choice(argv[initMethodFlagInd + 1]);
            if (initmethod == -1) {
                printf("Incorrect parameter for -i (initialization method)\n");
                exit(1);
            }
        }
    } else {
        initmethod = RANDOM_INIT;
    }
    
    // Getting centers from file or K from console.
    double* clusterCenters = 0;
    
    int centerFlagInd = get_option_location(argc, argv, "-ci");
    // Search for center flag and get centers from file, if given.
    if (centerFlagInd > 0) {
        if (centerFlagInd >= argc - 1) {
            printf("Give input cluster centers filename as -ci parameter!\n");
            exit(1);
        } else {
            printf("Loading cluster centers from a file\n");
            clusterCenters = get_centers_from_file(argv[centerFlagInd + 1], *d, k);
            if (!clusterCenters) {
                printf("Failed to load cluster centers from file! Exiting.\n");
                exit(1);
            }
        }
    } else {
        // If k is described as a program param, then use it.
        int kFlagInd = get_option_location(argc, argv, "-k");
        if (kFlagInd > 0) {
            if (kFlagInd >= argc - 1) {
                printf("Number of centers not given. Use -k parameter!\n");
                exit(1);
            } else {
                *k = atoi(argv[kFlagInd+1]);
                if (*k < 1) {
                    printf("-k parameter must be integer!\n");
                    exit(1);
                }
                if (*n < *k) {
                    printf("n (datapoints count) must greater or equal to k (centers count)!\n");
                    exit(1);
                }
                
                // Use seed as an input if given from console.
                int isSeed = get_option_location(argc, argv, "-s");
                unsigned int seed = 0;
                if (isSeed > 0) {
                    if (isSeed >= argc - 1) {
                        printf("Give seed nr as -s parameter!\n");
                        exit(1);
                    } else {
                        seed = atoi(argv[isSeed+1]);
                        printf("Using seed: %u\n", seed);
                    }
                } else {
                    seed = (unsigned int) time(NULL);
                    printf("Using random seed: %u\n", seed);
                }
                
                // Which center init method to use
                begin = clock();
                if (initmethod == RANDOM_INIT) {
                    printf("Using Forgy center initialization\n");
                    clusterCenters = random_cluster_init(datapoints, *n, *d, *k, seed);
                } else if (initmethod == KPP_INIT) {
                    printf("Using k-means++ initialization\n");
                    clusterCenters = kmeans_pp(datapoints, *n, *d, *k, seed, m_metric);
                } else if (initmethod == FIRST_N) {
                    printf("Using pick first k points as centers initialization\n");
                    clusterCenters = pick_n_first(datapoints, *n, *d, *k);
                } else if (initmethod == RANDOM_PARTITION) {
                    printf("Using random partition for centers initialization\n");
                    clusterCenters = random_partition(datapoints, *n, *d, *k, seed);
                } else if (initmethod == FURTHEST_FIRST) {
                    printf("Using furthest first for centers initialization\n");
                    clusterCenters = furthest_first(datapoints, *n, *d, *k, seed, m_metric);
                }
                end = clock();
                time_spent = (end - begin) * 1000 / CLOCKS_PER_SEC;
                printf("Center initialization time: %d seconds %d milliseconds\n", time_spent/1000, time_spent%1000);
                
            }
        } else {
            // No K or CI, 
            printf("Error: Either give cluster count or k. See -h for more\n");
            exit(1);
        }
    }
    
    printf("k: %zu\n", *k);
    // Iteration count flag
    int iterationFlagInd = get_option_location(argc, argv, "-n");
    if (iterationFlagInd > 0) {
        if (iterationFlagInd >= argc - 1) {
            printf("Give argument for iteration count -n flag,\n");
            exit(1);
        } else {
            iterationCount = atoi(argv[iterationFlagInd+1]);
            if (iterationCount < 1) {
                printf("Invalid iteration count given with -n flag. Must be an integer greater than 0\n");
                exit(1);
            } 
        }
    }
    
   
    size_t* assignments = malloc(sizeof(size_t) * (*n));
    ALCHECK(assignments, "assignments");
   
    // Begin timing.
    begin = clock();
    // Which algorithm used.
    switch(algorithm) {
        case LLOYD:
            printf("Doing Lloyd clustering\n");
            lloyd_clustering(datapoints, *n, *d, *k, iterationCount, clusterCenters, assignments, m_metric);
            break;
        case ELKAN:
            printf("Doing Elkan clustering\n");
            elkan_clustering(datapoints, *n, *d, *k, iterationCount, clusterCenters, assignments, m_metric);
            break;
        case MCQUEEN:
            printf("Doing MacQueen clustering\n");
            macqueen_clustering(datapoints, *n, *d, *k, iterationCount, clusterCenters, assignments, m_metric);
            break;
        case HAMERLY:
            printf("Doing Hamerly clustering\n");
            hamerly_clustering(datapoints, *n, *d, *k, iterationCount, clusterCenters, assignments, m_metric);
            break;
        case FIND_CLOSEST:
            printf("Finding closest clusters.\n");
            find_closest_centers(datapoints, *n, *d, *k, clusterCenters, assignments, m_metric);
            break;
        case HARTIGAN:
            printf("Doing Hartigan-Wong clustering\n");
            hartigan_clustering(datapoints, *n, *d, *k, iterationCount, clusterCenters, assignments, m_metric);
            break;
        default:
            printf("Chosen algorithm is not implemented\n");
            exit(1);
    }
    end = clock();
    time_spent = (end - begin) * 1000 / CLOCKS_PER_SEC;
    
    
    
    
    // Print output data.
    printf("Clustering time: %d seconds %d milliseconds\n", time_spent/1000, time_spent%1000);
    printf("Total Within Group Distance (TWGD) for clustering is: %lf\n", calculate_e_value(datapoints, *n, *d, clusterCenters, assignments, m_metric));
    printf("Within Cluster Sum of Squares (WCSS) for clustering is: %lf\n", calculate_sqrd_value(datapoints, *n, *d, clusterCenters, assignments, m_metric));
    printf("Done with clustering, start writing output to a file\n");
    
    
    
    
    if (outputFilename[0] == 0) {
        print_clustering(assignments, *n);
    } else {
        int success = write_assignments_to_file(outputFilename, assignments, *n);
        if (success == -1) {
            printf("Failed to write objects assignments to file\n");
            exit(1);
        } else {
            printf("Successfully wrote objects assignments to file\n");
        }
        success = write_centers_to_file(outputCentersFilename, clusterCenters, *d, *k);
        if (success == -1) {
            printf("Failed to write cluster centers info to file\n");
            exit(1);
        } else {
            printf("Successfully wrote cluster centers to a file\n");
        }
    }
    
    
    ALFREE(n);
    ALFREE(d);
    ALFREE(k);
    ALFREE(datapoints);
    ALFREE(outputFilename);
    ALFREE(inputPointsFilename);
    ALFREE(outputCentersFilename);
   
    return 0 ;
}
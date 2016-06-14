#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include "auxfunctions.h"
#include "commonmacros.h"
#include <stddef.h>
/** 
Strcmp implementation http://www.cplusplus.com/forum/general/17683/
*/
static int mystrcmp(const char *s1, const char *s2) {
        int i = 0;
        do {
            if(s1[i] != s2[i])
                return 1;
        } while(s1[i++] != '\0');
        return 0;
}
int get_option_location(int argc, char *argv[], char *choice) {

    int opt_nr;
    for (opt_nr = 1; opt_nr < argc; ++opt_nr) {
        if ((strcmp(choice, argv[opt_nr]) == 0) || (mystrcmp(choice, argv[opt_nr]) == 0)) {
            return opt_nr;
        }
    }
    // No option found
    return -1;

}

double* get_points_from_file(const char *fileName, size_t *n, size_t *d) {
   
    FILE *ifp;
    ifp = fopen(fileName, "r");
    
    if (ifp == NULL) {
        return 0;
    }
    // In our file format, data points count and dimensions are on the first row. 
    fscanf(ifp, "%zu %zu\n", n, d);
    if (*n <= 0 || *d <= 0) {
        printf("n and d must be greater than 0! Incorrect datapoints file\n");
        exit(1);
    }
    
    // According to the first row, we allocate memory.
    
    double* datapoints = malloc(sizeof(double) * (*n) * (*d));
    ALCHECK(datapoints, "datapoints");
    
    
    size_t i, j, d2 = *d;
    for (i = 0; i < *n; ++i) {
        for (j = 0; j < d2; ++j) {
            if (feof(ifp)) {
                printf("Incorrect datapoints file, not enough points!\n");
                exit(1);
            }
            fscanf(ifp, "%lf", &datapoints[i*d2+j]);
        }
    }
/*
    double test;
    fscanf(ifp, "%lf", &test);
    
    if (!feof(ifp)) {
        printf("More data in file. Incorrect datapoints file!\n");
        exit(1);
    } */
        
    fclose(ifp);
    
    return datapoints;
}


void print_clustering(size_t *assignments, size_t n) {
    size_t i;
    for (i = 0; i < n; i++) {
        printf("%zd\n", assignments[i]);
    }
}

int write_assignments_to_file(char *fileName, size_t *assignments, size_t n)  {
    FILE *ofp;
    
    ofp = fopen(fileName, "w");
    if (ofp == NULL) {
        printf("Error opening output file!");
        return -1;
    }
    
    size_t i;
    for (i = 0; i < n; ++i) {
        fprintf(ofp, "%zu\n", assignments[i]);
    }
    
    fclose(ofp);
    return 0;
    
}

int write_centers_to_file(char *fileName, double *centers, size_t d, size_t k) {
    FILE *ofp;
    
    ofp = fopen(fileName, "w");
    if (ofp == NULL) {
        printf("Error opening centers output file!");
        return -1;
    }
    
    // Write cluster count to file's first row.
    fprintf(ofp, "%zu\n", k);
   
    size_t i, j;
    // For each cluster center.
    for (i = 0; i < k; ++i) {
        // For each dimension of the cluster senter
        for (j = 0; j < d; ++j) {
            fprintf(ofp, "%.17g ", centers[i*d+j]);
        }
        fprintf(ofp, "\n");
    }
    
    fclose(ofp);
    return 0;
    
}

double* get_centers_from_file(char *fileName, size_t d, size_t *k) {
    
    
    FILE *ifp;
    ifp = fopen(fileName, "r");
    
    if (ifp == NULL) {
        return 0;
    }
    // In our file format, cluster count is on the first row.
   fscanf(ifp, "%zu\n", k);
    
    // According to the first row, we allocate memory.
    double* clusterCenters = malloc(sizeof(double) * (*k) * d);
    ALCHECK(clusterCenters, "clusterCenters");
    
    size_t i, j;
    for (i = 0; i < *k; ++i) {
        for (j = 0; j < d; ++j) {
            if (feof(ifp)) {
                printf("Incorrect centers file! Not enough centers in file.\n");
                exit(1);
            }
            fscanf(ifp, "%lf", &clusterCenters[i*d+j]);
        }
    }
 
    /*
    knr = fscanf(ifp, "%zu", &knr);
    if (!feof(ifp)) {
        printf("More data in file. Incorrect centers file!\n");
        exit(1);
    }
    */ 
     
    fclose(ifp);
    
    return clusterCenters;
    
}

// Source: http://stackoverflow.com/questions/28115724/getting-big-random-numbers-in-c-c (date: 27.04.2016)
unsigned long long llrand() {
    unsigned long long r = 0;

    for (int i = 0; i < 5; ++i) {
        r = (r << 15) | (rand() & 0x7FFF);
    }

    return r & 0xFFFFFFFFFFFFFFFFULL;
}


void find_closest_centers(double *points, size_t n, size_t d, size_t k, double* clusterCenters, size_t* assignments, metricType metric) {
    
    // For each point find closest center.
    size_t i, centerNr;
    double distance;
    for (i = 0; i < n; ++i) {
        double lowest_distance = DBL_MAX;
        for (centerNr = 0; centerNr < k; ++centerNr) {
            distance = metric(d, points + i*d, clusterCenters + centerNr * d);
            if (distance < lowest_distance) {
                assignments[i] = centerNr;
                lowest_distance = distance;
            }
        }
    }
}

int algorithm_choice(char* algorithm) {
    if (strcmp(algorithm, "1") == 0 || strcmp(algorithm, "lloyd") == 0) {
        return LLOYD;
    }
    else if (strcmp(algorithm, "2") == 0 || strcmp(algorithm, "elkan") == 0) {
        return ELKAN;
    }
    else if (strcmp(algorithm, "3") == 0 || strcmp(algorithm, "macqueen") == 0) {
        return MCQUEEN;
    }
    else if (strcmp(algorithm, "4") == 0 || strcmp(algorithm, "hamerly") == 0) {
        return HAMERLY;
    } 
    else if (strcmp(algorithm, "5") == 0 || strcmp(algorithm, "closest") == 0) {
        return FIND_CLOSEST;
    }
    else if (strcmp(algorithm, "6") == 0 || strcmp(algorithm, "hartigan") == 0) {
        return HARTIGAN;
    }
    return -1;
}

int metric_choice(char* metric) {
    if (strcmp(metric, "1") == 0 || strcmp(metric, "euclidean") == 0 || mystrcmp(metric, "euclidean")) {
        return EUCLIDEAN_METRIC;
    } else if (strcmp(metric, "2") == 0 || strcmp(metric, "manhattan") == 0 || mystrcmp(metric, "euclidean")) {
        return MANHATTAN_METRIC;
    } 
    return -1;
}


int initmethod_choice(char* method) {
    if (strcmp(method, "1") == 0 || strcmp(method, "forgy") == 0) {
        return RANDOM_INIT;
    }
    else if (strcmp(method, "2") == 0 || strcmp(method, "kpp") == 0) {
        return KPP_INIT;
    }
    else if (strcmp(method, "3") == 0 || strcmp(method, "firstn") == 0) {
        return FIRST_N;
    }
    else if (strcmp(method, "4") == 0 || strcmp(method, "furthest") == 0) {
        return FURTHEST_FIRST;
    }
    else if (strcmp(method, "5") == 0 || strcmp(method, "partition") == 0) {
        return RANDOM_PARTITION;
    }    
    return -1;
    
}


double calculate_e_value(double *points, size_t n, size_t d, double *clusterCenters, size_t *assignments, metricType metric) {
        size_t i;
        
        double eval = 0, dist;
        for (i = 0; i < n; ++i) {
                dist = metric(d, points + i*d, clusterCenters + assignments[i] * d);
                eval += dist;
        }
    return eval;
}

double calculate_sqrd_value(double *points, size_t n, size_t d, double *clusterCenters, size_t *assignments, metricType metric) {
        size_t i;
        
        double eval = 0, dist;
        for (i = 0; i < n; ++i) {
                dist = metric(d, points + i*d, clusterCenters + assignments[i]*d);
                dist *= dist;
                eval += dist;
        }
    return eval;
}

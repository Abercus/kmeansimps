#include <float.h>
#include <math.h>
#include "metrics.h"
#define PI 3.14159265358979323846

double euclid(size_t d, double* data1, double* data2) {
    
    size_t i;
    double distance = 0;
    
    for (i = 0; i < d; ++i) {
        distance += pow(data1[i] - data2[i], 2);
    }
    
    return sqrt(distance);
}


 double manhattan(size_t d, double* data1, double* data2) {
     
     size_t i;
     double distance = 0;
     
     for (i = 0; i < d; ++i) {
        distance += fabs(data1[i] - data2[i]);
     }
     
     return distance;
 }
 
double angular(size_t d, double* data1, double* data2) {
    
    double dot_p = 0;
    size_t i;
    // Calculate dot product
    for (i = 0; i < d; ++i) {
        dot_p += data1[i]*data2[i];
    }
    // Return 1 - dot_p, because we need it as a distance measure.
 
    return 1 - dot_p;
    //return sqrt(2-2*cos(dot_p));
    //return (acos(dot_p)/PI);
    
    
}

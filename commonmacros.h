#ifndef COMMONMACROS_H
#define COMMONMACROS_H
#include <stddef.h>

typedef double (*metricType)(size_t, double*, double*);

// Min and max macros
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

#define DEBUG 1

// Macro to check if managed to allocate memory.
#define ALCHECK(name, s) if (!name) { \
                            fprintf(stderr, "Failed to allocate memory for %s\n", s); \
                            exit(1); \
}


// Macro to free memory
#define ALFREE(name) if (name) {free(name);}


#endif // COMMONMACROS_H
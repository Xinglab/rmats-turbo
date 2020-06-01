#include "stdio.h"
#include "string.h"
#include "../include/rMATS_wrapper.h"


// TODO len of output.
void r_likelihood_test(char** inputf, double* cutoff, int* nthread, double *output) {
    int len;
    char* filename = *inputf;
    printf("%g, %g\n", output[0], output[1]);
    likelihood_test_wrapper(filename, *cutoff, *nthread, output, &len);

    return;
}

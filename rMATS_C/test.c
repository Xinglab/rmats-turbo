#include "time.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_vector.h"


int main(int argc, char *argv[])
{
    double dur;
    clock_t start,end;
    int times = 300, i = 0;
    size_t len = 200;
    gsl_vector *vec1 = gsl_vector_alloc(len), *vec2 = gsl_vector_alloc(len), *vec3 = NULL;

    start = clock();
    for (i = 0; i < times; ++i) {
        vec3 = gsl_vector_alloc(vec1->size);
        gsl_vector_memcpy(vec3, vec1);
        gsl_vector_add(vec3, vec2);
        gsl_vector_free(vec3);
    }
    end = clock();
    dur = (double)(end - start);
    printf("using vec3: %f\n", dur/CLOCKS_PER_SEC);

    start = clock();
    for (i = 0; i < times; ++i) {
        gsl_vector_add(vec1, vec2);
        gsl_vector_sub(vec1, vec2);
    }
    end = clock();
    dur = (double)(end - start);
    printf("using vec1: %f\n", dur/CLOCKS_PER_SEC);

    start = clock();
    for (i = 0; i < 100; ++i) {
        gsl_cdf_gaussian_Pinv(0.5, 1);
    }
    end = clock();
    dur = (double)(end - start);
    printf("time for ppf: %f\n", dur/CLOCKS_PER_SEC);

    printf("sizeof double: %ld\n", sizeof(double));

    return 0;
}
